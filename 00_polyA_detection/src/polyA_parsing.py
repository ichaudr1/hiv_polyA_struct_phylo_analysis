'''
@author ichaudr

The goal of this script is to parse out the putative polyA hairpin sequences for a set of HIV-1 genome sequences.

Workflow:
    1. Parse TSL4 from the full genome: TSL4 is defined as the start of the sequence to the gag start codon.
    2. Find the AAUAAA signal sequence and pull a segment that is ~100 nt upstream and downstream of the signal.
    3. Iteratively clip off bases from the 5' and 3' end and predict the secondary structure following each clipping. 
    4. Restrain resulting structures.
    5. Determine the representative strains for each subtype. 
    6. Make an alignment of the polyA hairpins from the representative sequences. 
    7. Write the aligned polyA hairpin sequences and the corresponding vstring to json 
'''

from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.Applications import PhymlCommandline

from tqdm import tqdm

#from fetch_pol_gene import pol
#from target_strain import TargetStrain
from sim_filter import get_rep_strains

from matplotlib import pyplot as plt
import numpy as np
#matplotlib.use('Qt5Agg')

import os
import sys
import json
import uuid
import statistics



##############
# Parameters #
##############

#--> Entrez authentication
Entrez.email = 'ichaudr1@umbc.edu'
Entrez.api_key = 'cd2f5bf7a67d086647ec33da2c985e018d08'

#--> Paths
### Path to full genomes
full_genomes_path = ''

### Vienna RNAFold executable path
rna_fold_path = ''

#--> Parameters for the minimization of polyA segment (see function below)

#Filter out segments shorter than min_len (30 by default until read in from input)
min_len = 30

#Filter out segments longer than max_len (60 by default until read in from input)
max_len = 60

#Filter out segments with a ratio of unpaired residues greater than max_unbp_ratio (0.45 by default until read in from input)
max_unbp_ratio = 0.45

#Filter out segments where the position of the AAUAAA is greater than hwd_limit from the middle of the segment (0.3 by default until read in from input)
hwd_limit = 0.3

#--> Parameters for ClustalW alignment and Tree formation
clustal_path = ''

#--> Output

#Generate a unique id for the output folder
folder = '../output/'
output_id = 'polyA_parse_' + str(uuid.uuid1())


def fold(seq):
    '''
    Folds the RNA sequence that is passed in and returns the vienna string.

    Parameters
    ----------
    seq: str
        RNA sequence
    
    Returns
    -------
    (vienna_str, free_energy): (str, float)
        The vienna string of the predicted structure and its free energy.
    '''

    #vrna_path='/Users/ichaudr/Documents/Lab/misc/viennaRNA/bin/RNAfold'
    #vrna_path='/data/DataRAID/mfslab/users/issac/laptopsync/Lab/misc/viennaRNA_2/ViennaRNA-2.4.16/src/bin/RNAfold'
    outfile='out_' + str(output_id) + '.temp'

    command = 'echo {sequence} | {path} --noPS --outfile={outfile}'

    os.system(command.format(sequence=seq, path=rna_fold_path, outfile=outfile))

    to_return = None
    with open(outfile, 'r+') as file:
        lines = [l.strip() for l in file.readlines() if len(l.strip()) != 0]
        vienna_str = lines[-1].split(' ')[0]
        free_energy = float(lines[-1].split(' ')[-1].replace('(', '').replace(')', ''))
        to_return = (vienna_str, free_energy)
        file.truncate(0)
        
    return to_return

def polyA_hp_seq(seq, ext_range):
    '''
    Finds the most likely sequence for the polyA HP by finding the hexameric signal (AAUAAA) and going up/downstream a certain range.

    Parameters
    ----------
    seq: str
        The full leader sequence
    ext_range: int
        The distance to pull upstream and downstream.
    
    Returns
    -------
    polyA_hp_seq: str
        The sequence of the polyA hairpin
    '''

    polyA_index = seq.index('AAUAAA')
    return seq[max(polyA_index - ext_range, 0) : min(polyA_index + 5 + ext_range, len(seq)-1)]

#Not used
def clip(rna_seq, v_str):
    '''
    Clips off residues from the 5' and 3' ends that are not predicted to be in a hairpin.

    Parameters
    ----------
    rna_seq: str
        The RNA sequence.
    v_str: str
        The vienna string for the predicted secondary strucuture.
    
    Return
    -------
    (clip_rna_seq, clip_v_str): (str, str)
        The clipped RNA sequence and Vienna String. 
    '''

    if len(rna_seq) != len(v_str):
        print('RNA and Vienna String are not the same length.')
        return
    
    if len(v_str.replace('.', '')) == 0:
        print('All unstructured, no hairpin present to clip.')
        return

    clip5 = 0
    clip3 = 0

    #Find number of bases to clip off 5' end
    i = 0
    while v_str[i] == '.':
        clip5 += 1
        i += 1
    
    #Find number of bases to clip off 3' end
    i = len(v_str) - 1
    while v_str[i] == '.':
        clip3 += 1
        i -= 1
    
    return (rna_seq[clip5:len(rna_seq) - clip3],v_str[clip5:len(v_str) - clip3] )

def single_hp(vstr):
    '''
    Determines if there is only a single hp in a vienna string. 

    Parameters
    ----------
    vstr: str
        The vienna string
    
    Returns
    -------
        True if there is a single hairpin, false otherwise. 
    '''

    num_hp = 0

    opp_brak = {'(':')', ')':'('}
    curr_brak = '('

    for v in vstr:
        if v not in opp_brak.keys():
            continue
        if v == opp_brak[curr_brak]:
            num_hp += 1
            curr_brak = opp_brak[curr_brak]

    return num_hp == 1

def minimize_structure(seq, min_len, max_len, max_unbp_ratio, hwd_limit):
    '''
    Aims to determine the lowest energy hairpin structure from the given sequence. Works by iteratively removing bases from the 5' and 3' ends, and keeping track of the free energy of 
    the predicted secondary structure after each clipping. 

    Parameters
    ----------
    seq: str
        The sequence of interest.
    min_len: int
        The shortest length of the final sequence that is acceptable. 
    max_unbp_ratio: float
        The maximum ratio of unbase paired residues allowed in the minimized segment.
    hwd_limit: float
        The half-way deviation limit (hwd) is how far the beginning of the AAUAAA can be from the center of the segment.
    
    Returns
    -------
    min_seq: str
        The sequence of the structure with the lowest free energy. 
    '''

    tqdm.write('-'*10 + 'Minimizing')

    #The range of clippings from the 5' end
    five_clips = list(range(0, seq.index('AAUAAA')))

    #The range of clippings from the 3' end
    three_clips = list(range(0, len(seq) - seq.index('AAUAAA') - 6))

    #Make all the combinations of 5' and 3' clippings
    all_combos = []

    for c5 in five_clips:
        for c3 in three_clips:
            all_combos.append((c5, c3))
    
    #Filter out any of the combos that would result in a sequence that is smaller than the min_len
    filtered_combos = []
    for c in all_combos:

        if seq.index('AAUAAA') > len(seq) - c[1] or seq.index('AAUAAA') < c[0]:
            continue
        if len(seq) - sum(c) < min_len or len(seq) - sum(c) > max_len:
            continue

        filtered_combos.append(c)
    
    tqdm.write('Testing ' + str(len(filtered_combos)) + ' combos.')

    #Dictionary to hold the free energy for each of the clipping combos
    combo_fe = {}

    for c in filtered_combos:

        sub_seq = seq[c[0] : len(seq) - c[1]]
        fold_results = fold(sub_seq)

        if not single_hp(fold_results[0]) or fold_results[0][0] == '.' or fold_results[0][-1] == '.':
            continue
        
        '''
        if 'AAUAAA' not in sub_seq:
            continue

        if len(sub_seq) < min_len or len(sub_seq) > max_len:
            continue'''

        tqdm.write('-' * 15 + 'Testing ' + str(c))
        tqdm.write((' ' * 25) + sub_seq)
        tqdm.write((' ' * 25) + fold_results[0] + '\t' + str(fold_results[1]))

        if fold_results[0].count('.')/len(fold_results[0]) > max_unbp_ratio:
            continue

        if abs(sub_seq.index('AAUAAA')/len(sub_seq) - .5) > hwd_limit:
            continue
        
        combo_fe[c] = float(fold_results[1])
        tqdm.write('-' * 15 + 'Added ' + str(c))


    if len(combo_fe.keys()) == 0:
        return None

    #Return the sequence with the lowest free energy
    lowest_combo = list(combo_fe.keys())[0]

    for c in combo_fe.keys():
        if combo_fe[c] <= combo_fe[lowest_combo]:
            lowest_combo = c

    min_seq = seq[lowest_combo[0]:len(seq) - lowest_combo[1]]
    tqdm.write('-' * 15 + 'Returning' + str(lowest_combo))
    return min_seq

def gag_position(accession, cache='../gag_pos_cache/'):
    '''
    Will determine the gag protein start position. This is used to parse out the initial fragment used for the polyA id pipeline.

    Parameters
    ----------
    accesssion: str
        The nucleotide accession for the strain of interest.
    cache: str
        Directory to store all collected gag positions.
    
    Returns
    --------
    gag_pos: int
        The gag start position. Returns -1 if no gag annotation is found in the gb file. 
    '''
    
    #Check if the gag position for the accession has already been stored in cache.
    for file in os.listdir(cache):
        if file == accession:
            with open(cache + file, 'r') as data_file:
                gag_pos = int(data_file.readline())
                return gag_pos

    genebank_data = None
    try:
        genebank_data = Entrez.read(Entrez.efetch(db='nuccore', id=accession, rettype='gbwithparts', retmode='xml'))[0]
    except:
        print('Genebank data for ' + accession + ' could not be downloaded.')
        return -1
        
    if genebank_data == None:
        print('Genebank data for ' + accession + ' could not be downloaded.')
        return -1
    
    gag_interval = None
    try:

        for feat in genebank_data['GBSeq_feature-table']:
            if feat['GBFeature_key'] == 'CDS':
                is_gag = False
                for quality in feat['GBFeature_quals']:
                    if quality['GBQualifier_name'] == 'gene':
                        if 'gag' in str(quality['GBQualifier_value']).lower().strip():
                            is_gag = True
                    if quality['GBQualifier_name'] == 'note':
                        if 'gag' in str(quality['GBQualifier_value']).lower().strip():
                            is_gag = True
                    if quality['GBQualifier_name'] == 'product':
                        if 'gag' in str(quality['GBQualifier_value']).lower().strip():
                            is_gag = True
                    
                if is_gag:
                    gag_interval = feat['GBFeature_intervals']

    except:
        print('Genebank data for ' + accession + ' could not be parsed')
        gag_interval = None
    

    if gag_interval == None:
        print('Gag interval was none.')
        return -1
    else:
        gag_pos = gag_interval[0]['GBInterval_from']
    
    with open(cache + accession, 'w+') as data_file:
        data_file.write(gag_pos)
    
    return int(gag_pos)

def accession(record_id):
    '''
    Returns the genome accession given the record id. 

    Parameters
    ----------
    record_id: str
        In the format: SUBTYPE.COUNTRY.SAMPLING_YR.NAME.ACCESSION.PATIENT_CODE
    '''

    return str(record_id.split('.')[-2])

##########################################################################################
#Parse the input file and set all parameters
##########################################################################################

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Invalid arguments. Usage: polyA_parsing.py [input_file.json]')
        quit()

    input_file = sys.argv[1]
    input_file = '../input/' + input_file
    inputs = json.load(open(input_file, 'r'))
    print('Reading input from ', str(input_file), '...')

    #--> Paths
    full_genomes_path = inputs['paths']['full_genomes_path']
    rna_fold_path = inputs['paths']['rna_fold_path']
    clustal_path = inputs['paths']['clustal_path']

    #--> Parameters for the pulling initial polyA segment
    ext_range = inputs['initial_seg']['ext_range']

    #--> Parameters for the minimization of polyA segment
    min_len = inputs['minimization']['min_len']
    max_len = inputs['minimization']['max_len']
    max_unbp_ratio = inputs['minimization']['max_unbp_ratio']
    hwd_limit = inputs['minimization']['hwd_limit']
    
    #--> Output
    folder = inputs['output']['folder']
    os.mkdir(folder + output_id)
    output_root = folder + output_id + '/{item}'

    #A dictionary holding peices of info that will be logged at the end. 
    STATS = {'no_ltr':0, 'no_pa_signal':0, 'no_pa_discovered':0, 'starting_seqs':0}
    ##########################################################################################

    print('Logging parameters...')
    json.dump(inputs, open(output_root.format(item='param.json'), 'w+'))

    ##########################################################################################
    
    #Load all genomes from JSON file
    whole_genomes = json.load(open(full_genomes_path, 'r'))
    STATS['starting_seqs'] = len(whole_genomes)

    #Pull out the 5' LTR for each genome (i.e position 0 to start of GAG) and determine all instances of the AAUAAA polyA signal
    all_five_ltrs = []
    for genome in tqdm(whole_genomes, desc='Getting 5\' LTRs'):
        gag_pos = gag_position(genome['accession'])
         
        if gag_pos < min_len:
            STATS['no_ltr'] += 1
            tqdm.write('Excluding: ' + genome['subtype'] + ':' + genome['accession'] + ' - GAG position invalid.')
            continue

        five_ltr = str(genome['sequence'][0:gag_pos]).upper().replace('T', 'U')
        polyA_signal_positions =  [i for i in range(len(five_ltr)) if five_ltr.startswith('AAUAAA', i) and i < len(five_ltr) - ext_range + 1]
        
        if len(polyA_signal_positions) < 1:
            STATS['no_pa_signal'] += 1
            tqdm.write('Excluding: ' + genome['subtype'] + ':' + genome['accession'] + ' - polyA signal invalid')
            continue

        all_five_ltrs.append(genome)
        all_five_ltrs[-1]['full_genome_sequence'] = all_five_ltrs[-1]['sequence']
        all_five_ltrs[-1]['sequence'] = five_ltr
        all_five_ltrs[-1]['polyA_signal_positions'] = polyA_signal_positions

    #Attempt to parse out polyA hairpin for each leader
    polyA_cache = '../discovered_polyA_cache/'
    for ltr in tqdm(all_five_ltrs):
        tqdm.write('***' + ltr['accession'] + '***')
        ltr["pa_ext_seq"] = None
        ltr["pa_seq"] = None
        ltr["pa_vstr"] = None
        ltr["pa_free_energy"] = None
        ltr['representative'] = 0

        tqdm.write('\t Checking cache...')

        #Check if the polyA hp for the accession has already been stored in cache.
        
        if ltr['accession'] in os.listdir(polyA_cache):
            with open(polyA_cache + ltr['accession'], 'r') as data_file:
                ltr["pa_ext_seq"] = str(data_file.readline()).strip()
                ltr["pa_seq"] = str(data_file.readline()).strip()
                ltr["pa_vstr"] = str(data_file.readline()).strip()
                ltr["pa_free_energy"] = float(data_file.readline())
            continue

        discovered_pas = []

        for pos in ltr['polyA_signal_positions']:
            tqdm.write('Testing position ' + str(pos))
            
            pa_ext_seq = ltr['sequence'][max(pos - ext_range, 0) : min(pos + ext_range, len(ltr['sequence']) - 1)]
            if len(pa_ext_seq) < min_len or 'AAUAAA' not in pa_ext_seq:
                tqdm.write('-'*10 + 'Failed - no polyA sequence found')
                continue

            vstr, free_energy = fold(seq=pa_ext_seq)
            tqdm.write(pa_ext_seq)
            tqdm.write(vstr + '\t' + str(free_energy))
            tqdm.write('\t\t|')
            tqdm.write('\t\t|')
            tqdm.write('\t\tv')

            #Minimize the segment    
            pa_seq = minimize_structure(seq=pa_ext_seq, min_len=min_len, max_len=max_len, max_unbp_ratio=max_unbp_ratio, hwd_limit=hwd_limit)

             #Get the fold for the minimized segment
            if pa_seq != None:
                vstr, free_energy = fold(seq=pa_seq)
                discovered_pas.append([pa_ext_seq, pa_seq, vstr, free_energy])
            else:
                continue
        
        if len(discovered_pas) == 0:
            STATS['no_pa_discovered'] += 1
            print('-'*10 + 'Failed - no polyA sequence found after minimization of all instances')
            print('***' * 40 + '\n')
            continue
    
        discovered_pas = sorted(discovered_pas, key=lambda x: x[3])

        ltr["pa_ext_seq"] = discovered_pas[0][0]
        ltr["pa_seq"] = discovered_pas[0][1]
        ltr["pa_vstr"] = discovered_pas[0][2]
        ltr["pa_free_energy"] = discovered_pas[0][3]

        with open(polyA_cache + ltr['accession'], 'w+') as data_file:
            data_file.write(str(discovered_pas[0][0])+'\n')
            data_file.write(str(discovered_pas[0][1])+'\n')
            data_file.write(str(discovered_pas[0][2])+'\n')
            data_file.write(str(discovered_pas[0][3])+'\n')

        #0 by default. Will be 1 if determined to be representative of subtype
        ltr['representative'] = 0
    
    #Check if all required strains have a polyA discovered
    strains_force_include = []
    for a in tqdm(inputs['strains_force_include'], desc='Checking required strains'):
        if a not in [l['accession'] for l in all_five_ltrs if l['pa_seq'] != None]:
            tqdm.write(str('*** THE REQUIRED STRAIN HAS NO DISOCOVERED POLYA - IT WILL BE IGNORED: ' + a))
            continue
        strains_force_include.append(a)


    #Determine representative strains and denote as such - only considering strains where a polyA was detectable. 
    rep_strains = get_rep_strains(input_ltrs=[l for l in all_five_ltrs if l['pa_seq'] != None])

    for ltr in all_five_ltrs:
        if ltr["accession"] in strains_force_include:
            ltr['representative'] = 1
            continue
        if ltr['subtype'] not in rep_strains.keys():
            continue
        ltr['representative'] = int(ltr['accession'] in [l['accession'] for l in rep_strains[ltr['subtype']]])

    json.dump(all_five_ltrs, open(output_root.format(item='raw_output.json'), 'w+'))
    
    os.system('rm out_' + str(output_id) + '.temp')

    #Alignment of the representative polyA
    ### |-> write all the sequences upstream of AAUAAA to a fasta file in the output directory
    ### |-> write all the sequences downstream of AAUAAA to a fasta file in the output directory
    ### |-> align the upstream and downstream sequences seperately 
    ### |-> combine the aligned sequnces
    ### |-> write out the json with the aligned represenative polyAs


    ### |-> write all the sequences upstream of AAUAAA to a fasta file in the output directory
    print('Writing sequences upstream of AAUAAA ...')
    to_write = [SeqRecord(seq=Seq(ltr['pa_seq'][0:ltr['pa_seq'].index('AAUAAA')]), id=(ltr['subtype'] + '-' + ltr['accession'])) for ltr in all_five_ltrs if ltr['representative']]
    with open(output_root.format(item='upstream_aauaaa_seqs.fasta'), 'w+') as file:
        SeqIO.write(to_write, file, 'fasta')
    
    ### |-> write all the sequences downstream of AAUAAA to a fasta file in the output directory
    print('Writing sequences downstream of AAUAAA ...')
    to_write = [SeqRecord(seq=Seq(ltr['pa_seq'][ltr['pa_seq'].index('AAUAAA') + 6:]), id=(ltr['subtype'] + '-' + ltr['accession'])) for ltr in all_five_ltrs if ltr['representative']]
    with open(output_root.format(item='downstream_aauaaa_seqs.fasta'), 'w+') as file:
        SeqIO.write(to_write, file, 'fasta')
    
    ### |-> align the upstream and downstream sequences seperately
    print('Aligning sequences upstream of AAUAAA ...')
    clustal_cmd = ClustalwCommandline(cmd=clustal_path, infile=output_root.format(item='upstream_aauaaa_seqs.fasta'), outfile=output_root.format(item='aligned_upstream_aauaaa_seqs.fasta'), type='DNA', output='FASTA')
    clustal_cmd()

    ### |-> align the upstream and downstream sequences seperately
    print('Aligning sequences downstream of AAUAAA ...')
    clustal_cmd = ClustalwCommandline(cmd=clustal_path, infile=output_root.format(item='downstream_aauaaa_seqs.fasta'), outfile=output_root.format(item='aligned_downstream_aauaaa_seqs.fasta'), type='DNA', output='FASTA')
    clustal_cmd()

    ### |-> combine the aligned sequnces and outputting to fasta
    print('Combining aligned upstream and downstream polyA hairpin sequences...')
    to_write = []
    for up_record in SeqIO.parse(open(output_root.format(item='aligned_upstream_aauaaa_seqs.fasta'),'r'),'fasta'):
        up_seq = str(up_record.seq)
        down_seq = ''
        for down_record in SeqIO.parse(open(output_root.format(item='aligned_downstream_aauaaa_seqs.fasta'),'r'),'fasta'):
            if up_record.id != down_record.id:
                continue
            down_seq = str(down_record.seq)
        
        if down_seq == '':
            print('Error while combining: ' + up_record.id)
        
        complete_aligned_seq = Seq(up_seq + 'AAUAAA' + down_seq)
        to_write.append(SeqRecord(complete_aligned_seq, id=up_record.id.split('-')[1], description=up_record.id))

    with open(output_root.format(item='aligned_pa_seqs.fasta'), 'w+') as file:
        SeqIO.write(to_write, file, 'fasta')
    
    #Writing all the aligned sequences to a JSON file so that information about the vienna string and free energy can be noted. 
    print('JSON FILE: Combining aligned upstream and downstream polyA hairpin sequences...')
    to_write_json = []
    for record in to_write:
        to_append = {'accession': record.id, 'subtype': record.description.split('-')[0], 'aligned_seq': str(record.seq), 'pa_vstr': None, 'pa_free_energy':0}

        for ltr in all_five_ltrs:
            if ltr['accession'] == record.id:
                to_append['pa_vstr'] = ltr['pa_vstr']
                to_append['pa_free_energy'] = ltr['pa_free_energy']

        to_write_json.append(to_append)
    
    with open(output_root.format(item='aligned_pa_seqs.json'), 'w+') as file:
        json.dump(to_write_json, file)
    
    #Write out all strains with identified polyAs to a fasta file to be used in locarna analysis
    records_to_write = [SeqRecord(Seq(s['pa_seq']), id=s['subtype'] + '-' + s['accession']) for s in all_five_ltrs if s['pa_seq'] != None]
    SeqIO.write(records_to_write, open(output_root.format(item='all_pas.fasta'), 'w+'), 'fasta')

    #Write the stats file
    with open(output_root.format(item='stats.txt'), 'w+') as file:
        file.write(str(STATS))




    

    






        

