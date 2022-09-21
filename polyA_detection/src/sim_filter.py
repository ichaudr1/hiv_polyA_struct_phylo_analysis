'''
@author ichaudr
August 2022

Objective: The main function of this file is to filter out the representative strains for each subtype. 

'''

from statistics import mean
from Bio import SeqIO
from Bio import pairwise2
from Bio import Entrez
from tqdm import tqdm

Entrez.email = 'ichaudr1@umbc.edu'
Entrez.api_key = 'cd2f5bf7a67d086647ec33da2c985e018d08'

def get_percent_matches(seq1, seq2, scoring={'match':2, 'mismatch':-1, 'gap_open':-1, 'gap_extd': -0.5}):
    '''
    Returns the percent matches between two sequences. 
    Parameters
    ----------
    seq1, seq2: str
        The two sequences to calculate the percent matches for.
    scoring: dictionary
        The scoring to use for the alignment in the following format: scoring={'match':__, 'mismatch':__, 'gap_open':__, 'gap_extd': __}
        Default: scoring={'match':2, 'mismatch':-1, 'gap_open':-1, 'gap_extd': -0.5}
    
    Returns
    -------
    average_percent_matches: float
        The average percent match between seq1 and seq2. Calculated as: #matches / (#matches + #mismataches); gaps are not included
    '''

    #A list of percent matches of all the alignments - used to calculate the average 
    percent_matches = []

    for a in pairwise2.align.globalms(seq1, seq2, scoring['match'], scoring['mismatch'], scoring['gap_open'], scoring['gap_extd']):
        #The format of an alignment:
            # [0] : first sequence aligned
            # [1] : second sequenced aligned
            # [2] : alignment score
            # [3] : starting position of the alignment
            # [4] : ending position of the alignment 
    
        seq1_aligned = a[0]
        seq2_aligned = a[1]
        alignment_length = a[4]
        num_gaps = 0
        num_matches = 0

        for n in range(alignment_length):

            #Adjust the gap count
            if seq1_aligned[n] == '-' or seq2_aligned[n] == '-':
                num_gaps += 1
                continue

            #Adjust match count
            if seq1_aligned[n] == seq2_aligned[n]:
                num_matches += 1
        
        #Calculate the percent matches for the current alignment and append it to the list of percent matches
        percent_match = num_matches / (alignment_length - num_gaps)
        percent_matches.append(percent_match)
        
    #Return the average percent matches
    if len(percent_matches) <= 0:
        return 0

    average_percent_matches = mean(percent_matches)

    return average_percent_matches

def get_rep_strains(input_ltrs, target_subtypes=['A', 'B', 'C', 'D', 'F1', 'F2', '01_AE', '02_AG', 'G', 'H', 'J', 'K', 'L', 'P'], target_count_per_subtype=5):
    '''
    Identifies a set of representative strains for each subtype.

    Parameters
    ----------
    input_ltrs: [{ltr_dictionary}}]
        The set of ltrs to pull the representative strains from
    target_subtypes: [str]
        The target subtypes
    target_count_per_subtype: int
        The max number of representative strains per subtype to return
    
    Returns
    -------
    rep_strains: {subtype: [target_strains]}
        A dictionary that holds the representative strains per subtype
    '''

    target_subtypes=['A1', 'A2', 'A3', 'A4', 'A6', 'A7', 'A8', 'B', 'C', 'D', 'F1', 'F2', 'G', 'H', 'J', 'K', 'L', '01_AE', '02_AG', 'P', 'N', 'O']

    #Sort the input target strains based of subtypes
    ssts = {} #ssts stands for subtype sorted target strains
    for subtype in target_subtypes:
        ssts[subtype] = [s for s in input_ltrs if s['subtype'] == subtype]
    
    #Get the representative strains
    rep_strains = {}
    for subtype in tqdm(ssts.keys(), desc=('Getting representative strains for')):

        if len(ssts[subtype]) == 0:
            continue
        
        if len(ssts[subtype]) <= target_count_per_subtype:
            rep_strains[subtype] = ssts[subtype]
            continue
        
        #For each subtype, calculate the average % similarity of each strain to the rest of the strains in the subtype
        rep_strains[subtype] = []
        current_subtype_strains = ssts[subtype]
        avg_pcent_sims_for_each_strain = {}
        for i in range(len(current_subtype_strains)):
            st_i_percent_sims = []
            for j in range(len(current_subtype_strains)):
                if i == j:
                    continue
                st_i_percent_sims.append(get_percent_matches(current_subtype_strains[i]['pa_seq'], current_subtype_strains[j]['pa_seq']))
            avg_pcent_sims_for_each_strain[i] = mean(st_i_percent_sims)
        
        avg_pcent_sims_for_each_strain = dict(sorted(avg_pcent_sims_for_each_strain.items(), key=lambda item: item[1]))
        for num in range(target_count_per_subtype):
            rep_strains[subtype].append(current_subtype_strains[list(avg_pcent_sims_for_each_strain.keys())[-1 - num]])

    return rep_strains

