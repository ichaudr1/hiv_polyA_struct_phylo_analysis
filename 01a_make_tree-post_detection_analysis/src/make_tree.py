'''
@author ichaudr

This script takes a set of sequences of a specific hairpin structure from a set of strains and constructs a phylogeny for the sequences that are denoted to come from strains that are
representative of their subtype. 

1. Pull the polymerase gene sequence to construct the phylogeny. 
2. Draw the tree. 
3. Annotate the tree.

'''

from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.Applications import PhymlCommandline
from tqdm import tqdm

from ete3 import Tree, faces, TreeStyle, TextFace
from ete3.treeview.main import NodeStyle
import svgwrite

import os
import sys
import json
import uuid

##############
# Parameters #
##############

#--> Entrez authentication
Entrez.email = ''
Entrez.api_key = ''

#--> input JSON formatted strain data
strain_file = ''

#--> Parameters for ClustalW alignment and Tree formation
clustal_path = ''
phyml_path = ''
bootstrap = 1
phylo_target = 'env'


#--> Output

#Generate a unique id for the output folder
folder = '../output/'
output_id = 'make_tree_' + str(uuid.uuid1())
os.mkdir(folder + output_id)
output_root = folder + output_id + '/{item}'


def fetch_phylo(accession, target='env'):
    '''
    Will download the genebank record for the accession and return the sequence for the target gene.

    Parameters
    ----------
    accession: str
        The genome of interest
    target: str
        The name of the gene that should be retrieved. Must be CDS target. 
    
    Returns
    -------
    target_seq: str
        The sequence found in the target genome for the target gene. None if not found.
    '''

    std_aa = ['G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T']
    target_seq = None
    try:
    
        results = Entrez.read(Entrez.efetch(db='nuccore', id=accession, rettype='gbwithparts', retmode='xml'))
    
    except:
        print('None for ', accession + '. There was an Entrez error.')
        return None

    for feat in results[0]['GBSeq_feature-table']:
        if feat['GBFeature_key'] == 'CDS':
            target_found = False
            for quality in feat['GBFeature_quals']:
                if quality['GBQualifier_name'] == 'gene':
                    if target in str(quality['GBQualifier_value']).lower().strip():
                        target_found = True
                if quality['GBQualifier_name'] == 'note':
                    if target in str(quality['GBQualifier_value']).lower().strip():
                        target_found = True
                if quality['GBQualifier_name'] == 'product':
                    if target in str(quality['GBQualifier_value']).lower().strip():
                        target_found = True
                if not target_found:
                    continue
                if quality['GBQualifier_name'] == 'translation':
                    target_seq = str(quality['GBQualifier_value']).upper()
    
    if target_seq == None:
        return None
    
    to_return = ''
    for a in target_seq:
        if a in std_aa:
            to_return = to_return + a
    return to_return

def get_target_seqs(targets, output_file=output_root.format(item='phylo_target_seqs.fasta')):
    '''
    Saves the polymerase sequences for all the accessions that are passed in to 1 FASTA file. This is used for the phylogenetic tree. 

    Parameters
    ----------
    targets: [(accession, id)]
        A list of tuples where the first element is the target accession and the second element is the target id to be saved in the fasta
    
    Returns
    -------
    None - output is written to the output directory
    '''
    
    to_write = []

    for t in tqdm(targets, desc='Getting polymerase seqs'):
        phylo_seq = fetch_phylo(t[0])
        if phylo_seq == None:
            print('No polymerase found for ', t[1])
            continue
        
        to_write.append(SeqRecord(Seq(phylo_seq), id=t[0]))
    
    with open(output_file, 'w+') as file:
        SeqIO.write(to_write, file, 'fasta')

def make_face(strain_id, pa_seq, vienna_str, path):
    '''
    Constructs and exports an SVG for the annotationg for a specific strain.

    Parameters
    ----------
    strain_id: str
        The name of the strain to be presented on the tree
    pa_seq: str
        The polyA hairpin sequence
    vienna_str: str
        The vienna string for the predicted secondary structure of the polyA hairpin
    path: str
        The output directory where the SVG is saved. 
    
    Returns
    -------
    None
    '''

    LETTER_COLOR = {'A':'green', 'G':'black', 'C': 'blue', 'U':'red'}

    ext_pa_seq = pa_seq
    #ext_pa_seq = tstrain_info['locorna_pa']
    vstr = vienna_str

    font_size = 50
    y = 55
    x = 0
    dwg = svgwrite.Drawing(path + strain_id + '_face.svg', profile='full', size=(5000, 60))
    dwg.add(dwg.text('    ', x=[x], y=[y], font_size=font_size))
    x += 40



    for i in range(len(ext_pa_seq)):
        color = 'black'
        if ext_pa_seq[i] in LETTER_COLOR:
            color = LETTER_COLOR[ext_pa_seq[i]]
        if ext_pa_seq[i] != '-':
            pos = i - str(ext_pa_seq[0:i]).count('-')
            size=(font_size, font_size)

            if(vstr[pos] == '.'):
                #dwg.add(dwg.rect(insert=(x, y-(font_size)), size=(font_size, font_size), fill=color, fill_opacity=0.25, stroke=color, stroke_width=5))
                dwg.add(dwg.rect(insert=(x, y-(font_size)), size=size, fill=color))
            else:
                #dwg.add(dwg.rect(insert=(x, y-(font_size)), size=(font_size, font_size), fill=color))
                dwg.add(dwg.rect(insert=(x, y-(font_size)), size=size, fill=color, fill_opacity=0.25, stroke=color, stroke_width=0))
        else:
            dwg.add(dwg.line(start=(x + (.25*font_size), .5*font_size), end=(x+(.75*font_size),.5*font_size), stroke=color, stroke_width=15))
            #dwg.add(dwg.rect(insert=(x, y-(.5*font_size)), size=(font_size, font_size*.15), stroke=color))
            #dwg.add(dwg.text('-', x=[x + .25*font_size], y=[y], font_size=font_size,  fill='black'))
            #dwg.add(dwg.rect(insert=(x, y-(font_size)), size=(font_size, font_size), fill='white'))

        #dwg.add(dwg.text(ext_pa_seq[i], x=[x], y=[y], font_size=font_size,  fill=LETTER_COLOR[ext_pa_seq[i]]))
        

        #if(vstr[pos] == '.'):
        #    dwg.add(dwg.rect(insert=(x, y-font_size), size=(font_size, font_size), fill='yellow'))
        #else:
        #    dwg.add(dwg.text(ext_pa_seq[i], x=[x], y=[y], font_size=font_size))
        x += font_size*1.3    
    x += font_size

    #dwg.add(dwg.text(tstrain_info['free_energy'], x=[x], y=[y], font_size=font_size, fill='black'))
    dwg.save()

def draw_tree_freq_plot(aligned_records, output_file=output_root.format(item='freq_bars.svg'), LETTER_COLOR = {'A':'green', 'G':'black', 'C': 'blue', 'U':'red', '-':'white'}):
    '''
    Draw the frquency plot for the alinged records. This svg plot should be placed under the rendered tree.

    Parameters
    ----------
    aligned_records_path: [str]
        List of the aligned sequences annotated on the tree
    output_file: str
        The path where the svg will be saved.
    LETTER_COLOR
        A dictionary holding the colors for each nucleotide symbol
    
    Returns
    -------
    None - output is saved to svg file. 
    '''

    #Parse out aligned records
    seqs = aligned_records

    #Calculate frequencies
    psfm = {}
    for pos in range(len(seqs[0])):
        col = [pa[pos] for pa in seqs]
        psfm[pos] = {
            'A': float(col.count('A'))/len(seqs),
            'C': float(col.count('C'))/len(seqs),
            'G': float(col.count('G'))/len(seqs),
            'U': float(col.count('U'))/len(seqs),
            '-':  float(col.count('-'))/len(seqs),
        }
    
    #Draw figure
    bar_h = 500
    bar_w = 100

    x = 0

    dwg = svgwrite.Drawing(output_file, profile='full', size=(7500, 1000))

    for freqs in psfm.values():


        #Proportional heights for each nucleotide

        prop_heights = {
            '-': freqs['-'] * bar_h,
            'A':freqs['A'] * bar_h,
            'G':freqs['G'] * bar_h,
            'C':freqs['C'] * bar_h,
            'U':freqs['U'] * bar_h
        }

        #prop_heights = dict(sorted(prop_heights.items(), key = lambda x: x[1], reverse=True))

        y = 0
        #dwg.add(dwg.rect(insert=(x, y), size=(bar_w, bar_h), fill='white', stroke='black', stroke_width=5))

        for c in prop_heights.keys():
            if prop_heights[c] == 0:
                continue

            color = LETTER_COLOR[c]

            insert_pos = (x, y)
            curr_size = (bar_w, prop_heights[c])

            dwg.add(dwg.rect(insert=insert_pos, size=curr_size, fill=color))

            y += prop_heights[c]

        x += bar_w * 1.3

    dwg.save()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Invalid arguments. Usage: make_tree.py [input_file.json]')
        quit()

    input_file = sys.argv[1]
    input_file = '../input/' + input_file
    inputs = json.load(open(input_file, 'r'))
    print('Reading input from ', str(input_file), '...')

    #--> input JSON formatted strain data
    strain_file = inputs['paths']['strain_info']

    #--> Parameters for ClustalW alignment and Tree formation
    clustal_path = inputs['alignment_tree']['clustal_path']
    phyml_path = inputs['alignment_tree']['phyml_path']
    bootstrap = inputs['params']['bootstrap']
    phylo_target = inputs['params']['phylo_gene']

    #Load target strains
    target_strains = json.load(open(strain_file, 'r'))
    
    #Get polymerase sequences for each strain
    phylo_target_output_file = output_root.format(item='phylo_target_seqs.fasta')
    get_target_seqs([(t['accession'] , str(t['subtype'] + '-' + t['accession'])) for t in target_strains], phylo_target_output_file)

    #Algin polymerase sequences
    print('Aligning' + phylo_target + ' sequences...')
    aligned_phylo_target_output_file = output_root.format(item='aligned_phylo_target_seqs.fasta')
    clustalw_cmd = ClustalwCommandline(cmd=clustal_path, infile=phylo_target_output_file, outfile=aligned_phylo_target_output_file, type='PROTEIN', output='FASTA')
    clustalw_cmd()

    aligned_phylo_target_output_file_phlip = output_root.format(item='aligned_phylo_target_seqs.phylip')
    clustalw_cmd = ClustalwCommandline(cmd=clustal_path, infile=phylo_target_output_file, outfile=aligned_phylo_target_output_file_phlip, type='PROTEIN', output='PHYLIP')
    clustalw_cmd()

    '''
    #Convert aligned polyermase sequence file to PHYLIP
    aligned_polymerase_output_file_phylip = output_root.format(item='aligned_polymerase_seqs.phylip')
    with open(aligned_polymerase_output_file_phylip,'w+') as file:
        SeqIO.write(SeqIO.parse(open(aligned_polymerase_output_file, 'r'), 'fasta'), file, 'phylip')'''

    #Make phylogenetic tree
    print('Constructing Maximum Likelihood Tree...')
    tree_file = output_root.format(item='tree.nh')
    phyml_cmd = PhymlCommandline(cmd=phyml_path, datatype='aa', input=aligned_phylo_target_output_file_phlip, bootstrap=inputs['params']['bootstrap'])
    phyml_cmd()
    os.rename(src=aligned_phylo_target_output_file_phlip +'_phyml_tree.txt', dst=tree_file)

    #Replacing the accessions to include the subtype in the tree
    tree_renaming = open(tree_file, 'r').readline()
    for t in target_strains:
        if t['accession'] in tree_renaming:
            tree_renaming = tree_renaming.replace(t['accession'], str(t['subtype'] + '-' + t['accession']))

    with open(tree_file, 'w') as file:
        file.write(tree_renaming)


    #--> Generate the figure for the tree

    #Make the faces for each strain
    faces_path = output_root.format(item='faces/')
    os.mkdir(faces_path)

    for t in target_strains:
        make_face(strain_id=str(t['subtype'] + '-' + t['accession']), pa_seq=t['aligned_seq'], vienna_str=t['pa_vstr'], path=faces_path)
     
    #Open the tree
    tree = Tree(tree_file, format=0)
    print(tree)

    #Add the image faces to each leaf
    for n in tree.traverse():
        if not n.is_leaf():
            continue

        nface = TextFace(n.name, fsize=100)
        nface.hz_align = 0
        nface.margin_left = 0
        nface.margin_top = 5
        nface.margin_bottom = 5
        nface.tight_text = True
        n.add_face(nface, column=0, position='aligned')

        f = faces.ImgFace(faces_path + n.name + '_face.svg')
        n.add_face(f, column=1, position='aligned')

        for strain in target_strains:
            if strain['accession'] in n.name:
                fe_face = TextFace(str(strain['pa_free_energy']), fsize=100)
                n.add_face(fe_face, column=2, position='aligned')
    

    #Set the distances for the tree and set the thickeness to represent the bootstrap value
    init_dist = 75
    for n in tree.traverse():
        if n.is_leaf():
            continue
        n.dist = init_dist

    max_dist = (tree.get_distance(tree.get_farthest_leaf()[0], topology_only=True) * init_dist) + init_dist + 3
    for n in tree.get_leaves():
        n.dist = max_dist - sum([an.dist for an in n.get_ancestors()])

    max_thick = 30
    for n in tree.traverse():

        ns = NodeStyle()
        ns['vt_line_color'] = 'black'
        ns['hz_line_color'] = 'black'
        ns['vt_line_type'] = 0
        ns['hz_line_type'] = 0
        ns['vt_line_width'] = max(max_thick * (n.support/1000), 5)
        ns['hz_line_width'] =  max(max_thick * (n.support/1000), 5)
        n.img_style = ns
    
    #Style the tree and render
    ts = TreeStyle()
    #ts.force_topology = True
    ts.show_leaf_name = False
    ts.complete_branch_lines_when_necessary = True
    #ts.draw_guiding_lines = True
    #ts.guiding_lines_type = 0
    #ts.guiding_lines_color = 'black'
    ts.branch_vertical_margin = 50
    ts.min_leaf_separation = 150
    ts.scale = 1
    #t.render('../output/' + OUPUT_ID + '/tree3_locorna_blocks_switch-color_large.svg', w=5000, h=5000, tree_style=ts)
    tree.render(output_root.format(item='tree.svg'), w=5000, h=5000, tree_style=ts)

    #Draw the frequency plot to be annotated at the bottom of the tree
    draw_tree_freq_plot(aligned_records=[t['aligned_seq'] for t in target_strains], output_file=output_root.format(item='freq_bars.svg'))
