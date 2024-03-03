'''
Running on hal does not produce the final svg file due to to graphics interface. This script will take the faces and nh output to make the tree. 
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




folder = '../output/'
output_id = 'hal_make_tree_aca2a150-3d2e-11ed-825f-d57ceadc3c97'
output_root = folder + output_id + '/{item}'


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



#Open the target strains
strain_file = '/Users/ichaudr/Documents/programming_projects/hiv_rna_struct_phylo_analysis/polyA_detection/output/polyA_parse_08120dfe-3899-11ed-9239-a683e73331cb/aligned_pa_seqs.json'
target_strains = json.load(open(strain_file, 'r'))

#Open the tree
tree_file = output_root.format(item='tree.nh')
tree = Tree(tree_file, format=0)
print(tree)

faces_path = output_root.format(item='faces/')

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