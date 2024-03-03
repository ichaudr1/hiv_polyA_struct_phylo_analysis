'''
@author ichaudr

Objective: Align each of the polyA hairpin structures to the consensus output by LOCARNA. The frequencies of each bp type is recorded for each position and used to plot the output graph.

'''

from hp_align import *
import json
import svgwrite
import uuid
import os
import sys

#--> Output

#Generate a unique id for the output folder
folder = '../output/'
output_id = 'hp_align_' + str(uuid.uuid1())
os.mkdir(folder + output_id)
output_root = folder + output_id + '/{item}'



if __name__ == "__main__":

    #Parse the input
    if len(sys.argv) != 2:
        print('Invalid arguments. Usage: compare_to_concensus.py [input_file.json]')
        quit()

    input_file = sys.argv[1]
    input_file = '../input/' + input_file
    inputs = json.load(open(input_file, 'r'))
    print('Reading input from ', str(input_file), '...')

    #Concensus info:
    conc_seq = inputs['concensus_sequence']
    conc_vstr = inputs['concensus_vstr']
    conc_pairs = get_pairs(conc_seq, conc_vstr)
    #print(len(conc_pairs))


    #Path to strains.json
    all_strains_path = inputs['target_strains_path']

    #Dictionary holding all the allignments
    all_alignments = {}
    all_pairs = {}
    for strain in json.load(open(all_strains_path, 'r')):
        name = strain['subtype'] + '-' + strain['accession']
        pa_seq = strain['pa_seq']
        vstr = strain['pa_vstr']

        if pa_seq == None or vstr == None:
            continue

        pairs = get_pairs(pa_seq, vstr)
        alignment = rhelix_align(conc_pairs, pairs)[0]
        print(name)
        print(alignment)
        print_rhelix_aln(conc_pairs, pairs, alignment)
        print('*' * 50)
        all_alignments[name] = alignment
        all_pairs[name] = pairs

    #Get frequences for each pos
    per_conc_pos_frqs = {}
    for i in range(len(conc_pairs)):
        per_conc_pos_frqs[i] = {}

    for strain in all_alignments.keys():
        alignment = all_alignments[strain]
        for a in alignment:
            conc_pos = a[0]
            target_pos = a[1]

            if conc_pos == None:
                continue

            if a[1] == None:
                target_pos = '-'
            else:
                target_pos = ''.join(all_pairs[strain][target_pos])

            if target_pos in per_conc_pos_frqs[conc_pos]:
                per_conc_pos_frqs[conc_pos][target_pos] += 1
            else:
                per_conc_pos_frqs[conc_pos][target_pos] = 1

    #Draw the bars for each position in the helix

    BP_COLORS = {
        'AU': str(inputs['colors']['AU']),
        'GC':str(inputs['colors']['GC']),
        'GU':str(inputs['colors']['GU'])
    }

    BOX_H = 75
    BOX_W = 10
    dwg = svgwrite.Drawing(output_root.format(item='freq_bars.svg'), profile='full', size=(10*BOX_W *(len(conc_pairs)), 2*BOX_H))

    x = 0
    y = 1.2*BOX_H
    for i in range(len(conc_pairs)):
        dwg.add(dwg.text(str(i), x=[x], y=[1.8*BOX_H]))

        #Calculate freqs
        freqs = {}
        freqs['invalid'] = 0
        for bp in BP_COLORS:
            freqs[bp] = 0

        curr_con_pos_counts = per_conc_pos_frqs[i]
        denom = sum(curr_con_pos_counts.values())

        for bp in freqs.keys():
            count = 0
            for cbp in curr_con_pos_counts.keys():
                if bp[0] in cbp and bp[1] in cbp:
                    count += curr_con_pos_counts[cbp]
            freqs[bp] = count / denom
        
        freqs['invalid'] = 1 - sum(freqs.values())

        print(i, conc_pairs[i], '\t', freqs)

        #Calculate the proportional heights for each bp type
        prop_heights = {}
        for bp in freqs:
            prop_heights[bp] = BOX_H * freqs[bp]
        
        y = 0
        for c in prop_heights.keys():
            if prop_heights[c] == 0:
                continue
            
            color = 'gray'
            if c in BP_COLORS.keys():
                color = BP_COLORS[c]

            insert_pos = (x, y)
            curr_size = (BOX_W, prop_heights[c])

            dwg.add(dwg.rect(insert=insert_pos, size=curr_size, fill=color))

            y += prop_heights[c]

        x += BOX_W * 2.5

    dwg.save()