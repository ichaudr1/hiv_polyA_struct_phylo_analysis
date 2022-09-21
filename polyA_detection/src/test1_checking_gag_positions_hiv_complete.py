from polyA_parsing import gag_position
import json
from tqdm import tqdm


path = '/Users/ichaudr/Documents/programming_projects/hiv_rna_struct_phylo_analysis/datasets/hiv_complete.json'

all_genomes = json.load(open(path, 'r'))

count = 0
a = []
for g in tqdm(all_genomes):
    i = gag_position(g['accession'])
    if i < 0:
        tqdm.write(g['accession'])
        a.append(g['accession'])
        count += 1

print(count)
with open('gag_position_detection.txt', 'w+') as file:
    file.write('Gag not found in ' + str(count))
    file.write(a)
