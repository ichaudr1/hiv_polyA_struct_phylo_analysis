from polyA_parsing import gag_position
from tqdm import tqdm


path = 'zarudnaya_raw.txt'

all_genomes = []

with open(path, 'r') as file:
    all_lines = file.readlines()
    for l in all_lines:
        cl = l.split(',')
        all_genomes.extend(s.strip() for s in cl if len(s.strip()) > 0)

print('Testing ' , str(len(all_genomes)) + ' genomes...')

count = 0
a = []
found = []
for g in tqdm(all_genomes):
    tqdm.write(g)
    i = gag_position(g)
    if i < 0:
        tqdm.write('---> Not found for: ' + g)
        a.append(g)
        count += 1
    else:
        found.append((g, i))

print(count)
with open('gag_position_detection_zarudnaya.txt', 'w+') as file:
    file.write('Gag not found in ' + str(count) + '\n') 
    file.write(str(a))

with open('gag_position_detection_zarudnaya_found.txt', 'w+') as file:
    for item in found:
        file.write(str(item[0]) + ', ' + str(item[1]) + '\n')
