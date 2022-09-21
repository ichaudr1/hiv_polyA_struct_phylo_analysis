import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from tqdm import tqdm


strains = json.load(open('raw_output.json', 'r'))

to_write = []

for s in tqdm(strains):
    if s['pa_seq'] == None:
        continue
    to_write.append(SeqRecord(Seq(s['pa_seq']), id=s['subtype'] + '-' + s['accession']))

SeqIO.write(to_write, open('all_pas.fasta', 'w+'), 'fasta')