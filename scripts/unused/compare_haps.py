'''
Last update: 2018/11/28 by Nae-Chyun Chen

Compares two haplotypes (from .fa files) and 
reports variants.
'''

import sys

# inputs
hap1_fn = '/scratch/groups/blangme2/naechyun/relaxing/na12878/na12878-chr9-phase3_hapA.fa'
hap2_fn = '/scratch/groups/blangme2/naechyun/relaxing/na12878/na12878-chr9-phase3_hapB.fa'

# laod fasta
def read_fasta(fn):
    seq = ''
    with open(fn, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq += line[: line.find('\\')]
    return seq
hap1_seq = read_fasta(hap1_fn)
hap2_seq = read_fasta(hap2_fn)

# check length
if len(hap1_seq) != len(hap2_seq):
    sys.stderr.write('WARNING: different lengths')

# report variants
# pos is 1-based, as 1KG vcf representation
print ('pos', 'hap1', 'hap2')
for i, nuc in enumerate(hap1_seq):
    if hap2_seq[i] != nuc:
        print (i + 1, nuc, hap2_seq[i])
