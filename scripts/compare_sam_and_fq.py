'''
This script compares different sets of sequences and reports the number of intersection/union.
Supported formats:
    .sam
    .fq (.fastq)
    .names (just a list of read names)
'''
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '-set1', '--set1_fn',
    help='set1, can be .sam, .fq .or .names'
)
parser.add_argument(
    '-set2', '--set2_fn',
    help='set2, can be .sam, .fq .or .names'
)
args = parser.parse_args()

set1_fn = args.set1_fn
set2_fn = args.set2_fn
set1_f = open(set1_fn, 'r')
set2_f = open(set2_fn, 'r')

list_set1_reads = []
if set1_fn.endswith('.fq') or set1_fn.endswith('.fastq'):
    for line in set1_f:
        if line.startswith('@'):
            list_set1_reads.append(line.split()[0][1:])
elif set1_fn.endswith('.sam'):
    for line in set1_f:
        list_set1_reads.append(line.split()[0])
elif set1_fn.endswith('.names'):
    for line in set1_f:
        list_set1_reads.append(line.rstrip())
else:
    print ('unrecognized file format')
    exit()
set_set1_reads = set(list_set1_reads)

list_set2_reads = []
if set2_fn.endswith('.fq'):
    for line in set2_f:
        if line.startswith('@'):
            list_set2_reads.append(line.split()[0][1:])
elif set2_fn.endswith('.sam'):
    for line in set2_f:
        list_set2_reads.append(line.split()[0])
elif set2_fn.endswith('.names'):
    for line in set2_f:
        list_set2_reads.append(line.rstrip())
else:
    print ('unrecognized file format')
    exit()
set_set2_reads = set(list_set2_reads)

print ('size of set1', len(set_set1_reads))
print ('size of set2', len(set_set2_reads))
print ('size of intersection', len(set_set1_reads.intersection(set_set2_reads)))
print ('size of union', len(set_set1_reads.union(set_set2_reads)))
