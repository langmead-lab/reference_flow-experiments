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
    '-set1', '--fn_set1',
    help='set1, can be .sam, .fq .or .names'
)
parser.add_argument(
    '-set2', '--fn_set2',
    help='set2, can be .sam, .fq .or .names'
)
parser.add_argument(
    '-o', '--fn_output',
    help='output filename, extracted intersection based on set2. \n \
    When -o is set, -set2 needs to be a sam file [None]'
)
args = parser.parse_args()

fn_set1 = args.fn_set1
fn_set2 = args.fn_set2
fn_output = args.fn_output
f_set1 = open(fn_set1, 'r')
f_set2 = open(fn_set2, 'r')

if fn_output != None and (fn_set2.endswith('.sam') is False):
    print ('ERROR: -set2 needs to be a sam file when -o is set.')
    exit()
elif fn_output != None:
    f_output = open(fn_output, 'w')

list_set1_reads = []
if fn_set1.endswith('.fq') or fn_set1.endswith('.fastq'):
    for line in f_set1:
        if line.startswith('@'):
            list_set1_reads.append(line.split()[0][1:])
elif fn_set1.endswith('.sam'):
    for line in f_set1:
        list_set1_reads.append(line.split()[0])
elif fn_set1.endswith('.names'):
    for line in f_set1:
        list_set1_reads.append(line.rstrip())
else:
    print ('unrecognized file format')
    exit()
set_set1_reads = set(list_set1_reads)

list_set2_reads = []
if fn_set2.endswith('.fq'):
    for line in f_set2:
        if line.startswith('@'):
            list_set2_reads.append(line.split()[0][1:])
elif fn_set2.endswith('.sam'):
    for line in f_set2:
        list_set2_reads.append(line.split()[0])
        if fn_output != None:
            if line.split()[0] in set_set1_reads:
                f_output.write(line)
elif fn_set2.endswith('.names'):
    for line in f_set2:
        list_set2_reads.append(line.rstrip())
else:
    print ('unrecognized file format')
    exit()
set_set2_reads = set(list_set2_reads)

print ('size of set1', len(set_set1_reads))
print ('size of set2', len(set_set2_reads))
print ('size of intersection', len(set_set1_reads.intersection(set_set2_reads)))
print ('size of union', len(set_set1_reads.union(set_set2_reads)))
