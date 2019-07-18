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
    help='output filename, extracted intersection based on set2. \
    When -o is set, -set2 needs to be a sam/fq file [None]'
)
parser.add_argument(
    '--mason2',
    type=int, default=0,
    help='set to ignore texts after "/" of a read name. \
            0: use full read name, \
            1: ignore for set1, \
            2: ignore for set2, \
            12, ignore for set1 and set2 [(INT) 0]'
)
args = parser.parse_args()

fn_set1 = args.fn_set1
fn_set2 = args.fn_set2
fn_output = args.fn_output
mason2 = args.mason2
assert mason2 in [0, 1, 2, 12]

f_set1 = open(fn_set1, 'r')
f_set2 = open(fn_set2, 'r')

if fn_output != None:
    assert fn_set2.endswith('.sam') or fn_set2.endswith('.fq')
    f_output = open(fn_output, 'w')

list_set1_reads = []
if fn_set1.endswith('.fq') or fn_set1.endswith('.fastq'):
    for line in f_set1:
        if line.startswith('@'):
            read_name = line.split()[0][1:]
            if mason2 in [1, 12]:
                read_name = read_name.split('/')[0]
            list_set1_reads.append(read_name)
elif fn_set1.endswith('.sam'):
    for line in f_set1:
        #: skips header
        if line[0] == '@':
            continue

        read_name = line.split()[0]
        if mason2 in [1, 12]:
            read_name = read_name.split('/')[0]
        list_set1_reads.append(read_name)
elif fn_set1.endswith('.names'):
    for line in f_set1:
        read_name = line.rstrip()
        if mason2 in [1, 12]:
            read_name = read_name.split('/')[0]
        list_set1_reads.append(read_name)
else:
    print ('Error: unrecognized file format: {}'.format(fn_set1))
    exit()
set_set1_reads = set(list_set1_reads)

list_set2_reads = []
if fn_set2.endswith('.fq'):
    fq_flag = 0
    for line in f_set2:
        if line.startswith('@'):
            read_name = line.split()[0][1:]
            if mason2 in [2, 12]:
                read_name = read_name.split('/')[0]
            list_set2_reads.append(read_name)
            if fn_output != None:
                if read_name in set_set1_reads:
                    fq_flag = 1
        if fq_flag > 0:
            f_output.write(line)
            fq_flag += 1
            if fq_flag > 4:
                fq_flag = 0
elif fn_set2.endswith('.sam'):
    for line in f_set2:
        #: skips header
        if line[0] == '@':
            continue

        read_name = line.split()[0]
        if mason2 in [2, 12]:
            read_name = read_name.split('/')[0]
        list_set2_reads.append(read_name)
        if fn_output != None:
            if read_name in set_set1_reads:
                f_output.write(line)
elif fn_set2.endswith('.names'):
    for line in f_set2:
        list_set2_reads.append(line.rstrip())
else:
    print ('Error: unrecognized file format: {}'.format(fn_set2))
    exit()
set_set2_reads = set(list_set2_reads)

print ('Size of set1', len(set_set1_reads))
print ('Size of set2', len(set_set2_reads))
print ('Size of intersection', len(set_set1_reads.intersection(set_set2_reads)))
print ('Size of union', len(set_set1_reads.union(set_set2_reads)))
