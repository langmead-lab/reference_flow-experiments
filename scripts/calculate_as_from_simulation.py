import pandas as pd
import numpy as np
import argparse
import re
import math

parser = argparse.ArgumentParser()
parser.add_argument(
    '-n', '--fn_sam',
    help='target sam file'
)
parser.add_argument(
    '-g', '--fn_gold',
    help='golden .fq.sam file'
)
parser.add_argument(
    '-read_len', '--read_len', type=int, default=100,
    help='read length (100)'
)
args = parser.parse_args()
read_len = args.read_len
fn_gold = args.fn_gold
fn_sam = args.fn_sam

df = pd.read_csv(fn_gold, sep='\t', header=None)
df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NM', 'MD', 'XE', 'XS', 'XI']

# processes MD string
split_md = []
for i in range(df.shape[0]):
    split_md.append(re.findall('([0-9]+|[a-zA-Z]+)', df['MD'].iloc[i].split(':')[-1]))
df['sMD'] = split_md

# processes CIGAR string
split_cigar = []
for i in range(df.shape[0]):
    split_cigar.append(re.findall('([0-9]+[a-zA-Z]+)', df['CIGAR'].iloc[i]))
# df['cigar'] = split_cigar

reduced_cigar = []
for i in range(df.shape[0]):
    # ci = df['cigar'][i]
    ci = split_cigar[i]
    for c in ci:
        if c.count('P') > 0:
            ci.remove(c)
    
    ci_flag = True
    while ci_flag:
        ci_flag = False
        for j, c in enumerate(ci):
            if j == len(ci)-1:
                break
            elif ci[j+1][-1] == c[-1]:
                ci_flag = True
                count = int(c[:-1]) + int(ci[j+1][:-1])
                ci[j] = str(count) + c[-1]
                del ci[j+1]

    #: checks if length of cigar (I+M) matches read length
    len_ci = 0
    for c in ci:
        assert c[-1] in ['I','D','M']
        if c[-1] != 'D':
            len_ci += int(c[:-1])
    assert len_ci == read_len
    reduced_cigar.append(ci)
df['rCIGAR'] = reduced_cigar

list_as = [None] * df.shape[0]
for i in range(df.shape[0]):
    #: reads with no error
    if int(df['XE'].iloc[i].split(':')[-1]) == 0:
        list_as[i] = 0
    else:
        score = 0
        vec_aln = []
        for c in df['rCIGAR'].iloc[i]:
            #: "0" = match/mismatch
            if c[-1] == 'M':
                for z in range(int(c[:-1])):
                    vec_aln.append(0)
            elif c[-1] in ['I', 'D']:
                len_indel = int(c[:-1])
                score -= (5 + 3 * len_indel)
                #: "1" = insertion
                if c[-1] == 'I':
                    for z in range(int(c[:-1])):
                        vec_aln.append(1)
                #: "2" = deletion
                else:
                    for z in range(int(c[:-1])):
                        vec_aln.append(2)
        current_pos = 0
        current_count_del = 0
        current_cigar_pos = current_pos + current_count_del
        for c in df['sMD'].iloc[i]:
            if c.isalpha():
                for a in c:
                    if vec_aln[current_cigar_pos] == 2:
                        current_count_del += 1
                        # print ('D at {}'.format(current_pos))
                    else:
                        qual = df['QUAL'].iloc[i][current_pos]
                        qual = ord(qual) - 33
                        score -= (2 + math.floor(4 * min(qual, 40)/40))
                        current_pos += 1
                    current_cigar_pos = current_pos + current_count_del
            else:
                #: checks INS
                overlapping_ins = vec_aln[current_cigar_pos: current_cigar_pos + int(c) + 1].count(1)
                if overlapping_ins > 0:
                    overlapping_ins = vec_aln[current_cigar_pos: current_cigar_pos + overlapping_ins + int(c) + 1].count(1)
                current_pos += int(c) + overlapping_ins
            current_cigar_pos = current_pos + current_count_del
        assert current_pos == read_len
        list_as[i] = score
df['AS'] = list_as

#: TODO: tags may differ, need parse line by line
#: compares with another sam file
df_sam = pd.read_csv(fn_sam, sep='\t', header=None)
df_sam_mapq_and_score = df_sam[[0,4,11,12]]
df_sam_mapq_and_score.columns = ['QNAME', 'MAPQ', 'AS', 'XS']
list_align_as = [int(i.split(':')[-1]) for i in df_sam_mapq_and_score['AS']]
df_sam_mapq_and_score['AS'] = list_align_as
list_align_xs = [int(i.split(':')[-1]) for i in df_sam_mapq_and_score['XS']]
df_sam_mapq_and_score['XS'] = list_align_xs
df_merge = df.merge(df_sam_mapq_and_score, on='QNAME')

print (sum(df_merge['AS_x'] <= df_merge['AS_y']))
print (sum(df_merge['AS_x'] > df_merge['AS_y']))