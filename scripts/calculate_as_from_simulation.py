import pandas as pd
import argparse
import re
import math

parser = argparse.ArgumentParser()
parser.add_argument(
    '-sam', '--fn_sam',
    help='target sam file'
)
parser.add_argument(
    '-read_len', '--read_len', type=int, default=100,
    help='read length (100)'
)
args = parser.parse_args()
read_len = args.read_len
fn_sam = args.fn_sam

df = pd.read_csv(fn_sam, sep='\t', header=None)

# processes CIGAR string
split_cigar = []
for i in range(df.shape[0]):
    split_cigar.append(re.findall('([0-9]+[a-zA-Z]+)', df[5].iloc[i]))
df['cigar'] = split_cigar

# processes MD string
split_md = []
for i in range(df.shape[0]):
    split_md.append(re.findall('([0-9]+|[a-zA-Z]+)', df[12].iloc[i].split(':')[-1]))
df['md'] = split_md

reduced_cigar = []
for i in range(df.shape[0]):
    ci = df['cigar'][i]
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
df['reduced_cigar'] = reduced_cigar

list_as = [None] * df.shape[0]
for i in range(df.shape[0]):
    #: reads with no error
    if int(df[13].iloc[i].split(':')[-1]) == 0:
        list_as[i] = 0
    else:
        score = 0
        for c in df['reduced_cigar'].iloc[i]:
            if c[-1] in ['I', 'D']:
                len_indel = int(c[:-1])
                score -= (5 + 3 * len_indel)
        current_pos = 0
        for c in df['md'].iloc[i]:
            if c.isalpha():
                for a in c:
                    qual = df[10].iloc[i][current_pos]
                    qual = ord(qual) - 33
                    score -= (2 + math.floor(4 * min(qual, 40)/40))
                    current_pos += 1
                #: ignores QUAL
                # score -= (len(c) * 6)
                # current_pos += len(c)
            else:
                current_pos += int(c)
        list_as[i] = score
df['AS'] = list_as

#: TODO
#: comparison with another sam file