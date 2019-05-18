import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument(
    '-sam', '--fn_sam',
    help='target sam file'
)
args = parser.parse_args()
fn_sam = args.fn_sam

df = pd.read_csv(fn_sam, sep='\t', header=None)
split_cigar = []
for i in range(df.shape[0]):
    split_cigar.append(re.findall('([0-9]+[a-zA-Z]+)', df[5].iloc[i]))
df['cigar'] = split_cigar
split_md = []
for i in range(df.shape[0]):
    split_md.append(re.findall('([0-9]+|[a-zA-Z]+)', df[12].iloc[i].split(':')[-1]))
df['md'] = split_md
