import pandas as pd
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--raw_rmsk_file',
        default='rmsk.txt',
        help='Path to rmsk.txt'
    )
    parser.add_argument(
        '-l', '--level',
        help='RepeatMasker annotation level [repName, repClass, repFamily]'
    )
    parser.add_argument(
        '-r', '--repeat_label',
        help='Queried repeat label'
    )
    parser.add_argument(
        '-o', '--output_file', default=sys.stdout,
        help='Path to output BED file [stdout]'
    )
    args = parser.parse_args()
    return args

def extract_repeats_from_rmsk(args):
    df = pd.read_csv(args.raw_rmsk_file, sep='\t', header=None)
    df.columns = ['bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'genoName', 'genoStart', 'genoEnd', 'genoLeft', 'strand', 'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'repLeft', 'id']
    df_extracted = df[df[args.level] == args.repeat_label][['genoName', 'genoStart', 'genoEnd', 'repName', 'strand', 'repClass', 'repFamily']]
    df_extracted.to_csv(args.output_file, sep='\t', header=None, index=None)

if __name__ == '__main__':
    args = parse_args()
    extract_repeats_from_rmsk(args)
