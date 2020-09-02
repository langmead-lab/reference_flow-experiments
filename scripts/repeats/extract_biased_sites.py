'''Reads an allele.bias file and extracts biased sites.
'''
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--bias_file',
        help='A *.allele.bias file'
    )
    parser.add_argument(
        '-o', '--output_file',
        help='Output file in BED format'
    )
    parser.add_argument(
        '-t', '--threshold', type=float, default=0.3,
        help='A float. A site with `abs(v-0.5) >= threshold` \
                is defined as biased. v is its allelic balance'
    )
    parser.add_argument(
        '--min_reads', type=int, default=15,
        help='An integer. The min read coverage for each site. Sites with lower depths are ignored'
    )
    args = parser.parse_args()
    return args

def extract_biased_sites(arg):
    df = pd.read_csv(args.bias_file, sep='\t')
    df_high_coverage = df[df['NUM_READS'] >= args.min_reads]
    df_high_coverage_biased = df_high_coverage[
        abs(df_high_coverage['REFERENCE_BIAS'] - 0.5) >= args.threshold]
    with open(args.output_file, 'w') as f:
        for i in range(df_high_coverage_biased.shape[0]):
            site = df_high_coverage_biased.iloc[i]
            chromosome = site['CHR']
            position = site['HET_SITE']
            f.write(f'{chromosome}\t{position-1}\t{position}\n')

if __name__ == '__main__':
    args = parse_args()
    extract_biased_sites(args)

