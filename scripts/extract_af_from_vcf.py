'''
Reads a VCF file and extracts the AF fields in the output file
'''
import argparse

def extract_af_from_vcf(vcf, af):
    f = open(vcf, 'r')
    records = []
    for line in f:
        if line[0] != '#':
            records.append(line)
    infos = [r.split()[7].split(';') for r in records]
    afs = []
    for info in infos:
        for i in info:
            if i.startswith('AF='):
                afs.append(float(i.split('=')[1]))
    fw = open(af, 'w')
    for i in afs:
        fw.write(f'{i}\n')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fn_vcf',
        help='VCF file.')
    parser.add_argument(
        '--fn_af',
        help='Output filename for allele frequencies.')
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_args()
    extract_af_from_vcf(args.fn_vcf, args.fn_af)
