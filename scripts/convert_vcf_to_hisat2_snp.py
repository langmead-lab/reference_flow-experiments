"""Convert variants in a VCF file to HISAT2 SNP format.

Example:
convert_vcf_to_hisat2_snp.py -v wg_filtered.q0.1.vcf \
    -s wg_filtered.q0.1.hisat2.snp
"""
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='Path to the VCF file to convert.'
    )
    parser.add_argument(
        '-s', '--snp',
        help='Output snp file.'
    )
    return parser.parse_args()

def convert_vcf_to_hisat2_snp_line(line):
    if line[0] != '#':
        line = line.split()
        if len(line[0].split(',')) > 1:
            print('multi-allelic', line[0], line[1])
            return None
        vid = line[0] + '_' + line[1]
        if len(line[3]) == len(line[4]):
            vtype = 'single'
            alt = line[4]
        elif len(line[3]) > len(line[4]):
            vtype = 'deletion'
            alt = len(line[3]) - len(line[4])
        elif len(line[3]) < len(line[4]):
            vtype = 'insertion'
            alt = line[4][len(line[3]):]
        else:
            raise ValueError('Invalid type', line)
        chrom = line[0]
        pos = int(line[1]) - 1
        return f'{vid}\t{vtype}\t{chrom}\t{pos}\t{alt}\n'
    return None


def convert_vcf_to_hisat2_snp(snp, vcf):
    fw = open(snp, 'w')
    f = open(vcf, 'r')
    for line in f:
        output = convert_vcf_to_hisat2_snp_line(line)
        if output:
            fw.write(output)

if __name__ == '__main__':
    args = parse_args()
    convert_vcf_to_hisat2_snp(args.snp, args.vcf)

