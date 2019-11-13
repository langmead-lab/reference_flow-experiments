'''
Reads a VCF and reports a filtered VCF including only the variants
carried by a given individual and a given haplotype (e.g NA12878_hapA)
'''

import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf', type=str, help="Path to VCF file [None]"
    )
    parser.add_argument(
        '-s', '--indiv', type=str, help="Name of individual in VCF to process; cannot be blank [None]"
    )
    parser.add_argument(
        '-hap', '--haplotype', type=str, help="Which haplotype (A or B)? [None]"
    )
    parser.add_argument(
        '-o', '--output', type=str, help="Path to output VCF file [None]"
    )
    args = parser.parse_args()
    return args

def filter_vcf_by_haplotype(fn_vcf, indiv, haplotype, fn_output):
    f = open(fn_vcf, 'r')
    f_out = open(fn_output, 'w')

    for line in f:
        if line[0] == '#' and line[1] == '#':
            f_out.write(line)
            continue
        if line[0] == '#':
            labels = line.rstrip().split('\t')
            #: if "indiv" is set, select corresponding columns
            col = None
            for i in range(9, len(labels)):
                if labels[i] == indiv:
                    col = i
            if not col:
                print('Error! Could not find individual {} in VCF'.format(indiv))
                exit(1)
            reduced_labels = labels[:9]
            reduced_labels.append(labels[col])
            reduced_line = '\t'.join(reduced_labels) + '\n'
            # print (reduced_line)
            f_out.write(reduced_line)
            continue

        fields = line.rstrip().split('\t')
        gt = fields[col].split('|')
        out_fields = fields[:9]
        out_fields.append(fields[col])
        if (haplotype == 'A' and gt[0] != '0') or \
            (haplotype == 'B' and gt[1] != '0'):
            out_line = '\t'.join(out_fields) + '\n'
            f_out.write(out_line)


if __name__ == '__main__':
    args = parse_args()
    fn_vcf = args.vcf
    indiv = args.indiv
    haplotype = args.haplotype
    fn_output = args.output

    assert fn_vcf != None
    assert indiv != None
    assert haplotype in ['A', 'B']
    assert fn_output != None

    filter_vcf_by_haplotype(fn_vcf, indiv, haplotype, fn_output)
