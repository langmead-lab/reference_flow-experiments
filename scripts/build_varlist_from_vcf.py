'''
Parses a VCF (V4.2) file and builds a list of variants.

Also supports the comparison between a golden var list.
'''
import argparse
from build_erg import read_var

class VCF4_2:
    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
    v_chrom = ''
    v_pos = 0
    v_id = ''
    ref_allele = ''
    alt_alleles = []
    num_alt_alleles = 1
    v_qual = 0.0
    v_filter = ''
    v_info = ''
    v_format = ''
    v_type = ''
    
    def __init__(self, line):
        row = line.rstrip().split('\t')
        self.v_chrom = row[0]
        self.v_pos = int(row[1])
        self.v_id = row[2]
        self.ref_allele = row[3]
        self.alt_alleles = row[4].split(',')
        self.v_qual = float(row[5])
        self.v_filter = row[6]
        self.v_info = row[7]
        self.v_format = row[8]
        assert len(self.ref_allele.split(',')) == 1
        self.num_alt_alleles = len(self.alt_alleles)
        # duplicates if there are > 1 alt alleles
        for alt_allele in self.alt_alleles:
            # SNP
            if len(alt_allele) == len(self.ref_allele):
                if len(self.v_type) > 0:
                    self.v_type += ','
                self.v_type += 'SNP'
            # INS
            elif len(alt_allele) > len(self.ref_allele):
                if len(self.v_type) > 0:
                    self.v_type += ','
                self.v_type += 'INS'
            # DEL
            elif len(alt_allele) < len(self.ref_allele):
                if len(self.v_type) > 0:
                    self.v_type += ','
                self.v_type += 'DEL'
            else:
                print ('Error: unexpected variant type :', self.ref_allele, alt_allele)
                exit ()
        

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='target vcf file'
    )
    parser.add_argument(
        '--min_qual',
        default=0.0, type=float,
        help='vcf quality threshold [0.0]'
    )
    parser.add_argument(
        '-vg', '--golden_vcf',
        default='',
        help='golden vcf file'
    )

    args = parser.parse_args()
    return args

def build_varlist_from_vcf(args):
    '''
    ##fileformat=VCFv4.2
    '''
    target_vcf_fn = args.vcf
    golden_vcf_fn = args.golden_vcf
    min_qual = float(args.min_qual)
    target_vcf_f = open(target_vcf_fn, 'r')
    
    var_list = read_var(golden_vcf_fn, remove_conflict=True, remove_coexist=False)
    print ('len of golden list:', len(var_list))

    num_vars = 0
    num_true_pos = 0

    golden_idx = 0
    for line in target_vcf_f:
        # Skip header lines (##...) and info line (#...)
        if line[0] == '#':# and line[1] == '#'
            continue
        vcf = VCF4_2(line)

        # filtered by QUAL
        if vcf.v_qual < min_qual:
            continue
        
        for i, alt_allele in enumerate(vcf.alt_alleles):
            # a diploid variant is counted twice
            num_vars += 1
            # v_type = vcf.v_type.split(',')[i]
            # print (vcf.v_chrom, v_type, vcf.v_pos, vcf.ref_allele, alt_allele)
            for j in range(golden_idx, len(var_list)):
                if var_list[j].chrm == vcf.v_chrom and \
                var_list[j].ref_pos == vcf.v_pos and \
                var_list[j].ref_allele == vcf.ref_allele and \
                var_list[j].alt_allele == alt_allele:
                    golden_idx = j
                    num_true_pos += 1
                    break
                if var_list[j].chrm == vcf.v_chrom and \
                    var_list[j].ref_pos > vcf.v_pos:
                    golden_idx = max(j - 1, 0)
                    break
    
    # Stats
    precision = 100 * float(num_true_pos) / num_vars
    sensitivity = 100 * float(num_true_pos) / len(var_list)
    
    print ('PRECISION = %.2f%%' % precision)
    print ('SENSITIVITY = %.2f%%' % sensitivity)
    

if __name__ == '__main__':
    args = parse_args()
    build_varlist_from_vcf(args)
