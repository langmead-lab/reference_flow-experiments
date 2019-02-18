'''
Fuctions to process vcf file
'''
import argparse
import sys

class VCF:
    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
    v_chrom = ''
    v_pos = 0
    v_id = ''
    ref_allele = ''
    alt_allele = ''
    v_qual = 0.0
    v_filter = ''
    v_info = ''
    v_format = ''
    v_type = ''
    v_af = None
    
    def __init__(self, row, alt_id=0):
        self.v_chrom = row[0]
        self.v_pos = int(row[1])
        self.v_id = row[2]
        self.ref_allele = row[3]
        assert len(self.ref_allele.split(',')) == 1
        self.alt_allele = row[4].split(',')[alt_id]
        self.v_qual = float(row[5])
        self.v_filter = row[6]
        self.v_info = row[7]
        self.v_format = row[8]
        
        for info in self.v_info.split(';'):
            if info.startswith('AF'):
                af = info[info.find('=') + 1:]
                self.v_af = float(af.split(',')[alt_id])

        # get var type
        if len(self.alt_allele) == len(self.ref_allele):
            if len(self.v_type) > 0:
                self.v_type += ','
            self.v_type += 'SNP'
        # INS
        elif len(self.alt_allele) > len(self.ref_allele):
            if len(self.v_type) > 0:
                self.v_type += ','
            self.v_type += 'INS'
        # DEL
        elif len(self.alt_allele) < len(self.ref_allele):
            if len(self.v_type) > 0:
                self.v_type += ','
            self.v_type += 'DEL'
        else:
            print ('Error: unexpected variant type :', self.ref_allele, self.alt_allele)
            exit ()
    
    def print(self, chrm=False, pos=True, ref_allele=True, alt_allele=True, af=True):
        msg = ''
        if chrm:
            msg += self.v_chrom
            msg += ' '
        if pos:
            msg += str(self.v_pos)
            msg += ' '
        if ref_allele:
            msg += self.ref_allele
            msg += ' '
        if alt_allele:
            msg += str(self.alt_allele)
            msg += ' '
        if af:
            msg += str(self.v_af)
            msg += ' '
        print (msg)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--vcf',
        help='vcf file'
    )
    parser.add_argument(
        '--out_var_loc',
        default=None,
        help='output file storing filtered variant locations; no output if this argument is not specified [None]'
    )
    parser.add_argument(
        '--out_no_conflicting_vcf',
        default=None,
        help='output vcf file with conflicting variants removed; no output if this argument is not specified [None]'
    )
    parser.add_argument(
        '--min_af',type=float,default=0.0,
        help='min allele frequency to be considered [0.0]'
    )
    args = parser.parse_args()
    return args

def build_vcf(line):
    if line.startswith('##'):
        return
    #: header
    if line.startswith('#'):
        return
    row = line.split()
    num_alt_alleles = len(row[4].split(','))
    list_vcf = []
    for i in range(num_alt_alleles):
        list_vcf.append(VCF(row, i))
    return list_vcf

def comp_var_with_list(var, list_var):
    list_comp = [None] * len(list_var)
    for list_idx, v in enumerate(list_var):
        offset = var.v_pos - v.v_pos
        #: check if ref allele is the same
        ref_allele_check = True 
        # if ref_a != v.ref_allele[offset:]:
        #     ref_allele_check = False
        for idx, ref_a in enumerate(var.ref_allele):
            try:
                if ref_a != v.ref_allele[offset + idx]:
                    ref_allele_check = False
            except:
                break
        assert ref_allele_check == True
        #: check if alt allele is the same
        alt_allele_check = True
        for idx, alt_a in enumerate(var.alt_allele):
            try:
                if alt_a != v.alt_allele[offset + idx]:
                    alt_allele_check = False
            except:
                break
        list_comp[list_idx] = alt_allele_check
    return list_comp

def remove_conflicting_vars(vcf_fn, out_vcf_fn, min_af):
    vcf_f = open(vcf_fn, 'r')
    prev_range = range(0)
    prev_var = []
    for line in vcf_f:
        if line.startswith('#'):
            continue
        list_vcf = build_vcf(line)
        for vcf in list_vcf:
            #: ignore var with allele freq less than given threshold
            #: keep var if allele freq is not specified in vcf
            if vcf.v_af != None and vcf.v_af < min_af:
                continue
            
            if vcf.v_pos in prev_range:
                for pv in prev_var:
                    pv.print()
                vcf.print()
                comp = comp_var_with_list(vcf, prev_var)
                input(comp)
                prev_var.append(vcf)
            if len(prev_range) == 0 or vcf.v_pos > max(prev_range):
                prev_range = range(vcf.v_pos, vcf.v_pos + len(vcf.ref_allele))
                prev_var = [vcf]


if __name__ == '__main__':
    args = parse_args()
    vcf_fn = args.vcf
    out_vcf_fn = args.out_no_conflicting_vcf
    min_af = args.min_af
    remove_conflicting_vars(vcf_fn, out_vcf_fn, min_af)
