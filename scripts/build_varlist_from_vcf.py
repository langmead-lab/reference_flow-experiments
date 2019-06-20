'''
Parses a VCF (V4.2) file and builds a list of variants.

Also supports the comparison between a golden var list.
'''
import argparse
import sys
from build_erg import read_var, read_genome

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
    def print(self, chrm=False, pos=True, ref_allele=True, alt_allele=True):
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
            msg += str(self.alt_alleles)
            msg += ' '
        print (msg)

def write_erg_wrt_ref(var_list, ref_genome, f_len):
    '''
    Write one ERG seq
    '''
    erg = ''
    start_pos = var_list[0].v_pos - f_len
    end_v = var_list[len(var_list) - 1]
    
    # only looks at the first allele at a locus
    # TODO
    end_alt_allele = end_v.alt_alleles[0]
    # if len(end_v.alt_alleles) > 1:
    # else:
    #     end_alt_allele = end_v.alt_alleles
    end_pos = end_v.v_pos + f_len + len(end_alt_allele)
    erg_ref = ref_genome[start_pos: end_pos]
    
    prev_bound = [0,0]
    nxt_start_pos = start_pos
    for v in var_list:
        # skips overlapped variants
        if v.v_pos >= prev_bound[0] and \
            v.v_pos < prev_bound[1]:
            continue
        prev_bound[0] = v.v_pos
        prev_bound[1] = v.v_pos + len(v.alt_alleles[0])
        erg += ref_genome[nxt_start_pos : v.v_pos]
        erg += v.alt_alleles[0]
        nxt_start_pos = v.v_pos + len(v.alt_alleles[0])
    if len(erg) == 0:
        erg = erg_ref
    else:
        erg += ref_genome[nxt_start_pos : end_pos]
    
    assert end_pos - start_pos == len(erg)
    WRITE_ERG = True
    if WRITE_ERG:
        CHRM = 9
        print (
            '>%s-erg-%s-%s' % 
            (CHRM, start_pos, end_pos)
        )
        print (erg)
    return len(erg)

def build_erg_wrt_ref(ref_genome, var_list, args):
    '''
    Reads a set of variants under VCF4.2 format and constructs its ERG
    '''
    f_len = args.flanking_len
    tmp_var_list = []
    num_erg = 0
    total_len_erg = 0

    for vinfo in var_list:
        prev_var = 0
        if len(tmp_var_list) > 0:
            prev_var = tmp_var_list[len(tmp_var_list) - 1]
        if len(tmp_var_list) > 0 and \
            vinfo.v_pos > prev_var.v_pos + 2 * f_len:
            # write previous vars
            len_erg = write_erg_wrt_ref(tmp_var_list, ref_genome, f_len)
            num_erg += 1
            total_len_erg += len_erg
            # reset
            tmp_var_list = [vinfo]
        elif len(tmp_var_list) > 0:
            tmp_var_list.append(vinfo)
        else: # len == 0
            tmp_var_list = [vinfo]
    len_erg = write_erg_wrt_ref(tmp_var_list, ref_genome, f_len)
    num_erg += 1
    total_len_erg += len_erg
    SHOW_SUMMARY = True
    if SHOW_SUMMARY:
        sys.stderr.write ('Num ERGs: %d\n' % num_erg)
        sys.stderr.write ('Avg ERG len: %.2f\n' % (float(total_len_erg) / num_erg))

def build_varlist_from_vcf(args):
    '''
    ##fileformat=VCFv4.2
    '''
    target_vcf_fn = args.vcf
    golden_vcf_fn = args.golden_vcf
    min_qual = float(args.min_qual)
    target_vcf_f = open(target_vcf_fn, 'r')
    
    gvar_list = read_var(golden_vcf_fn, remove_conflict=True, remove_homo_alt=False, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    var_list = []

    num_vars = 0
    num_snps = 0
    num_indels = 0
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
        
        var_list.append(vcf)
        
        for i, alt_allele in enumerate(vcf.alt_alleles):
            # a diploid variant is counted twice
            num_vars += 1
            v_type = vcf.v_type.split(',')[i]
            if v_type == 'SNP':
                num_snps += 1
            elif v_type in ['INS', 'DEL']:
                num_indels += 1
            
            # var_list.append([vcf.v_chrom, v_type, vcf.v_pos, vcf.ref_allele, alt_allele])

            for j in range(golden_idx, len(gvar_list)):
                if gvar_list[j].chrm == vcf.v_chrom and \
                gvar_list[j].ref_pos == vcf.v_pos and \
                gvar_list[j].ref_allele == vcf.ref_allele and \
                gvar_list[j].alt_allele == alt_allele:
                    golden_idx = j
                    num_true_pos += 1
                    break
                if gvar_list[j].chrm == vcf.v_chrom and \
                    gvar_list[j].ref_pos > vcf.v_pos:
                    golden_idx = max(j - 1, 0)
                    break
    
    SHOW_STATS = True
    if SHOW_STATS:
        # Stats
        precision = 100 * float(num_true_pos) / num_vars
        sensitivity = 100 * float(num_true_pos) / len(gvar_list)
        
        sys.stderr.write ('NUM GOLDEN VARS: %d\n' % len(gvar_list))
        # sys.stderr.write ('NUM VARS', num_vars)
        sys.stderr.write ('NUM SNPS: %d\n' % num_snps)
        sys.stderr.write ('NUM INDELS: %d\n' % num_indels)
        sys.stderr.write ('NUM TRUE POS: %d\n' % num_true_pos)
        
        sys.stderr.write ('PRECISION = %.2f%%\n' % precision)
        sys.stderr.write ('SENSITIVITY = %.2f%%\n' % sensitivity)

    return var_list

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
        help='golden vcf file ['']'
    )
    parser.add_argument(
        '-r', '--ref_genome',
        default='',
        help='ref genome file ['']'
    )
    parser.add_argument(
        '-f', '--flanking_len', type=int,
        default=100,
        help='the length of flanking regions, total length of an erg is 2*f+1 [100]'
    )
    parser.add_argument(
        '-erg', '--build_erg',
        default=None,
        help='Specify to build ERG based on input variants. The ERG will be output through stdout. `ref_genome` must be specified as well. [None]'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    var_list = build_varlist_from_vcf(args)

    if args.build_erg != None:
        ref_genome = read_genome(args.ref_genome, print_main=True)
        sys.stderr.write ('Warning: ignore HET vars\n')
        build_erg_wrt_ref(ref_genome, var_list, args)
