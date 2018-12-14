import argparse
from analyze_multimapped_regions import VarInfo, build_var_dic

MAIN_STRAND = 'A'
ALT_STRAND = 'B'

def read_genome(fn):
    with open(fn, 'r') as f:
        # seq[0] is empty to fit vcf coordinate (starting from 1)
        seq = '^'
        for line in f:
            if line.startswith('>') == False:
                seq += line[:line.find('\\')]
    return seq
'''
def write_erg(var_list, genome, start_p, f_len, main_var_list):
    print ('write_erg (%d vars)' % len(var_list))
    print (var_list[0].ref_pos)
    for i, v in enumerate(var_list):
        main_allele_len = len(v.ref_allele)
        for j, vm in enumerate(main_var_list):
            if v.samepos_diffvar(vm):
                main_allele_len = len(vm.alt_allele)
        first_end_p = start_p + f_len
        second_start_p = first_end_p + main_allele_len
        second_end_p = second_start_p + f_len
        erg = genome[start_p: first_end_p]
        erg += v.alt_allele
        erg += genome[second_start_p: second_end_p]
        print (erg)
    input ()
'''
def write_erg(var_list, main_genome, f_len):
    '''
    Write one ERG seq
    '''
    erg = ''
    main_pos_list = []
    erg_start_pos = 0
    offset = 0
    for i, v in enumerate(var_list):
        if v.strand == MAIN_STRAND:
            # v: ref hapA
            main_pos = v.alt_pos
            main_pos_list.append(main_pos)
            erg_alt_allele = v.ref_allele
            offset = -v.offset + v.cor_offset
        else:
            # v: ref hapB
            main_pos = v.alt_pos - v.offset + v.cor_offset
            main_pos_list.append(main_pos)
            erg_alt_allele = v.alt_allele
            offset = v.offset - v.cor_offset
#        print (main_pos, erg_alt_allele)
        if i == 0:
            erg_start_pos = main_pos - f_len
            erg += main_genome[main_pos - f_len : main_pos]
        else:
            erg += main_genome[main_pos_list[i-1] + 1 : main_pos]
        erg += erg_alt_allele
    erg += main_genome[main_pos_list[-1] + 1 : main_pos_list[-1] + 1 + f_len]
    print (
        '>hapB-%s-%s' % 
        (erg_start_pos + offset, main_pos_list[-1] + 1 + f_len + offset)
    )
    print (erg)


def build_erg(args):
    mg_fn = args.main_genome
    var_fn = args.var
    hap_mode = args.hap
    f_len = args.flanking_len
    # In var file:
    # the strand occurs earlier is the main_strand

    main_genome = read_genome(mg_fn)

    with open(var_fn, 'r') as var_f:
        var_list = []
        for line in var_f:
            vinfo = VarInfo(line)
            if len(var_list) > 0:
                prev_var = var_list[len(var_list) - 1]
            else:
                prev_var = 0
            if __debug__:
                print (prev_var)
                print (vinfo.line)
            if vinfo.strand == MAIN_STRAND:
                if len(var_list) > 0 and \
                    vinfo.alt_pos > prev_var.alt_pos + 2 * f_len:
                    # write previous
                    write_erg(var_list, main_genome, f_len)
                    # reset
                    var_list = [vinfo]
                elif len(var_list) > 0:
                    var_list.append(vinfo)
                else: #len == 0
                    var_list = [vinfo]
                continue
            # var not main strand
            if vinfo.samevar(prev_var):
                # delete if a var occurs in both strands
                if __debug__:
                    print ('samevar')
                    print ('A*', prev_var.line)
                    print ('B*', vinfo.line)
                del var_list[len(var_list) - 1]
                continue
            else:
                if prev_var == 0:
                    var_list.append(vinfo)
                elif vinfo.alt_pos > prev_var.alt_pos + 2 * f_len:
                    write_erg(var_list, main_genome, f_len)
                    var_list = [vinfo]
                else:
                    var_list.append(vinfo)
        write_erg(var_list, main_genome, f_len)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m', '--main_genome',
        help='main genome file'
    )
    parser.add_argument(
        '-v', '--var',
        help='var file'
    )
    parser.add_argument(
        '--hap', type=int,
        default=1,
        help='set 1 to enable haplotype-based erg [1]'
    )
    parser.add_argument(
        '-f', '--flanking_len', type=int,
        default=50,
        help='the length of flanking regions, total length of an erg is 2*f+1 [50]'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    build_erg(args)
