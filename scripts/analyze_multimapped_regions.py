'''
Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse
from analyze_sam import SamInfo, parse_line, load_golden_dic, compare_sam_info, Summary
from build_erg import build_erg, read_var

MAIN_CHRM = '9A'
MAIN_HAP = 'hapA'
ALT_CHRM = '9B'
ALT_HAP = 'hapB'
STEP = 1000

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n', '--sam',
        help='aligned sam file'
    )
    parser.add_argument(
        '-g', '--golden',
        help='golden sam file'
    )
    parser.add_argument(
        '-t', '--threshold', type=int,
        default=10,
        help='max allowed distance for a correct mapping [10]'
    )
    parser.add_argument(
        '-r', '--read_len', type=int,
        default=100,
        help='read length [100]'
    )
    parser.add_argument(
        '-d', '--diff',
        help='the snp file specifying the differences between two haplotypes'
    )
    parser.add_argument(
        '--var',
        default=None,
        help='the file specifying the variants'
    )
    args = parser.parse_args()
    return args

# def build_var_dic(var_fn):
#     '''
#     Build a dictionary for the .var file
    
#     var file format:
#     STRAND CHRM TYPE REFPOS ALTPOS REF ALT OFFSET

#     key:
#         STRAND_CHRM_ALTPOS (str)
#     value:
#         TYPE, REFPOS, REF, ALT, OFFSET (list)
#     '''
#     var_dic = {}
#     with open(var_fn, 'r') as var_f:
#         for line in var_f:
#             line = line.split()
#             strand = line[0]
#             chrm = line[1]
#             vtype = line[2]
#             ref_pos = line[3]
#             alt_pos = line[4]
#             ref_allele = line[5]
#             alt_allele = line[6]
#             offset = line[7]
#             c_offset = line[8]
#             v_key = strand + '_' + chrm + '_' + alt_pos
#             var_dic[v_key] = vtype, ref_pos, ref_allele, alt_allele, offset, c_offset
#             # input (var_dic)
#     return var_dic

def build_offset_index(var_list):
    '''
    MAIN/ALT-offset indexes are dictionaries with
        key: pos on MAIN/ALT
        value: pos on ALT/MAIN
    
    MAIN/ALT indexes are dictionaries storing
    variants based on MAIN/ALT coordinates
    '''
    # dict storing the diff from main to alt
    # main_pos + main_offset_index[i] = alt_pos
    main_offset_index = [0]
    # dict storing the diff from alt to main
    # alt_pos + alt_offset_index[i] = main_pos
    alt_offset_index = [0]
    for v in var_list:
        if v.chrm == MAIN_CHRM:
            main_pos = v.alt_pos
            alt_pos = v.alt_pos + v.offset - v.cor_offset
            main_offset = -v.offset + v.cor_offset
            alt_offset = v.offset - v.cor_offset
        else:
            main_pos = v.alt_pos + v.offset - v.cor_offset
            alt_pos = v.alt_pos
            main_offset = v.offset - v.cor_offset
            alt_offset = -v.offset + v.cor_offset
        while main_pos > STEP * len(main_offset_index):
            # print ('MAIN_POS')
            # print (main_pos)
            if main_pos > STEP * (len(main_offset_index) - 1):
                main_offset_index.append(alt_offset)
            else:
                main_offset_index.append(main_offset_index[len(main_offset_index) - 1])
        while alt_pos > STEP * len(alt_offset_index):
            # print ('ALT_POS')
            # print (alt_pos)
            if alt_pos > STEP * (len(alt_offset_index) - 1):
                alt_offset_index.append(main_offset)
            else:
                alt_offset_index.append(alt_offset_index[len(alt_offset_index) - 1])
        # print (v.line)
        # print ('len(main_index)', len(main_index))
        # print ('main', main_index)
        # print ('len(alt_index)', len(alt_index))
        # print ('alt ', alt_index)
        # input ()
    main_index, alt_index = \
        build_erg(
            main_genome = '', 
            ref_genome = '',
            test_genome = '',
            var_list = var_list,
            hap_mode = 1, 
            f_len = 100, 
            mode = 'index'
        )
    # print (len(main_index))
    # print (len(alt_index))
    return main_index, alt_index, main_offset_index, alt_offset_index
'''
def build_diff_dic(diff_snp_fn):
    # build dictionary
    diff_snp_dic = {}
    with open(diff_snp_fn, 'r') as diff_snp_f:
        for line in diff_snp_f:
            if line.startswith('pos'):
                continue
            line = line.split()
            pos = int(line[0])
            diff_snp_dic[pos] = line[1], line[2]
    return diff_snp_dic
'''
def diploid_compare(
    info, 
    g_info, 
    threshold, 
    dip_flag, 
    main_offset_index = {}, 
    alt_offset_index = {}
):
    # don't check the other strand
    if dip_flag in ['c', 'n']:
        return compare_sam_info(info, g_info, threshold)
    # check the other strand
    elif dip_flag == 'i':
        if info.chrm == MAIN_CHRM:
            info.chrm = ALT_CHRM
            i = int(info.pos / STEP)
            if i >= len(main_offset_index):
                info.pos += main_offset_index[len(main_offset_index) - 1]
            else:
                info.pos += main_offset_index[i]
            return compare_sam_info(info, g_info, threshold)
        elif info.chrm == ALT_CHRM:
            info.chrm = MAIN_CHRM
            i = int(info.pos / STEP)
            if i >= len(alt_offset_index):
                info.pos += alt_offset_index[len(alt_offset_index) - 1]
            else:
                info.pos += alt_offset_index[i]
            return compare_sam_info(info, g_info, threshold)
        else:
            print ('Error: invalid chrm', info.chrm)
            exit()
        # info.chrm = MAIN_CHRM
        # comp1 = compare_sam_info(info, g_info, threshold)
        # info.chrm = ALT_CHRM
        # comp2 = compare_sam_info(info, g_info, threshold)
#        input ('i: %s' % (comp1 | comp2))
        # return comp1 | comp2
    else:
        print ('Error: undistinguished dip_flag: %s' % dip_flag)
        return False

def show_info(name, info, g_info, dip_flag):
    if 1: #__debug__:
        print (name)
        print (dip_flag)
        print ('info')
        info.print(flag=False, mapq=False, score=False)
        print (' ')
        print ('golden')
        g_info.print(flag=False, mapq=False, score=False)
        input ()

def analyze_mutimapped_regions(args):
    sam_fn = args.sam
    golden_fn = args.golden
    threshold = args.threshold
    read_len = args.read_len
    var_fn = args.var

    var_list = read_var(var_fn, remove_redundant=True)
    # print (var_list)
    # var_dic = build_var_dic(var_fn)
    main_index, alt_index, main_offset_index, alt_offset_index = build_offset_index(var_list)
    '''
    MAIN/ALT indexes:
        key: pos on MAIN/ALT
        value: pos on ALT/MAIN, VTYPE, ALLELES
    '''
    golden_dic = load_golden_dic(golden_fn, 1)
    summary = Summary(has_answer=True)
    with open(sam_fn, 'r') as sam_f:
        for line in sam_f:
            name, info = parse_line(line, 0)
            # headers
            if name == 'header':
#                print (line[:line.find('\\')])
                continue
#            if info.mapq < 10:
#                print (line[:line.find('\\')])
#            continue
            summary.add_one()
            # aligned to incorrect haplotype
            if info.is_unaligned():
                summary.add_unaligned()
                dip_flag = 'u'
                comp = False
            # TODO: this condition should be more generic
            elif (name.find(MAIN_HAP) > 0 and info.chrm != MAIN_CHRM) \
            or (name.find(ALT_HAP) > 0 and info.chrm != ALT_CHRM):
                v_id_hap = True
                for i in range(info.pos, info.pos + read_len):
                    if info.chrm == MAIN_CHRM:
                        if main_index.get(i) != None:
                            v_id_hap = False
                            # print ('var_dic', i, main_index.get(i))
                            # input ()
                            break
                    elif info.chrm == ALT_CHRM:
                        if alt_index.get(i) != None:
                            v_id_hap = False
                            # print (i, alt_index.get(i))
                            # input ()
                            break
                    else:
                        print ('Error: unexpected chrm', info.chrm)
                        exit()
                # aligned to incorrect haplotype and two haps are NOT equal
                if v_id_hap == False:
                    comp = diploid_compare(info, golden_dic[name], threshold, 'n')
                    dip_flag = 'n'
                    summary.add_diff_var(comp)
                else:
                    dip_flag = 'i'
                    comp = diploid_compare(info, golden_dic[name], threshold, 'i', main_offset_index, alt_offset_index)
                    summary.add_diff_id(comp)
            else:
                dip_flag = 'c'
                comp = diploid_compare(info, golden_dic[name], threshold, 'c')
                summary.add_same_strand(comp)
            if __debug__:
                if comp: continue
                print ('dip_flag =', dip_flag)
                print ('info')
                info.print(flag=False, mapq=False, score=False)
                print ()
                print ('golden')
                golden_dic[name].print(flag=False, mapq=False, score=False)
                input ('comp=%s\n' % comp)
    
    summary.show_summary(has_answer = True)

if __name__ == '__main__':
    args = parse_args()
    analyze_mutimapped_regions(args)
