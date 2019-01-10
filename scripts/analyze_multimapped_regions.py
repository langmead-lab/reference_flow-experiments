'''
Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse
from analyze_sam import SamInfo, parse_line, load_golden_dic, compare_sam_info, Summary
from build_erg import build_erg, read_var

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
        '-p', '--personalized', type=int,
        default=0,
        help='specify whether the ref seq(s) are standard (0), personalized-haploid (1), or personalized-diploid (2) sample [0]'
    )
    parser.add_argument(
        '--var',
        default=None,
        help='the file specifying the variants'
    )
    args = parser.parse_args()

    # Global variables
    global STEP, MAIN_CHRM, ALT_CHRM, MAIN_HAP, ALT_HAP, MAIN_STRAND, ALT_STRAND
    MAIN_HAP = 'hapA'
    ALT_HAP = 'hapB'
    MAIN_STRAND = 'A'
    ALT_STRAND = 'B'
    STEP = 1000
    if args.personalized == 2:
        MAIN_CHRM = '9A'
        ALT_CHRM = '9B'
    else:
        MAIN_CHRM = '9'
        ALT_CHRM = '9'
    return args

def build_offset_index(var_list, per):
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
        if v.strand == MAIN_STRAND:
            main_pos = v.alt_pos
            alt_pos = v.alt_pos + v.offset - v.cor_offset
            if per == 2:
                main_offset = -v.offset + v.cor_offset
                alt_offset = v.offset - v.cor_offset
            else:
                # offset: ref to hap
                main_offset = v.offset
                alt_offset = v.cor_offset
        elif v.strand == ALT_STRAND:
            main_pos = v.alt_pos + v.offset - v.cor_offset
            alt_pos = v.alt_pos
            if per == 2:
                main_offset = v.offset - v.cor_offset
                alt_offset = -v.offset + v.cor_offset
            else:
                # offset: ref to hap
                main_offset = v.cor_offset
                alt_offset = v.offset
        else:
            print ('Error: unspecified strand', v.strand)
            exit()
        while main_pos > STEP * len(main_offset_index):
            if main_pos > STEP * (len(main_offset_index) - 1):
                main_offset_index.append(main_offset)
            else:
                main_offset_index.append(main_offset_index[len(main_offset_index) - 1])
        while alt_pos > STEP * len(alt_offset_index):
            if alt_pos > STEP * (len(alt_offset_index) - 1):
                alt_offset_index.append(alt_offset)
            else:
                alt_offset_index.append(alt_offset_index[len(alt_offset_index) - 1])
        if __debug__:
            print (v.line)
            print (main_offset_index)
            print (alt_offset_index)
            input ()
    
    if per == 2:
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
    else:
        main_index = {}
        alt_index = {}
    return main_index, alt_index, main_offset_index, alt_offset_index

def diploid_compare(
    info, 
    g_info, 
    threshold, 
    dip_flag, 
    main_offset_index = {}, 
    alt_offset_index = {}
):
    # don't check the other strand
    if dip_flag in ['same_strand']:
        return compare_sam_info(info, g_info, threshold)
    elif dip_flag in ['same_strand_ref']:
        # neglect chrom name difference
        info.chrm = g_info.chrm
        i = int(info.pos / STEP)
        # try hapA
        if i >= len(main_offset_index):
            diff = main_offset_index[len(main_offset_index) - 1]
        else:
            diff = main_offset_index[i]
        info.pos += diff
        # try hapB
        comp1 = compare_sam_info(info, g_info, threshold)
        info.pos -= diff
        if i >= len(alt_offset_index):
            diff = alt_offset_index[len(alt_offset_index) - 1]
        else:
            diff = alt_offset_index[i]
        info.pos += diff
        comp2 = compare_sam_info(info, g_info, threshold)
        comp = comp1 | comp2
        return comp
    # check the other strand
    elif dip_flag in ['diff_id', 'diff_var']:
        if info.chrm == MAIN_CHRM:
            info.chrm = ALT_CHRM
            i = int(info.pos / STEP)
            if i >= len(main_offset_index):
                info.pos += main_offset_index[len(main_offset_index) - 1]
            else:
                info.pos += main_offset_index[i]
            comp = compare_sam_info(info, g_info, threshold)
            return comp
        elif info.chrm == ALT_CHRM:
            info.chrm = MAIN_CHRM
            i = int(info.pos / STEP)
            if i >= len(alt_offset_index):
                info.pos += alt_offset_index[len(alt_offset_index) - 1]
            else:
                info.pos += alt_offset_index[i]
            comp = compare_sam_info(info, g_info, threshold)
            return comp
        else:
            print ('Error: invalid chrm', info.chrm)
            exit()
    else:
        print ('Error: undistinguished dip_flag: %s' % dip_flag)
        return False

def analyze_mutimapped_regions(args):
    sam_fn = args.sam
    golden_fn = args.golden
    threshold = args.threshold
    read_len = args.read_len
    var_fn = args.var
    personalized = args.personalized

    var_list = read_var(var_fn, remove_redundant=True)
    # diploid personalized ref
    if personalized == 2:
        main_index, alt_index, main_offset_index, alt_offset_index = build_offset_index(var_list, per=2)
    # standard ref seq
    elif personalized == 0:
        main_index, alt_index, main_offset_index, alt_offset_index = build_offset_index(var_list, per=0)
    else:
        print ('Error: unsupported personalzed parameter', personalized)
        exit()
    '''
    MAIN/ALT indexes:
        key: pos on MAIN/ALT
        value: pos on ALT/MAIN, VTYPE, ALLELES
    '''
    golden_dic = load_golden_dic(golden_fn, 1)
    summary = Summary(has_answer=True)
    sam_f = open(sam_fn, 'r')
    for line in sam_f:
        name, info = parse_line(line, 0)
        # headers
        if name == 'header':
            '''
            # kept for low-q experiment
            print (line[:line.find('\\')])
            '''
            continue
        '''
        # kept for low-q experiment
        if info.mapq < 10:
           print (line[:line.find('\\')])
        continue
        '''
        summary.add_one()
        # aligned to incorrect haplotype
        if info.is_unaligned():
            summary.add_unaligned()
            comp = False
        elif (name.find(MAIN_HAP) > 0 and info.chrm != MAIN_CHRM) \
        or (name.find(ALT_HAP) > 0 and info.chrm != ALT_CHRM):
            v_id_hap = True
            for i in range(info.pos, info.pos + read_len):
                if info.chrm == MAIN_CHRM:
                    if main_index.get(i) != None:
                        v_id_hap = False
                        break
                elif info.chrm == ALT_CHRM:
                    if alt_index.get(i) != None:
                        v_id_hap = False
                        break
                else:
                    print ('Error: unexpected chrm', info.chrm)
                    exit()
            # aligned to incorrect haplotype and two haps are NOT equal
            if v_id_hap == False:
                comp = diploid_compare(info, golden_dic[name], threshold, 'diff_var', main_offset_index, alt_offset_index)
                summary.add_diff_var(comp)
            else:
                comp = diploid_compare(info, golden_dic[name], threshold, 'diff_id', main_offset_index, alt_offset_index)
                summary.add_diff_id(comp)
        else:
            if personalized == 2:
                comp = diploid_compare(info, golden_dic[name], threshold, 'same_strand')
            elif personalized == 0:
                comp = diploid_compare(info, golden_dic[name], threshold, 'same_strand_ref', main_offset_index, alt_offset_index)
            else:
                print ('Error: unsupported personalzed parameter', personalized)
                exit()
            summary.add_same_strand(comp)
        if __debug__:
            if comp: continue
            print ('info')
            info.print(flag=False, mapq=False, score=False)
            print ()
            print ('golden')
            golden_dic[name].print(flag=False, mapq=False, score=False)
            input ('comp=%s\n' % comp)

    summary.show_summary(has_answer = True)
    sam_f.close()

if __name__ == '__main__':
    args = parse_args()
    analyze_mutimapped_regions(args)
