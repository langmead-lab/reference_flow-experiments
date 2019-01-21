'''
Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse, math
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
        '--var',
        help='the file specifying the variants'
    )
    parser.add_argument(
        '-t', '--threshold', type=int,
        default=10,
        help='max allowed distance for a correct mapping [10]'
    )
    parser.add_argument(
        '--read_len', type=int,
        default=100,
        help='read length [100]'
    )
    parser.add_argument(
        '-p', '--personalized', type=int,
        default=0,
        help='specify whether the ref seq(s) are standard (0), personalized-haploid (1), or personalized-diploid (2) sample [0]'
    )
    parser.add_argument(
        '--step_size', type=int,
        default=1000,
        help='the step size for main/alt offset indexes [1000]'
    )
    args = parser.parse_args()

    # Global variables
    global STEP, MAIN_CHRM, ALT_CHRM, MAIN_HAP, ALT_HAP, MAIN_STRAND, ALT_STRAND, READ_LEN
    MAIN_HAP = 'hapA'
    ALT_HAP = 'hapB'
    MAIN_STRAND = 'A'
    ALT_STRAND = 'B'
    STEP = args.step_size
    READ_LEN = args.read_len
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
    SHOW_BUILD_INFO = False
    SHOW_DUP_WARN = False
    if SHOW_BUILD_INFO and __debug__:
        print ('DEBUG_INFO: build_offset_index ')
    tmp_v = 0
    for v in var_list:
        if v.strand == MAIN_STRAND:
            main_pos = v.alt_pos
            alt_pos = v.ref_pos + v.cor_offset
            if per == 2:
                main_offset = -v.offset + v.cor_offset
                alt_offset = v.offset - v.cor_offset
            else:
                # offset: ref to hap
                main_pos = v.ref_pos
                alt_pos = v.ref_pos
                main_offset = v.offset
                alt_offset = v.cor_offset
            if v.ref_pos == tmp_v:
                if SHOW_DUP_WARN:
                    print ('Warning: duplicated variant', v.line)
            tmp_v = v.ref_pos
        elif v.strand == ALT_STRAND:
            main_pos = v.ref_pos + v.cor_offset
            alt_pos = v.alt_pos
            if per == 2:
                main_offset = v.offset - v.cor_offset
                alt_offset = -v.offset + v.cor_offset
            else:
                # offset: ref to hap
                main_pos = v.ref_pos
                alt_pos = v.ref_pos
                main_offset = v.cor_offset
                alt_offset = v.offset
        else:
            print ('Error: unspecified strand', v.strand)
            exit()
        
        i_main = math.ceil(main_pos / STEP)
        while i_main >= len(main_offset_index):
            main_offset_index.append(main_offset_index[len(main_offset_index) - 1])
        main_offset_index[i_main] = main_offset
        i_alt = math.ceil(alt_pos / STEP)
        while i_alt >= len(alt_offset_index):
            alt_offset_index.append(alt_offset_index[len(alt_offset_index) - 1])
        alt_offset_index[i_alt] = alt_offset

        if __debug__ and SHOW_BUILD_INFO:
            print (v.line)
            print (len(main_offset_index), main_offset_index)
            print (len(alt_offset_index), alt_offset_index)
            input ()
    
    if per == 2:
        main_index, alt_index = \
            build_erg(
                main_genome = '', 
                ref_genome = '',
                test_genome = '',
                var_list = var_list,
                hap_mode = 1, 
                f_len = READ_LEN, 
                mode = 'index'
            )
    else:
        main_index = {}
        alt_index = {}
    return main_index, alt_index, main_offset_index, alt_offset_index

def print_near_aln(offsets, info, g_info, threshold):
    tmp = []
    for i in offsets:
        tmp.append(abs(info.pos + i - g_info.pos))
    diff = min(tmp)
    if diff < threshold:
        print ('offsets', offsets)
        print ('diff', diff)
        print ('info')
        info.print(flag=False, mapq=False, score=False)
        print ()
        print ('golden')
        g_info.print(flag=False, mapq=False, score=False)
        input()

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
        i_low = int(info.pos / STEP)
        i_high = math.ceil(info.pos / STEP)
        # try hapA
        if i_low >= len(main_offset_index):
            offset_lowA = main_offset_index[len(main_offset_index) - 1]
        else:
            offset_lowA = main_offset_index[i_low]
        info.pos += offset_lowA
        comp1 = compare_sam_info(info, g_info, threshold)
        info.pos -= offset_lowA
        if i_high >= len(main_offset_index):
            offset_highA = main_offset_index[len(main_offset_index) - 1]
        else:
            offset_highA = main_offset_index[i_high]
        info.pos += offset_highA
        comp1 = comp1 | compare_sam_info(info, g_info, threshold)
        info.pos -= offset_highA
        # try hapB
        if i_low >= len(alt_offset_index):
            offset_lowB = alt_offset_index[len(alt_offset_index) - 1]
        else:
            offset_lowB = alt_offset_index[i_low]
        info.pos += offset_lowB
        comp2 = compare_sam_info(info, g_info, threshold)
        info.pos -= offset_lowB
        if i_high >= len(alt_offset_index):
            offset_highB = alt_offset_index[len(alt_offset_index) - 1]
        else:
            offset_highB = alt_offset_index[i_high]
        info.pos += offset_highB
        comp2 = comp2 | compare_sam_info(info, g_info, threshold)
        info.pos -= offset_highB

        comp = comp1 | comp2
        
        if comp == False and __debug__:
            offsets = [offset_lowA, offset_highA, offset_lowB, offset_highB]
            print_near_aln(offsets, info, g_info, 1000)
        
        return comp
    # check the other strand
    elif dip_flag in ['diff_id', 'diff_var']:
        if info.chrm == MAIN_CHRM:
            info.chrm = ALT_CHRM
            i_low = int(info.pos / STEP)
            i_high = math.ceil(info.pos / STEP)
            if i_low >= len(main_offset_index):
                offset_low = main_offset_index[len(main_offset_index) - 1]
            else:
                offset_low = main_offset_index[i_low]
            if i_high >= len(main_offset_index):
                offset_high = main_offset_index[len(main_offset_index) - 1]
            else:
                offset_high = main_offset_index[i_high]
            comp = compare_sam_info(info, g_info, threshold, offset_low) | compare_sam_info(info, g_info, threshold, offset_high)
            if comp == False and __debug__:
                offsets = [offset_low, offset_high]
                print_near_aln(offsets, info, g_info, 1000)
            return comp
        elif info.chrm == ALT_CHRM:
            info.chrm = MAIN_CHRM
            i_low = int(info.pos / STEP)
            i_high = math.ceil(info.pos / STEP)
            if i_low >= len(alt_offset_index):
                offset_low = alt_offset_index[len(alt_offset_index) - 1]
            else:
                offset_low = alt_offset_index[i_low]
            if i_high >= len(alt_offset_index):
                offset_high = alt_offset_index[len(alt_offset_index) - 1]
            else:
                offset_high = alt_offset_index[i_high]
            comp = compare_sam_info(info, g_info, threshold, offset_low) | compare_sam_info(info, g_info, threshold, offset_high)
            if comp == False and __debug__:
                offsets = [offset_low, offset_high]
                print_near_aln(offsets, info, g_info, 1000)
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
    golden_dic = load_golden_dic(golden_fn, 1)
    summary = Summary(has_answer=True)
    PERFORM_LOWQ_EXP = False
    sam_f = open(sam_fn, 'r')
    for line in sam_f:
        name, info = parse_line(line, 0)
        # headers
        if name == 'header':
            if PERFORM_LOWQ_EXP:
                # kept for low-q experiment
                print (line[:line.find('\\')])
            continue
        if PERFORM_LOWQ_EXP:
            # kept for low-q experiment
            if info.mapq < 10:
                print (line[:line.find('\\')])
            continue
        summary.add_one()
        # aligned to incorrect haplotype
        if info.is_unaligned():
            summary.add_unaligned()
            comp = False
        elif (name.find(MAIN_HAP) > 0 and info.chrm != MAIN_CHRM) \
        or (name.find(ALT_HAP) > 0 and info.chrm != ALT_CHRM):
            v_id_hap = True
            for i in range(info.pos, info.pos + READ_LEN):
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
        if __debug__ and comp == False:
            print (name)
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
