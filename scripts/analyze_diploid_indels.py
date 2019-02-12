'''
Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse, math
from analyze_sam import SamInfo, parse_line, load_golden_dic, compare_sam_info, Summary
from build_erg import read_var

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n', '--sam',
        help='target sam file'
    )
    parser.add_argument(
        '-g', '--golden',
        help='golden sam file'
    )
    parser.add_argument(
        '-v', '--var',
        help='the file specifying variants'
    )
    parser.add_argument(
        '-t', '--threshold', type=int,
        default=10,
        help='(int) max allowed distance for a correct mapping [10]'
    )
    parser.add_argument(
        '--read_len', type=int,
        default=100,
        help='(int) read length [100]'
    )
    parser.add_argument(
        '-p', '--personalized', type=int,
        default=0,
        help='(int) specify whether the ref seq(s) are standard (0) or personalized-diploid (2) sample [0]'
    )
    parser.add_argument(
        '--step_size', type=int,
        default=1000,
        help='(int) the step size for main/alt offset indexes [1000]'
    )
    parser.add_argument(
        '--write_wrt_correctness', type=int,
        default=None,
        help='(int) If set, writes two files recording correct/incorrect alignments respectively. The output files use target sam prefix [None].'
    )
    parser.add_argument(
        '--write_wrt_mapq', type=int,
        default=None,
        help='(int) If specified, writes two files recording alignments with mapq >= t and mapq < t. This argument is the threshold. The output files use target sam prefix [None].'
    )
    args = parser.parse_args()
    return args

def build_index(var_list, per, MAIN_STRAND, ALT_STRAND):
    '''
    Reads var_list and records variants based on hapA/B coordinates
    '''
    SHOW_CONFLICTS = False
    # Conflicts are not always bad, e.g.,
    # A   9   SNP 121398204   121393883   A   C
    # B   9   SNP 121398204   121393566   A   T
    main_index = {}
    alt_index = {}
    for v in var_list:
        # only main_index is needed for ref-based alignment
        # not checking conflicts here
        if per == 0:
            pos = v.ref_pos
            for i in range(pos, pos + len(v.alt_allele)):
                main_index[i] = [v.alt_pos, v.vtype, v.ref_allele, v.alt_allele]
            continue
        # for personalized experiment
        pos = v.alt_pos
        c_pos = v.ref_pos + v.cor_offset
        if v.strand == MAIN_STRAND:
            if SHOW_CONFLICTS:
                if main_index.get(pos):
                    print (pos, main_index[pos])
                    print (v.line)
                if alt_index.get(c_pos):
                    print (c_pos, alt_index[c_pos])
                    print (v.line)
            for i in range(pos, pos + len(v.alt_allele)):
                main_index[i] = [c_pos, v.vtype, v.ref_allele, v.alt_allele]
            alt_index[c_pos] = [pos, v.vtype, v.ref_allele, v.alt_allele]
        else:
            if SHOW_CONFLICTS:
                if main_index.get(c_pos):
                    print (c_pos, main_index[c_pos])
                    print (v.line)
                if alt_index.get(pos):
                    print (pos, alt_index[pos])
                    print (v.line)
            main_index[c_pos] = [pos, v.vtype, v.ref_allele, v.alt_allele]
            for i in range(pos, pos + len(v.alt_allele)):
                alt_index[i] = [c_pos, v.vtype, v.ref_allele, v.alt_allele]
    return main_index, alt_index

def build_offset_index(var_list, per, MAIN_STRAND, ALT_STRAND):
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
    
    main_index, alt_index = build_index(var_list, per, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
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
    '''
    Uses variable 'dip_flag' to handle different cases of a diploid alignment and check if the alignment is correct.
    '''
    # don't check the other strand
    if dip_flag in ['same_strand']:
        return compare_sam_info(info, g_info, threshold)
    elif dip_flag in ['same_strand_ref']:
        # neglect chrom name difference
        # info.chrm = g_info.chrm
        i_low = int(info.pos / STEP)
        i_high = math.ceil(info.pos / STEP)
        # try hapA
        if i_low >= len(main_offset_index):
            offset_lowA = main_offset_index[len(main_offset_index) - 1]
        else:
            offset_lowA = main_offset_index[i_low]
        if i_high >= len(main_offset_index):
            offset_highA = main_offset_index[len(main_offset_index) - 1]
        else:
            offset_highA = main_offset_index[i_high]
        
        # try hapB
        if i_low >= len(alt_offset_index):
            offset_lowB = alt_offset_index[len(alt_offset_index) - 1]
        else:
            offset_lowB = alt_offset_index[i_low]
        if i_high >= len(alt_offset_index):
            offset_highB = alt_offset_index[len(alt_offset_index) - 1]
        else:
            offset_highB = alt_offset_index[i_high]
        
        comp = compare_sam_info(info, g_info, threshold, [offset_highA, offset_lowA, offset_highB, offset_lowB], ignore_chrm=True)
        if comp == False and __debug__:
            offsets = [offset_lowA, offset_highA, offset_lowB, offset_highB]
            print_near_aln(offsets, info, g_info, 1000)
        return comp
    # check the other strand
    elif dip_flag in ['diff_id', 'diff_var']:
        if info.chrm == MAIN_CHRM:
            # info.chrm = ALT_CHRM
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
            comp = compare_sam_info(info, g_info, threshold, [offset_low, offset_high], ignore_chrm=True)
            if comp == False and __debug__:
                offsets = [offset_low, offset_high]
                print_near_aln(offsets, info, g_info, 1000)
            return comp
        elif info.chrm == ALT_CHRM:
            # info.chrm = MAIN_CHRM
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
            comp = compare_sam_info(info, g_info, threshold, [offset_low, offset_high], ignore_chrm=True)
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

def check_var_in_region(info, main_index, alt_index, MAIN_CHRM, ALT_CHRM, READ_LEN):
    '''
    For an alignment segment, reads indexes to count the number of variants covered by the segment.
    '''
    num_var = 0
    for i in range(info.pos, info.pos + READ_LEN):
        if info.chrm == MAIN_CHRM:
            if main_index.get(i) != None:
                num_var += 1
                # break
        elif info.chrm == ALT_CHRM:
            if alt_index.get(i) != None:
                num_var += 1
                # break
        else:
            print ('Error: unexpected chrm', info.chrm)
            info.print()
            exit()
    return num_var

def analyze_diploid_indels(args):
    '''
    Handles I/O and different opperating modes of this script.
    '''
    sam_fn = args.sam
    golden_fn = args.golden
    threshold = args.threshold
    var_fn = args.var
    personalized = args.personalized
    write_wrt_correctness = args.write_wrt_correctness
    write_wrt_mapq = args.write_wrt_mapq

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

    # diploid personalized ref
    if personalized == 2:
        var_list = read_var(var_fn, remove_conflict=True, remove_coexist=True)
        main_index, alt_index, main_offset_index, alt_offset_index = build_offset_index(var_list, per=2, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # standard ref seq
    elif personalized == 0:
        var_list = read_var(var_fn, remove_conflict=True, remove_coexist=False)
        main_index, alt_index, main_offset_index, alt_offset_index = build_offset_index(var_list, per=0, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    else:
        print ('Error: unsupported personalized parameter', personalized)
        exit()
    
    sam_f = open(sam_fn, 'r')
    sam_prefix = sam_fn[: sam_fn.find('.')]
    golden_dic = load_golden_dic(golden_fn, 1)
    summary = Summary(has_answer=True)
    
    if write_wrt_correctness:
        correct_fn = sam_prefix + '-correct.sam'
        incorrect_fn = sam_prefix + '-incorrect.sam'
        print ('Write sam files %s and %s wrt to correctness...' % (correct_fn, incorrect_fn))
        correct_f = open(correct_fn, 'w')
        incorrect_f = open(incorrect_fn, 'w')
    
    if write_wrt_mapq:
        highmapq_fn = sam_prefix + '-mapqgeq' + str(write_wrt_mapq) + '.sam'
        lowmapq_fn = sam_prefix + '-mapql' + str(write_wrt_mapq) + '.sam'
        print ('Write sam files %s and %s wrt to mapq...' % (highmapq_fn, lowmapq_fn))
        highmapq_f = open(highmapq_fn, 'w')
        lowmapq_f = open(lowmapq_fn, 'w')
        for line in sam_f:
            name, info = parse_line(line)
            # headers
            if name == 'header':
                highmapq_f.write(line)
                lowmapq_f.write(line)
                continue
            if info.mapq >= write_wrt_mapq:
                highmapq_f.write(line)
            else:
                lowmapq_f.write(line)
            continue
        exit()
   
    CHECK_VAR_OVERLAPPING_REF = False

    for line in sam_f:
        if personalized == 0:
            # name, info = parse_line(line)
            name, info = parse_line(line, erg=True)
        elif personalized == 2:
            name, info = parse_line(line, erg=True)
        # headers
        if name == 'header':
            continue
        name_chrm_mismatch = (name.find(MAIN_HAP) > 0 and info.chrm != MAIN_CHRM) or (name.find(ALT_HAP) > 0 and info.chrm != ALT_CHRM)
        summary.add_one()
        if info.is_unaligned():
            summary.add_unaligned()
            comp = False
        # aligned to incorrect haplotype
        elif personalized == 2 and name_chrm_mismatch:
            num_var = check_var_in_region(info, main_index, alt_index,  MAIN_CHRM=MAIN_CHRM, ALT_CHRM=ALT_CHRM, READ_LEN=READ_LEN)
            # aligned to incorrect haplotype and two haps are equal
            if num_var == 0:
                comp = diploid_compare(info, golden_dic[name], threshold, 'diff_id', main_offset_index, alt_offset_index)
                summary.add_diff_id(comp)
            # aligned to incorrect haplotype and two haps are NOT equal
            else:
                comp = diploid_compare(info, golden_dic[name], threshold, 'diff_var', main_offset_index, alt_offset_index)
                summary.add_diff_var(comp)
        else:
            if personalized == 2:
                comp = diploid_compare(info, golden_dic[name], threshold, 'same_strand')
                num_var = check_var_in_region(info, main_index, alt_index,  MAIN_CHRM=MAIN_CHRM,ALT_CHRM=ALT_CHRM, READ_LEN=READ_LEN)
                # aligned to correct haplotype and two haps are equal
                if num_var == 0:
                    summary.add_same_id(comp)
                # aligned to correct haplotype and two haps are NOT equal
                else:
                    summary.add_same_var(comp)
            elif personalized == 0:
                if CHECK_VAR_OVERLAPPING_REF:
                    num_var = check_var_in_region(info, main_index, alt_index, MAIN_CHRM=MAIN_CHRM, ALT_CHRM=ALT_CHRM, READ_LEN=READ_LEN)
                else:
                    num_var = 0
                comp = diploid_compare(info, golden_dic[name], threshold, 'same_strand_ref', main_offset_index, alt_offset_index)
                
                # simply add results to the same_strand category
                # summary.add_same_strand(comp)

                # add results to same-diff and same-var, this doesn't actually
                # matter for ref-based alignemnt, but helps analysis
                # aligned to correct haplotype and two haps are equal
                if num_var == 0:
                    summary.add_same_id(comp)
                # aligned to correct haplotype and two haps are NOT equal
                else:
                    summary.add_same_var(comp)

        if write_wrt_correctness:
            if comp:
                correct_f.write(line)
            else:
                incorrect_f.write(line)
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
    analyze_diploid_indels(args)
