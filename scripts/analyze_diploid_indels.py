'''
Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse, math, sys
from analyze_sam import SamInfo, parse_line, load_golden_dic, compare_sam_info, Summary
from build_erg import read_var, read_genome

#: TODO LEV should support different scoring schemes
from get_levenshtein import levenshtein

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

def build_index(var_list, MAIN_STRAND, ALT_STRAND):
    '''
    Reads var_list and maps variants from hapA/B to the reference coordinate
    
    SHOW_MULT:
        A MULT here can be an ALT locus, such as 
            A   9   INDEL   10362       10362   C  CT
            B   9   INDEL   10362       10362   C  CT
        or a multi-allelic locus, such as
            A   9   SNP 121398204   121393883   A   C
            B   9   SNP 121398204   121393566   A   T
    '''
    SHOW_MULT = False
    main_index = {}
    alt_index = {}
    for v in var_list:
        pos = v.alt_pos
        c_pos = v.ref_pos + v.cor_offset
        if v.strand == MAIN_STRAND:
            if SHOW_MULT:
                if main_index.get(pos):
                    print (pos, main_index[pos])
                    print (v.line)
                if alt_index.get(c_pos):
                    print (c_pos, alt_index[c_pos])
                    print (v.line)
            for i in range(pos, pos + len(v.alt_allele)):
                main_index[i] = [v.ref_pos, v.vtype, v.ref_allele, v.alt_allele]
        else:
            if SHOW_MULT:
                if main_index.get(c_pos):
                    print (c_pos, main_index[c_pos])
                    print (v.line)
                if alt_index.get(pos):
                    print (pos, alt_index[pos])
                    print (v.line)
            for i in range(pos, pos + len(v.alt_allele)):
                alt_index[i] = [v.ref_pos, v.vtype, v.ref_allele, v.alt_allele]
    return main_index, alt_index

def build_offset_index(var_list, per, MAIN_STRAND, ALT_STRAND):
    '''
    MAIN/ALT-offset indexes are dictionaries with
        key: pos on MAIN/ALT
        value: pos on ALT/MAIN
    
    MAIN/ALT indexes are dictionaries storing
    variants based on MAIN/ALT coordinates
    '''
    #: dict storing the diff from main to alt
    # main_pos + main_offset_index[i] = alt_pos
    main_offset_index = [0]
    #: dict storing the diff from alt to main
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
    
    return main_offset_index, alt_offset_index

def build_offset_index_ref(var_list, MAIN_STRAND, ALT_STRAND):
    '''
    ALT1/ALT2-offset indexes are dictionaries with
        key: pos on ALT1/ALT2
        value: offset on ALT1/ALT2 reference sequence at 'key'

    Outputs:
        alt1_offset_index:
            dict storing the diff from alt1 to ref
            alt1_pos - alt1_offset_index[i] = ref_pos
        alt2_offset_index:
            dict storing the diff from alt2 to ref
            alt2_pos - alt2_offset_index[i] = ref_pos
    '''
    alt1_offset_index = [0]
    alt2_offset_index = [0]
    for v in var_list:
        #: offset: ref to hap
        offset = v.offset
        idx = math.ceil(v.alt_pos / STEP)
        if v.strand == MAIN_STRAND:
            while idx >= len(alt1_offset_index):
                alt1_offset_index.append(alt1_offset_index[len(alt1_offset_index) - 1])
            alt1_offset_index[idx] = offset
        elif v.strand == ALT_STRAND:
            while idx >= len(alt2_offset_index):
                alt2_offset_index.append(alt2_offset_index[len(alt2_offset_index) - 1])
            alt2_offset_index[idx] = offset
        else:
            print ('Error: unspecified strand', v.strand)
            exit()
    
    return alt1_offset_index, alt2_offset_index

def print_and_stop(name, offsets, diff, info, g_info):
    '''
    This is for debugging.
    Prints alignment info and stops the script.
    '''
    print ('name', name)
    print ('offsets', offsets)
    print ('diff', diff)
    print ('info')
    info.print(flag=False, mapq=False, score=False)
    print ()
    print ('golden')
    g_info.print(flag=False, mapq=False, score=False)
    input()

def print_near_aln(name, offsets, info, g_info, threshold):
    '''
    Compares alignment with the golden profile if they are near.
    If COMPARE_SEQ is specified, retrieves sequences from ref and haps and calculate the distance.
    '''
    tmp = []
    for i in offsets:
        tmp.append(abs(info.pos + i - g_info.pos))
    diff = min(tmp)
    if (diff < threshold) or (threshold < 0):
        if COMPARE_SEQ:
            global TOTALNEAR
            TOTALNEAR += 1
            seq_ref = REF_G[info.pos: info.pos + READ_LEN]
            seq_hapA = HAPA_G[g_info.pos: g_info.pos + READ_LEN]
            seq_hapB = HAPB_G[g_info.pos: g_info.pos + READ_LEN]
            leven_score_g = []
            for i in offsets:
                seq_ref_g = REF_G[g_info.pos - i: g_info.pos - i + READ_LEN]
                leven_score_g.append(levenshtein(seq_ref_g, seq_hapA))
                leven_score_g.append(levenshtein(seq_ref_g, seq_hapB))
            called_d = min(levenshtein(seq_ref, seq_hapA), levenshtein(seq_ref, seq_hapB))
            golden_d = min(leven_score_g)
            global HIGHC, CALL_D_ALT, SIM_D_ALT, CALL_D_ORIG, SIM_D_ORIG
            if called_d >= golden_d:
                CALL_D_ORIG.append(called_d)
                SIM_D_ORIG.append(golden_d)
            else:
                HIGHC += 1
                CALL_D_ALT.append(called_d)
                SIM_D_ALT.append(golden_d)
                if called_d > 5 or golden_d > 5 and __debug__:
                    print ('called distance', called_d)
                    print ('golden distance', golden_d)
                    print ('CALLED (%10d) = %s' % (info.pos, REF_G[info.pos : info.pos + 80]))
                    print ('ORIG1  (%10d) = %s' % (g_info.pos - offsets[0], REF_G[g_info.pos - offsets[0] : g_info.pos - offsets[0] + 80]))
                    if offsets[0] != offsets[1]:
                        print ('ORIG2  (%10d) = %s' % (g_info.pos - offsets[1], REF_G[g_info.pos - offsets[1] : g_info.pos - offsets[1] + 80]))
                    print ('PERSON (#%9d) = %s' % (g_info.pos, HAPA_G[g_info.pos : g_info.pos + 80]))
                    print_and_stop(name, offsets, diff, info, g_info)
            return
        if __debug__:
            print_and_stop(name, offsets, diff, info, g_info)

def diploid_compare(
    info, 
    g_info,
    name, 
    threshold, 
    dip_flag, 
    main_offset_index = {}, 
    alt_offset_index = {}
):
    '''
    Uses variable 'dip_flag' to handle different cases of a diploid alignment and check if the alignment is correct.
    '''
    #: don't check the other strand
    if dip_flag in ['same_strand']:
        return compare_sam_info(info, g_info, threshold)
    elif dip_flag in ['same_strand_ref']:
        #: neglect chrom name difference
        # info.chrm = g_info.chrm
        i_low = int(g_info.pos / STEP)
        i_high = math.ceil(g_info.pos / STEP)
        offsets = []
        #: check hapA
        if name.find(MAIN_HAP) > 0:
            if i_low >= len(main_offset_index):
                offsets.append(main_offset_index[len(main_offset_index) - 1])
            else:
                offsets.append(main_offset_index[i_low])
            if i_high >= len(main_offset_index):
                offsets.append(main_offset_index[len(main_offset_index) - 1])
            else:
                offsets.append(main_offset_index[i_high])
        #: check hapB
        if name.find(ALT_HAP) > 0:
            if i_low >= len(alt_offset_index):
                offsets.append(alt_offset_index[len(alt_offset_index) - 1])
            else:
                offsets.append(alt_offset_index[i_low])
            if i_high >= len(alt_offset_index):
                offsets.append(alt_offset_index[len(alt_offset_index) - 1])
            else:
                offsets.append(alt_offset_index[i_high])
        
        comp = compare_sam_info(info, g_info, threshold, offsets, ignore_chrm=True)
        if comp == False and __debug__:
            print_near_aln(name, offsets, info, g_info, 1000)
        return comp
    #: check the other haplotype
    elif dip_flag in ['diff_id', 'diff_var']:
        if info.chrm == MAIN_CHRM:
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
                print_near_aln(name, offsets, info, g_info, 1000)
            return comp
        elif info.chrm == ALT_CHRM:
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
                print_near_aln(name, offsets, info, g_info, 1000)
            return comp
        else:
            print ('Error: invalid chrm', info.chrm)
            exit()
    else:
        print ('Error: undistinguished dip_flag: %s' % dip_flag)
        return False

def count_overlapping_vars(name, info, g_info, main_index, alt_index, MAIN_CHRM, ALT_CHRM, READ_LEN):
    '''
    For an alignment, count the number of overlapping variants.
    Counting is based on raw read (look up golden dictionary)
    '''
    num_var = 0   
    for i in range(g_info.pos, g_info.pos + READ_LEN):
        if g_info.chrm == MAIN_CHRM:
            if main_index.get(i) != None:
                num_var += 1
        elif g_info.chrm == ALT_CHRM:
            if alt_index.get(i) != None:
                num_var += 1
        else:
            print ('Error: unexpected chrm', info.chrm)
            info.print()
            exit()
    return num_var

def analyze_diploid_indels(sam_fn, golden_fn, threshold, var_fn, personalized, write_wrt_correctness):
    '''
    Handles I/O and different opperating modes of this script.
    '''
    # Global variables
    global STEP, MAIN_CHRM, ALT_CHRM, MAIN_HAP, ALT_HAP, MAIN_STRAND, ALT_STRAND, READ_LEN
    MAIN_HAP = 'hapA'
    ALT_HAP = 'hapB'
    MAIN_STRAND = 'A'
    ALT_STRAND = 'B'
    STEP = args.step_size
    READ_LEN = args.read_len
    MAIN_CHRM = '9A'
    ALT_CHRM = '9B'

    var_list = read_var(var_fn, remove_conflict=True, remove_coexist=False)
    main_index, alt_index = build_index(var_list, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    #: diploid personalized ref
    if personalized == 2:
        var_list = read_var(var_fn, remove_conflict=True, remove_coexist=True)
        # main_index, alt_index, 
        main_offset_index, alt_offset_index = build_offset_index(var_list, per=2, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    #: standard ref seq
    elif personalized == 0:
        var_list = read_var(var_fn, remove_conflict=True, remove_coexist=False)
        # main_index, alt_index, 
        main_offset_index, alt_offset_index = build_offset_index_ref(var_list, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
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

    for line in sam_f:
        name, info = parse_line(line, erg=True)
        #: headers
        if name == 'header':
            continue
        summary.add_one()
        if info.is_unaligned():
            comp = False
            summary.add_by_categories(flag='unaligned', comp=comp)
        else:
            #: count the number of overlapping variants for all aligned reads
            num_var = count_overlapping_vars(
                name=name,
                info=info,
                g_info=golden_dic[name],
                main_index=main_index,
                alt_index=alt_index,
                MAIN_CHRM=MAIN_CHRM,
                ALT_CHRM=ALT_CHRM,
                READ_LEN=READ_LEN
            )
        
        #: alignment against personalized genomes
        if (not info.is_unaligned()) and (personalized == 2):
            name_chrm_mismatch = (name.find(MAIN_HAP) > 0 and info.chrm != MAIN_CHRM) or (name.find(ALT_HAP) > 0 and info.chrm != ALT_CHRM)
            #: aligned to incorrect haplotype
            if name_chrm_mismatch:
                if num_var == 0:
                    comp = diploid_compare(info, golden_dic[name], name, threshold, 'diff_id', main_offset_index, alt_offset_index)
                    summary.add_by_categories(flag='diff_id', comp=comp)
                else:
                    comp = diploid_compare(info, golden_dic[name], name, threshold, 'diff_var', main_offset_index, alt_offset_index)
                    summary.add_by_categories(flag='diff_var', comp=comp)
            else:
                comp = diploid_compare(info, golden_dic[name], name, threshold, 'same_strand')
                if num_var == 0:
                    summary.add_by_categories(flag='same_id', comp=comp)
                else:
                    summary.add_by_categories(flag='same_var', comp=comp)
        #: alignment against standard ref (and ERG)
        elif (not info.is_unaligned()) and (personalized == 0):
            comp = diploid_compare(info, golden_dic[name], name, threshold, 'same_strand_ref', main_offset_index, alt_offset_index)
            if num_var == 0:
                summary.add_by_categories(flag='same_id', comp=comp)
            else:
                summary.add_by_categories(flag='same_var', comp=comp)

        if write_wrt_correctness:
            if comp:
                correct_f.write(line)
            else:
                incorrect_f.write(line)

        if __debug__ and comp == False:
            print_and_stop(name, [], None, info, golden_dic[name])

    summary.show_summary(has_answer=True)
    sam_f.close()

def write_wrt_mapq(sam_fn, mapq_threshold):
    sam_f = open(sam_fn, 'r')
    sam_prefix = sam_fn[: sam_fn.find('.')]
    highmapq_fn = sam_prefix + '-mapqgeq' + str(mapq_threshold) + '.sam'
    lowmapq_fn = sam_prefix + '-mapql' + str(mapq_threshold) + '.sam'
    sys.stderr.write('Write sam files %s and %s wrt to mapq...\n' % (highmapq_fn, lowmapq_fn))
    
    highmapq_f = open(highmapq_fn, 'w')
    lowmapq_f = open(lowmapq_fn, 'w')
    for line in sam_f:
        name, info = parse_line(line)
        #: headers
        if name == 'header':
            highmapq_f.write(line)
            lowmapq_f.write(line)
            continue
        if info.mapq >= mapq_threshold:
            highmapq_f.write(line)
        else:
            lowmapq_f.write(line)
        continue
    return

if __name__ == '__main__':
    args = parse_args()
    sam_fn = args.sam
    golden_fn = args.golden
    threshold = args.threshold
    var_fn = args.var
    personalized = args.personalized
    write_wrt_correctness = args.write_wrt_correctness
    write_wrt_mapq = args.write_wrt_mapq

    if write_wrt_mapq:
        write_wrt_mapq(sam_fn, int(write_wrt_mapq))
        exit()

    global COMPARE_SEQ, HIGHC, TOTALNEAR, REF_G, HAPA_G, HAPB_G, CALL_D_ALT, SIM_D_ALT, CALL_D_ORIG, SIM_D_ORIG
    COMPARE_SEQ = False
    if COMPARE_SEQ:
        HIGHC = 0
        TOTALNEAR = 0
        CALL_D_ALT = []
        SIM_D_ALT = []
        CALL_D_ORIG = []
        SIM_D_ORIG = []
        #TODO
        fn_ref =  '/scratch/groups/blangme2/naechyun/relaxing/chr9/chr9_singleline.fa'
        fn_hapA = '/scratch/groups/blangme2/naechyun/relaxing/na12878/indels/hapA_single.fa'
        fn_hapB = '/scratch/groups/blangme2/naechyun/relaxing/na12878/indels/hapB_single.fa'
        REF_G = read_genome(fn_ref)
        HAPA_G = read_genome(fn_hapA)
        HAPB_G = read_genome(fn_hapB)

    analyze_diploid_indels(sam_fn, golden_fn, threshold, var_fn, personalized, write_wrt_correctness)

    if COMPARE_SEQ:
        print ('Number of alns have higher score than golden =', HIGHC)
        print ('Total number of near alignments =', TOTALNEAR)
        print ('Avg Lev. dist of called ALT alignments =', sum(CALL_D_ALT)/len(CALL_D_ALT))
        print (CALL_D_ALT)
        print ('Avg Lev. dist of simulated ALT alignments =', sum(SIM_D_ALT)/len(SIM_D_ALT))
        print (SIM_D_ALT)
        print ('Avg Lev. dist of called ORIG alignments =', sum(CALL_D_ORIG)/(TOTALNEAR-HIGHC))
        print ('Avg Lev. dist of simulated ORIG alignments =', sum(SIM_D_ORIG)/(TOTALNEAR-HIGHC))
