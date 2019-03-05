'''
Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse, math, sys, os
import pandas as pd
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

def build_offset_index(var_list, step, MAIN_STRAND, ALT_STRAND):
    '''
    MAIN/ALT-offset indexes are dictionaries with
        key: pos on MAIN/ALT
        value: pos on ALT/MAIN
    
    MAIN/ALT indexes are dictionaries storing
    variants based on MAIN/ALT coordinates
    '''
    #: dict storing the diff from main to alt
    #: main_pos + main_offset_index[i] = alt_pos
    main_offset_index = [0]
    #: dict storing the diff from alt to main
    #: alt_pos + alt_offset_index[i] = main_pos
    alt_offset_index = [0]
    SHOW_DUP_WARN = False
    tmp_v = 0
    for v in var_list:
        if v.strand == MAIN_STRAND:
            main_pos = v.alt_pos
            alt_pos = v.ref_pos + v.cor_offset
            main_offset = -v.offset + v.cor_offset
            alt_offset = v.offset - v.cor_offset
            if v.ref_pos == tmp_v:
                if SHOW_DUP_WARN:
                    print ('Warning: duplicated variant', v.line)
            tmp_v = v.ref_pos
        elif v.strand == ALT_STRAND:
            main_pos = v.ref_pos + v.cor_offset
            alt_pos = v.alt_pos
            main_offset = v.offset - v.cor_offset
            alt_offset = -v.offset + v.cor_offset
        else:
            print ('Error: unspecified strand', v.strand)
            exit()
        
        i_main = math.ceil(main_pos / step)
        while i_main >= len(main_offset_index):
            main_offset_index.append(main_offset_index[len(main_offset_index) - 1])
        main_offset_index[i_main] = main_offset
        i_alt = math.ceil(alt_pos / step)
        while i_alt >= len(alt_offset_index):
            alt_offset_index.append(alt_offset_index[len(alt_offset_index) - 1])
        alt_offset_index[i_alt] = alt_offset
    
    return main_offset_index, alt_offset_index

def build_offset_index_ref(var_list, step, MAIN_STRAND, ALT_STRAND):
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
        idx = math.ceil(v.alt_pos / step)
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

# TODO
def print_aln_within_distance(name, offsets, info, g_info, threshold, read_len):
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
            seq_ref = REF_G[info.pos: info.pos + read_len]
            seq_hapA = HAPA_G[g_info.pos: g_info.pos + read_len]
            seq_hapB = HAPB_G[g_info.pos: g_info.pos + read_len]
            leven_score_g = []
            for i in offsets:
                seq_ref_g = REF_G[g_info.pos - i: g_info.pos - i + read_len]
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
    step,
    main_offset_index = {}, 
    alt_offset_index = {}
):
    '''
    Uses variable 'dip_flag' to handle different cases of a diploid alignment and check if the alignment is correct.
    '''
    #: don't check the other strand
    if dip_flag in ['same_id', 'same_var']:
        return compare_sam_info(info, g_info, threshold)
    elif dip_flag in ['same_strand_ref']:
        #: neglect chrom name difference
        # info.chrm = g_info.chrm
        i_low = int(g_info.pos / step)
        i_high = math.ceil(g_info.pos / step)
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
        elif name.find(ALT_HAP) > 0:
            if i_low >= len(alt_offset_index):
                offsets.append(alt_offset_index[len(alt_offset_index) - 1])
            else:
                offsets.append(alt_offset_index[i_low])
            if i_high >= len(alt_offset_index):
                offsets.append(alt_offset_index[len(alt_offset_index) - 1])
            else:
                offsets.append(alt_offset_index[i_high])
    #: check the other haplotype
    elif dip_flag in ['diff_id', 'diff_var']:
        i_low = int(info.pos / step)
        i_high = math.ceil(info.pos / step)
        if info.chrm == MAIN_CHRM:
            if i_low >= len(main_offset_index):
                offset_low = main_offset_index[len(main_offset_index) - 1]
            else:
                offset_low = main_offset_index[i_low]
            if i_high >= len(main_offset_index):
                offset_high = main_offset_index[len(main_offset_index) - 1]
            else:
                offset_high = main_offset_index[i_high]
        elif info.chrm == ALT_CHRM:
            if i_low >= len(alt_offset_index):
                offset_low = alt_offset_index[len(alt_offset_index) - 1]
            else:
                offset_low = alt_offset_index[i_low]
            if i_high >= len(alt_offset_index):
                offset_high = alt_offset_index[len(alt_offset_index) - 1]
            else:
                offset_high = alt_offset_index[i_high]
        else:
            print ('Error: invalid chrm', info.chrm)
            exit()
        offsets = [offset_low, offset_high]
    else:
        print ('Error: undistinguished dip_flag: %s' % dip_flag)
        return False

    comp = compare_sam_info(info, g_info, threshold, offsets, ignore_chrm=True)
    if comp == False and __debug__:
        print_aln_within_distance(name, offsets, info, g_info, 1000, read_len)
    return comp

def count_overlapping_vars(name, info, g_info, main_index, alt_index, MAIN_CHRM, ALT_CHRM, read_len):
    '''
    For an alignment, count the number of overlapping variants.
    Counting is based on raw read (look up golden dictionary)
    '''
    num_var = 0   
    for i in range(g_info.pos, g_info.pos + read_len):
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

def analyze_diploid_indels(
    sam_fn,
    golden_fn,
    threshold,
    var_fn,
    personalized,
    step,
    read_len,
    write_wrt_correctness
):
    '''
    Handles I/O and different opperating modes of this script.
    '''
    # Global variables
    global MAIN_CHRM, ALT_CHRM, MAIN_HAP, ALT_HAP, MAIN_STRAND, ALT_STRAND
    MAIN_STRAND = 'A'
    ALT_STRAND = 'B'
    MAIN_HAP = 'hap' + MAIN_STRAND
    ALT_HAP = 'hap' + ALT_STRAND
    MAIN_CHRM = '9A'
    ALT_CHRM = '9B'

    var_list = read_var(var_fn, remove_conflict=True, remove_coexist=False)
    main_index, alt_index = build_index(var_list, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    #: diploid personalized ref
    if personalized == 2:
        var_list = read_var(var_fn, remove_conflict=True, remove_coexist=True)
        main_offset_index, alt_offset_index = build_offset_index(var_list, step, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    #: standard ref seq
    elif personalized == 0:
        # var_list = read_var(var_fn, remove_conflict=True, remove_coexist=False)
        main_offset_index, alt_offset_index = build_offset_index_ref(var_list, step, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    else:
        print ('Error: unsupported personalized parameter', personalized)
        exit()
    
    sam_f = open(sam_fn, 'r')
    sam_prefix = sam_fn[: sam_fn.find('.')]
    golden_dic = load_golden_dic(golden_fn, 1)
    summary = Summary(has_answer=True)
    results = []
    
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

        #: count the number of overlapping variants for all aligned reads
        num_var = count_overlapping_vars(
            name=name,
            info=info,
            g_info=golden_dic[name],
            main_index=main_index,
            alt_index=alt_index,
            MAIN_CHRM=MAIN_CHRM,
            ALT_CHRM=ALT_CHRM,
            read_len=read_len
        )

        if info.is_unaligned():
            dist = -3
            comp = False
            flag = 'unaligned'
            summary.add_by_categories(flag=flag, comp=comp)
        #: alignment against personalized genomes
        elif personalized == 2:
            name_chrm_mismatch = (name.find(MAIN_HAP) > 0 and info.chrm != MAIN_CHRM) or (name.find(ALT_HAP) > 0 and info.chrm != ALT_CHRM)
            #: aligned to incorrect haplotype
            if name_chrm_mismatch:
                if num_var == 0:
                    flag = 'diff_id'
                else:
                    flag = 'diff_var'
            else:
                if num_var == 0:
                    flag = 'same_id'
                else:
                    flag = 'same_var'
            dist = diploid_compare(info, golden_dic[name], name, threshold, flag, step, main_offset_index, alt_offset_index)
            if dist < 0 or dist > threshold:
                comp = False
            else:
                comp = True
            summary.add_by_categories(flag=flag, comp=comp)
        #: alignment against standard ref (and ERG)
        elif personalized == 0:
            dist = diploid_compare(info, golden_dic[name], name, threshold, 'same_strand_ref', step, main_offset_index, alt_offset_index)
            if dist < 0 or dist > threshold:
                comp = False
            else:
                comp = True
            if num_var == 0:
                flag = 'same_id'
            else:
                flag = 'same_var'
            summary.add_by_categories(flag=flag, comp=comp)
        
        results.append([name, dist, info.mapq, num_var, flag])

        if write_wrt_correctness:
            if comp:
                correct_f.write(line)
            else:
                incorrect_f.write(line)

        if __debug__ and comp == False:
            print_and_stop(name, [], None, info, golden_dic[name])

    summary.show_summary(has_answer=True)
    sam_f.close()

    results_df = pd.DataFrame(results, columns=['name', 'dist', 'mapq', 'numvar', 'category'])
    results_df.to_pickle(sam_fn + '-stats.pkl')

    return results_df

def print_df_stats(df, threshold, var_opt):
    print ()
    print ('--- Stats ---')
    unaligned = (df['dist'] == -3)
    correct = (df['dist'] >= 0) & (df['dist'] <= threshold)
    aligned_incorrect = (df['dist'] == -1) | (df['dist'] == -2)

    #: categories
    cat_same_id = (df['category'] == 'same_id')
    cat_same_var = (df['category'] == 'same_var')
    cat_diff_id = (df['category'] == 'diff_id')
    cat_diff_var = (df['category'] == 'diff_var')

    sensitivity_all = df[correct].shape[0] / df.shape[0]
    print ('sensitivity_all      = {0:.4%} ({1} / {2})'.format(sensitivity_all, df[correct].shape[0], df.shape[0]))
    precision_all = df[correct].shape[0] / (df.shape[0] - df[unaligned].shape[0])
    print ('precision_all        = {0:.4%} ({1} / {2})'.format(precision_all, df[correct].shape[0], df.shape[0] - df[unaligned].shape[0]))
    fdr_all = 1 - precision_all
    print ('fdr_all = {0:.4%}'.format(fdr_all))
    unaligned_rate = df[unaligned].shape[0] / df.shape[0]
    print ('unaligned            = {0:.4%} ({1} / {2})'.format(unaligned_rate, df[unaligned].shape[0], df.shape[0]))
    try:
        sensitivity_same_id = df[correct & cat_same_id].shape[0] / df[cat_same_id].shape[0]
        print ('sensitivity_same_id  = {0:.4%} ({1} / {2})'.format(sensitivity_same_id, df[correct & cat_same_id].shape[0], df[cat_same_id].shape[0]))
    except:
        print ('Warning: no element in "same_id"')
    try:
        sensitivity_same_var = df[correct & cat_same_var].shape[0] / df[cat_same_var].shape[0]
        print ('sensitivity_same_var = {0:.4%} ({1} / {2})'.format(sensitivity_same_var, df[correct & cat_same_var].shape[0], df[cat_same_var].shape[0]))
    except:
        print ('Warning: no element in "same_var"')
    try:
        sensitivity_diff_id = df[correct & cat_diff_id].shape[0] / df[cat_diff_id].shape[0]
        print ('sensitivity_diff_id  = {0:.4%} ({1} / {2})'.format(sensitivity_diff_id, df[correct & cat_diff_id].shape[0], df[cat_diff_id].shape[0]))
    except:
        print ('Warning: no element in "diff_id"')
    try:
        sensitivity_diff_var = df[correct & cat_diff_var].shape[0] / df[cat_diff_var].shape[0]
        print ('sensitivity_diff_var = {0:.4%} ({1} / {2})'.format(sensitivity_diff_var, df[correct & cat_diff_var].shape[0], df[cat_diff_var].shape[0]))
    except:
        print ('Warning: no element in "diff_var"')

    #: number of overlapping variants
    print ()
    var_all = (df['numvar'] >= 0)
    var0 = (df['numvar'] == 0)
    var1 = (df['numvar'] == 1)
    var2 = (df['numvar'] == 2)
    var3plus = (df['numvar'] >= 3)
    
    if var_opt == 'all': 
        v_filter = var_all
    elif var_opt == '0':
        v_filter = var0
    elif var_opt == '1':
        v_filter = var1
    elif var_opt == '2':
        v_filter = var2
    elif var_opt == '3+':
        v_filter = var3plus

    try:
        sensitivity_var0 = df[correct & var0].shape[0] / df[var0].shape[0]
        print ('sensitivity_var0     = {0:.4%} ({1} / {2})'.format(sensitivity_var0, df[correct & var0].shape[0], df[var0].shape[0]))
    except:
        print ('Warning: no read has 0 variants')
    try:
        sensitivity_var1 = df[correct & var1].shape[0] / df[var1].shape[0]
        print ('sensitivity_var1     = {0:.4%} ({1} / {2})'.format(sensitivity_var1, df[correct & var1].shape[0], df[var1].shape[0]))
    except:
        print ('Warning: no read has 1 variants')
    try:
        sensitivity_var2 = df[correct & var2].shape[0] / df[var2].shape[0]
        print ('sensitivity_var2     = {0:.4%} ({1} / {2})'.format(sensitivity_var2, df[correct & var2].shape[0], df[var2].shape[0]))
    except:
        print ('Warning: no read has 2 variants')
    try:
        sensitivity_var3plus = df[correct & var3plus].shape[0] / df[var3plus].shape[0]
        print ('sensitivity_var3plus = {0:.4%} ({1} / {2})'.format(sensitivity_var3plus, df[correct & var3plus].shape[0], df[var3plus].shape[0]))
    except:
        print ('Warning: no read has 3+ variants')
    
    #: mapq
    print ()
    mapq5plus = (df['mapq'] >= 5)
    mapq10plus = (df['mapq'] >= 10)
    mapq20plus = (df['mapq'] >= 20)
    mapq30plus = (df['mapq'] >= 30)
    mapq40plus = (df['mapq'] >= 40)

    sensitivity_mapqall = df[correct & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapqall    = {0:.4%} ({1} / {2})'.format(sensitivity_mapqall, df[correct & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq5plus = df[correct & mapq5plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq5plus  = {0:.4%} ({1} / {2})'.format(sensitivity_mapq5plus, df[correct & mapq5plus & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq10plus = df[correct & mapq10plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq10plus = {0:.4%} ({1} / {2})'.format(sensitivity_mapq10plus, df[correct & mapq10plus & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq20plus = df[correct & mapq20plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq20plus = {0:.4%} ({1} / {2})'.format(sensitivity_mapq20plus, df[correct & mapq20plus & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq30plus = df[correct & mapq30plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq30plus = {0:.4%} ({1} / {2})'.format(sensitivity_mapq30plus, df[correct & mapq30plus & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq40plus = df[correct & mapq40plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq40plus = {0:.4%} ({1} / {2})'.format(sensitivity_mapq40plus, df[correct & mapq40plus & v_filter].shape[0], df[v_filter].shape[0]))
    
    print ()
    try:
        precision_mapqall = df[correct & v_filter].shape[0] / df[v_filter].shape[0]
        print ('precision_mapqall    = {0:.4%} ({1} / {2})'.format(precision_mapqall, df[correct & v_filter].shape[0], df[v_filter].shape[0]))
        fdr_mapqall = 1 - precision_mapqall
        print ('fdr_mapqall = {0:.4%}'.format(fdr_mapqall))
    except:
        print ('Warning: no read is 0+ mapq')
    try:
        precision_mapq5plus = df[correct & mapq5plus & v_filter].shape[0] / df[mapq5plus & v_filter].shape[0]
        print ('precision_mapq5plus  = {0:.4%} ({1} / {2})'.format(precision_mapq5plus, df[correct & mapq5plus & v_filter].shape[0], df[mapq5plus & v_filter].shape[0]))
        fdr_mapq5plus = 1 - precision_mapq5plus
        print ('fdr_mapq5plus = {0:.4%}'.format(fdr_mapq5plus))
    except:
        print ('Warning: no read is 5+ mapq')
    try:
        precision_mapq10plus = df[correct & mapq10plus & v_filter].shape[0] / df[mapq10plus & v_filter].shape[0]
        print ('precision_mapq10plus = {0:.4%} ({1} / {2})'.format(precision_mapq10plus, df[correct & mapq10plus & v_filter].shape[0], df[mapq10plus & v_filter].shape[0]))
        fdr_mapq10plus = 1 - precision_mapq10plus
        print ('fdr_mapq10plus = {0:.4%}'.format(fdr_mapq10plus))
    except:
        print ('Warning: no read is 10+ mapq')
    try:
        precision_mapq20plus = df[correct & mapq20plus & v_filter].shape[0] / df[mapq20plus & v_filter].shape[0]
        print ('precision_mapq20plus = {0:.4%} ({1} / {2})'.format(precision_mapq20plus, df[correct & mapq20plus & v_filter].shape[0], df[mapq20plus & v_filter].shape[0]))
        fdr_mapq20plus = 1 - precision_mapq20plus
        print ('fdr_mapq20plus = {0:.4%}'.format(fdr_mapq20plus))
    except:
        print ('Warning: no read is 20+ mapq')
    try:
        precision_mapq30plus = df[correct & mapq30plus & v_filter].shape[0] / df[mapq30plus & v_filter].shape[0]
        print ('precision_mapq30plus = {0:.4%} ({1} / {2})'.format(precision_mapq30plus, df[correct & mapq30plus & v_filter].shape[0], df[mapq30plus & v_filter].shape[0]))
        fdr_mapq30plus = 1 - precision_mapq30plus
        print ('fdr_mapq30plus = {0:.4%}'.format(fdr_mapq30plus))
    except:
        print ('Warning: no read is 30+ mapq')
    try:
        precision_mapq40plus = df[correct & mapq40plus & v_filter].shape[0] / df[mapq40plus & v_filter].shape[0]
        print ('precision_mapq40plus = {0:.4%} ({1} / {2})'.format(precision_mapq40plus, df[correct & mapq40plus & v_filter].shape[0], df[mapq40plus & v_filter].shape[0]))
        fdr_mapq40plus = 1 - precision_mapq40plus
        print ('fdr_mapq40plus = {0:.4%}'.format(fdr_mapq40plus))
    except:
        print ('Warning: no read is 40+ mapq')
    
    x = [fdr_mapqall, fdr_mapq5plus, fdr_mapq10plus, fdr_mapq20plus, fdr_mapq30plus, fdr_mapq40plus]
    y = [sensitivity_mapqall, sensitivity_mapq5plus, sensitivity_mapq10plus, sensitivity_mapq20plus, sensitivity_mapq30plus, sensitivity_mapq40plus]

    return x, y

def plot_ROC(list_roc, title):
    q = [0, 5, 10, 20, 30, 40]

    fig, ax = plt.subplots()
    plt.title(title)
    plt.xscale('log')
    plt.xlabel('FDR (1-precision)')
    plt.ylabel('TPR (sensitivity)')
    
    marker_dict={0:'.', 1:'x', 2:'+'}
    
    for iroc, roc in enumerate(list_roc):
        x = roc[0]
        y = roc[1]
        # ax.scatter(x, y, label=roc[2])
        c = 'C' + str(iroc)
        ax.plot(x, y, color=c, label=roc[2], marker=marker_dict[iroc%3])
        for i, txt in enumerate(q):
            ax.annotate(txt, (x[i], y[i]), color=c)
            # ax.annotate(txt, (x[i], y[i]+0.005*(iroc==0)), color=c)
    ax.legend()
    # plt.ylim(0.8, 0.95)
    # plt.show()
    plt.savefig('roc_' + title.lower() + '.pdf')

if __name__ == '__main__':
    args = parse_args()
    sam_fn = args.sam
    golden_fn = args.golden
    threshold = args.threshold
    var_fn = args.var
    personalized = args.personalized
    write_wrt_correctness = args.write_wrt_correctness
    mapq_threshold = args.write_wrt_mapq
    step = args.step_size
    read_len = args.read_len

    if mapq_threshold:
        write_wrt_mapq(sam_fn, int(mapq_threshold))
        exit()

    if os.path.isfile(sam_fn + '-stats.pkl'):
        print ('Read stats from {0}-stats.pkl'.format(sam_fn))
        df = pd.read_pickle(sam_fn + '-stats.pkl')
        print_df_stats(df, threshold, 'all')
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

    results_df = analyze_diploid_indels(
        sam_fn=sam_fn,
        golden_fn=golden_fn,
        threshold=threshold,
        var_fn=var_fn,
        personalized=personalized,
        step=step,
        read_len=read_len,
        write_wrt_correctness=write_wrt_correctness
    )

    print_df_stats(results_df, threshold, 'all')

    if COMPARE_SEQ:
        print ('Number of alns have higher score than golden =', HIGHC)
        print ('Total number of near alignments =', TOTALNEAR)
        print ('Avg Lev. dist of called ALT alignments =', sum(CALL_D_ALT)/len(CALL_D_ALT))
        print (CALL_D_ALT)
        print ('Avg Lev. dist of simulated ALT alignments =', sum(SIM_D_ALT)/len(SIM_D_ALT))
        print (SIM_D_ALT)
        print ('Avg Lev. dist of called ORIG alignments =', sum(CALL_D_ORIG)/(TOTALNEAR-HIGHC))
        print ('Avg Lev. dist of simulated ORIG alignments =', sum(SIM_D_ORIG)/(TOTALNEAR-HIGHC))
