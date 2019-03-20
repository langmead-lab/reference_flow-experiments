'''
Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse, math, sys, os
import pandas as pd
import matplotlib.pyplot as plt
from analyze_sam import SamInfo, parse_line, load_golden_dic, Summary
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
        '-vr', '--var_reads',
        help='the file specifying variants for the synthetic reads'
    )
    parser.add_argument(
        '-vs', '--var_sample',
        help='the file specifying variants for the current target reference'
    )
    parser.add_argument(
        '-p', '--personalized', type=int,
        default=0,
        help='(int) specify whether the ref seq(s) are standard (0) or personalized-diploid (2) sample [0]'
    )
    parser.add_argument(
        '-c', '--chrm',
        help='(str) chromosome [None]'
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
    parser.add_argument(
        '--debug_ref',
        help='reference fasta for debug purpose [None]'
    )
    parser.add_argument(
        '--debug_hapA',
        help='hapA fasta for debug purpose [None]'
    )
    parser.add_argument(
        '--debug_hapB',
        help='hapB fasta for debug purpose [None]'
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
    CURRENTLY UNUSED
    
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

def print_and_stop(name, reads_offsets, sample_offsets, dist, info, g_info, flag):
    '''
    This is for debugging.
    Prints alignment info and stops the script.
    '''
    print ('name', name)
    print ('offsets on reads', reads_offsets)
    print ('offsets on sample', sample_offsets)
    print ('dist', dist)
    print ('flag', flag)
    print ('info')
    info.print(flag=False, mapq=False, score=False)
    print ()
    print ('golden')
    g_info.print(flag=False, mapq=False, score=False)
    input()

# TODO
def print_aln_within_distance(name, reads_offsets, sample_offsets, info, g_info, threshold, read_len):
    '''
    Compares alignment with the golden profile if they are near.
    If COMPARE_SEQ is specified, retrieves sequences from ref and haps and calculate the distance.
    '''
    if COMPARE_SEQ == False:
        return
    tmp = []
    for i in reads_offsets:
        tmp.append(abs(info.pos + i - g_info.pos))
    diff = min(tmp)
    if (diff < threshold) or (threshold < 0):
        global TOTALNEAR
        TOTALNEAR += 1
        seq_ref = REF_G[info.pos: info.pos + read_len]
        seq_hapA = HAPA_G[g_info.pos: g_info.pos + read_len]
        seq_hapB = HAPB_G[g_info.pos: g_info.pos + read_len]
        leven_score_g = []
        for i in reads_offsets:
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
                print ('ORIG1  (%10d) = %s' % (g_info.pos - reads_offsets[0], REF_G[g_info.pos - offsets[0] : g_info.pos - reads_offsets[0] + 80]))
                if reads_offsets[0] != reads_offsets[1]:
                    print ('ORIG2  (%10d) = %s' % (g_info.pos - reads_offsets[1], REF_G[g_info.pos - offsets[1] : g_info.pos - reads_offsets[1] + 80]))
                print ('PERSON (#%9d) = %s' % (g_info.pos, HAPA_G[g_info.pos : g_info.pos + 80]))
        return

def compare_sam_info(
        info,
        ginfo,
        threshold,
        sample_offsets=[0],
        reads_offsets=[0],
        ignore_chrm=False
    ):
    '''
    Inputs:
        info:
            info from alignment
        ginfo:
            info from simulation profile (golden)
        offset:
            positiontal offset
        ignore_chrm:
            set True to ignore alignment against different chromosomes
    
    Output:
        an INT representing alignment correctness
        if < 0:
            -1: unmatched chromosome
            -2: unmatched direction
        if >= 0:
            the difference in alignment position
            0 is a perfect match
    '''
    if (ignore_chrm is False) and (info.chrm != ginfo.chrm):
        #: diff chromosome
        if __debug__:
            print ("False: chr, mapq =", info.mapq)
        return -1
    if (info.is_rc() ^ ginfo.is_rc()) is True:
        #: diff direction
        if __debug__: 
            print ("False: direction (%s, %s)" % (info.is_rc(), ginfo.is_rc()), "mapq =", info.mapq)
        return -2
    dist = []
    for soff in sample_offsets:
        for roff in reads_offsets:
            dist.append(abs(info.pos - soff - ginfo.pos + roff))
    return min(dist)
    #return min([abs(info.pos + off - ginfo.pos) for off in sample_offsets])

def diploid_compare(
    info, 
    g_info,
    name, 
    threshold, 
    dip_flag, 
    step,
    reads_main_offset_index = {}, 
    reads_alt_offset_index = {},
    sample_main_offset_index = {},
    sample_alt_offset_index = {}
):
    '''
    Uses variable 'dip_flag' to handle different cases of a diploid alignment 
    and check if the alignment is correct.
    '''
    sample_offsets = [0]
    if dip_flag in ['same_strand_ref']:
        if sample_main_offset_index != {}:
            i_low = int(info.pos / step)
            i_high = math.ceil(info.pos / step)
            if i_low >= len(sample_main_offset_index):
                sample_offset_low = sample_main_offset_index[len(sample_main_offset_index) - 1]
            else:
                sample_offset_low = sample_main_offset_index[i_low]
            if i_high >= len(sample_main_offset_index):
                sample_offset_high = sample_main_offset_index[len(sample_main_offset_index) - 1]
            else:
                sample_offset_high = sample_main_offset_index[i_high]
            sample_offsets = [sample_offset_low, sample_offset_high]
        #pass
    elif dip_flag in ['same_id', 'same_var', 'diff_id', 'diff_var']:
        i_low = int(info.pos / step)
        i_high = math.ceil(info.pos / step)
        if info.chrm == MAIN_CHRM:
            if i_low >= len(sample_main_offset_index):
                sample_offset_low = sample_main_offset_index[len(sample_main_offset_index) - 1]
            else:
                sample_offset_low = sample_main_offset_index[i_low]
            if i_high >= len(sample_main_offset_index):
                sample_offset_high = sample_main_offset_index[len(sample_main_offset_index) - 1]
            else:
                sample_offset_high = sample_main_offset_index[i_high]
        elif info.chrm == ALT_CHRM:
            if i_low >= len(sample_alt_offset_index):
                sample_offset_low = sample_alt_offset_index[len(sample_alt_offset_index) - 1]
            else:
                sample_offset_low = sample_alt_offset_index[i_low]
            if i_high >= len(sample_alt_offset_index):
                sample_offset_high = sample_alt_offset_index[len(sample_alt_offset_index) - 1]
            else:
                sample_offset_high = sample_alt_offset_index[i_high]
        else:
            print ('Error: invalid chrm', info.chrm, MAIN_CHRM, ALT_CHRM)
            exit()
        sample_offsets = [sample_offset_low, sample_offset_high]
    else:
        print ('Error: undistinguished dip_flag: %s' % dip_flag)
        return False
    
    i_low = int(g_info.pos / step)
    i_high = math.ceil(g_info.pos / step)
    reads_offsets = []
    #: check hapA
    if name.find(MAIN_HAP) > 0:
        if i_low >= len(reads_main_offset_index):
            reads_offsets.append(reads_main_offset_index[len(reads_main_offset_index) - 1])
        else:
            reads_offsets.append(reads_main_offset_index[i_low])
        if i_high >= len(reads_main_offset_index):
            reads_offsets.append(reads_main_offset_index[len(reads_main_offset_index) - 1])
        else:
            reads_offsets.append(reads_main_offset_index[i_high])
    #: check hapB
    elif name.find(ALT_HAP) > 0:
        if i_low >= len(reads_alt_offset_index):
            reads_offsets.append(reads_alt_offset_index[len(reads_alt_offset_index) - 1])
        else:
            reads_offsets.append(reads_alt_offset_index[i_low])
        if i_high >= len(reads_alt_offset_index):
            reads_offsets.append(reads_alt_offset_index[len(reads_alt_offset_index) - 1])
        else:
            reads_offsets.append(reads_alt_offset_index[i_high])

    dist = compare_sam_info(
        info=info,
        ginfo=g_info,
        threshold=threshold,
        sample_offsets=sample_offsets,
        reads_offsets=reads_offsets,
        ignore_chrm=True
    )

    if (dist < 0 or dist > threshold) and __debug__:
        print_and_stop(name, reads_offsets, sample_offsets, dist, info, g_info, dip_flag)
        #print_aln_within_distance(name, reads_offsets, sample_offsets, info, g_info, 1000, read_len)
    if (dist < 0 or dist > threshold) and COMPARE_SEQ:
        print_aln_within_distance(
            name=name,
            reads_offsets=reads_offsets,
            sample_offsets=sample_offsets,
            info=info,
            g_info=g_info,
            threshold=1000,
            read_len=read_len
        )

    return dist

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

def build_all_indexes(
    var_reads_fn,
    var_sample_fn,
    personalized
):
    '''
    Reads two var files and builds all the indexes we use for computing correctness
    '''
    var_reads_list = read_var(var_reads_fn, remove_conflict=True, remove_coexist=False)
    main_index, alt_index = build_index(var_reads_list, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    #: diploid personalized ref
    if personalized == 2:
        var_reads_list = read_var(var_reads_fn, remove_conflict=True, remove_coexist=False)
        reads_main_offset_index, reads_alt_offset_index = build_offset_index_ref(var_reads_list, step, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
        var_sample_list = read_var(var_sample_fn, remove_conflict=True, remove_coexist=False)
        sample_main_offset_index, sample_alt_offset_index = build_offset_index_ref(var_sample_list, step, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    #: standard ref seq
    elif personalized == 0:
        reads_main_offset_index, reads_alt_offset_index = build_offset_index_ref(var_reads_list, step, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
        #: major allele reference with indels
        if var_sample_fn != None:
            var_sample_list = read_var(var_sample_fn, remove_conflict=True, remove_coexist=False)
            sample_main_offset_index, _ = build_offset_index_ref(var_sample_list, step, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
        else:
            sample_main_offset_index = {}
        sample_alt_offset_index = {}
    else:
        print ('Error: unsupported personalized parameter', personalized)
        exit()
    return main_index, alt_index, reads_main_offset_index, reads_alt_offset_index, sample_main_offset_index, sample_alt_offset_index

def analyze_diploid_indels(
    sam_fn,
    golden_dic,
    threshold,
    all_indexes,
    personalized,
    chrm,
    step,
    read_len,
    write_wrt_correctness
):
    '''
    Handles I/O and different opperating modes of this script.
    '''
    main_index, alt_index, reads_main_offset_index, reads_alt_offset_index, sample_main_offset_index, sample_alt_offset_index = all_indexes
    
    sam_f = open(sam_fn, 'r')
    sam_prefix = sam_fn[: sam_fn.find('.')]
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

        #: counts the number of overlapping variants for all aligned reads
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
            flag = 'unaligned'
        #: alignment against personalized genomes
        elif personalized == 2:
            #: aligned to incorrect haplotype
            name_chrm_mismatch = (name.find(MAIN_HAP) > 0 and info.chrm != MAIN_CHRM) or (name.find(ALT_HAP) > 0 and info.chrm != ALT_CHRM)
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
            dist = diploid_compare(
                info=info, 
                g_info=golden_dic[name],
                name=name, 
                threshold=threshold, 
                dip_flag=flag, 
                step=step,
                reads_main_offset_index = reads_main_offset_index, 
                reads_alt_offset_index = reads_alt_offset_index,
                sample_main_offset_index = sample_main_offset_index,
                sample_alt_offset_index = sample_alt_offset_index
            )
        #: alignment against standard ref (and ERG)
        elif personalized == 0:
            flag = 'same_strand_ref'
            dist = diploid_compare(
                info=info, 
                g_info=golden_dic[name],
                name=name, 
                threshold=threshold, 
                dip_flag=flag, 
                step=step,
                reads_main_offset_index = reads_main_offset_index, 
                reads_alt_offset_index = reads_alt_offset_index,
                sample_main_offset_index = sample_main_offset_index,
                sample_alt_offset_index = sample_alt_offset_index
            )
            if num_var == 0:
                flag = 'same_id'
            else:
                flag = 'same_var'

        #: converts "dist" to binary comparsion decision "comp" and adds to summary
        if dist < 0 or dist > threshold:
            comp = False
        else:
            comp = True
        summary.add_by_categories(flag=flag, comp=comp)
        results.append([name, dist, info.mapq, num_var, flag])

        if write_wrt_correctness:
            if comp:
                correct_f.write(line)
            else:
                incorrect_f.write(line)

        # if __debug__ and comp == False:
        #     print_and_stop(name, [], [], dist, info, golden_dic[name], flag)

    summary.show_summary(has_answer=True)
    sam_f.close()

    results_df = pd.DataFrame(results, columns=['name', 'dist', 'mapq', 'numvar', 'category'])
    results_df.to_pickle(sam_fn + '-stats.pkl')

    return results_df

def plot_hist(filename, df):
    '''
    Plots a histogram with the distance between an alignment and its golden postion as x-axis
    '''
    x = df['dist']
    x = x.replace([-3, -2, 0], [0.001, 0.01, 1])
    upper = 1000000000
    rr = [0.001, 0.01, 0.1, 1, 10,100,1000,10000,100000, 1000000, 10000000, 100000000, 1000000000]
    n, bins, patches = plt.hist(x=x, bins=rr, log=True)
    # n_upper_outliers = (x > upper).sum()
    # patches[-1].set_height(patches[-1].get_height() + n_upper_outliers)
    # patches[-1].set_facecolor('m')
    # patches[-1].set_label('including longer')
    
    # patches[0].set_facecolor('#e00000') #: red
    patches[0].set_facecolor('#ff9900')
    patches[0].set_label('unaligned')
    # patches[1].set_facecolor('#ff9900') #: orange
    patches[1].set_facecolor('#ffdb4d')
    patches[1].set_label('diff direction')
    patches[3].set_facecolor('g')
    patches[3].set_label('correct')
    plt.ylabel('counts')
    plt.ylim(300, 300000)
    plt.grid(True)
    plt.xlabel('distance')
    plt.xscale('log')
    plt.legend()
    plt.title(filename.split('/')[-1])
    plt.savefig(filename+'-hist.pdf', format='pdf')
    print (filename+'-hist.pdf is plotted')
    return 

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
    print ('fdr_all              = {0:.4%}'.format(fdr_all))
    unaligned_rate = df[unaligned].shape[0] / df.shape[0]
    print ('unaligned            = {0:.4%} ({1} / {2})'.format(unaligned_rate, df[unaligned].shape[0], df.shape[0]))
    
    print ()
    try:
        sensitivity_same_id = df[correct & cat_same_id].shape[0] / df[cat_same_id].shape[0]
        print ('sensitivity_same_id  = {0:.4%} ({1} / {2})'.format(sensitivity_same_id, df[correct & cat_same_id].shape[0], df[cat_same_id].shape[0]))
        fnr_same_id = 1 - sensitivity_same_id
        print ('fnr_same_id          = {0:.4%} ({1} / {2})'.format(fnr_same_id, df[cat_same_id].shape[0] - df[correct & cat_same_id].shape[0], df[cat_same_id].shape[0]))
    except:
        print ('Warning: no element in "same_id"')
    try:
        sensitivity_same_var = df[correct & cat_same_var].shape[0] / df[cat_same_var].shape[0]
        print ('sensitivity_same_var = {0:.4%} ({1} / {2})'.format(sensitivity_same_var, df[correct & cat_same_var].shape[0], df[cat_same_var].shape[0]))
        fnr_same_var = 1 - sensitivity_same_var
        print ('fnr_same_var         = {0:.4%} ({1} / {2})'.format(fnr_same_var, df[cat_same_var].shape[0] - df[correct & cat_same_var].shape[0], df[cat_same_var].shape[0]))
    except:
        print ('Warning: no element in "same_var"')
    try:
        sensitivity_diff_id = df[correct & cat_diff_id].shape[0] / df[cat_diff_id].shape[0]
        print ('sensitivity_diff_id  = {0:.4%} ({1} / {2})'.format(sensitivity_diff_id, df[correct & cat_diff_id].shape[0], df[cat_diff_id].shape[0]))
        fnr_diff_id = 1 - sensitivity_diff_id
        print ('fnr_diff_id          = {0:.4%} ({1} / {2})'.format(fnr_diff_id, df[cat_diff_id].shape[0] - df[correct & cat_diff_id].shape[0], df[cat_diff_id].shape[0]))
    except:
        print ('Warning: no element in "diff_id"')
    try:
        sensitivity_diff_var = df[correct & cat_diff_var].shape[0] / df[cat_diff_var].shape[0]
        print ('sensitivity_diff_var = {0:.4%} ({1} / {2})'.format(sensitivity_diff_var, df[correct & cat_diff_var].shape[0], df[cat_diff_var].shape[0]))
        fnr_diff_var = 1 - sensitivity_diff_var
        print ('fnr_diff_var         = {0:.4%} ({1} / {2})'.format(fnr_diff_var, df[cat_diff_var].shape[0] - df[correct & cat_diff_var].shape[0], df[cat_diff_var].shape[0]))
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
        print ('fdr_mapqall          = {0:.4%}'.format(fdr_mapqall))
    except:
        print ('Warning: no read is 0+ mapq')
        fdr_mapqall = 1
    try:
        precision_mapq5plus = df[correct & mapq5plus & v_filter].shape[0] / df[mapq5plus & v_filter].shape[0]
        print ('precision_mapq5plus  = {0:.4%} ({1} / {2})'.format(precision_mapq5plus, df[correct & mapq5plus & v_filter].shape[0], df[mapq5plus & v_filter].shape[0]))
        fdr_mapq5plus = 1 - precision_mapq5plus
        print ('fdr_mapq5plus        = {0:.4%}'.format(fdr_mapq5plus))
    except:
        print ('Warning: no read is 5+ mapq')
        fdr_mapq5plus = 1
    try:
        precision_mapq10plus = df[correct & mapq10plus & v_filter].shape[0] / df[mapq10plus & v_filter].shape[0]
        print ('precision_mapq10plus = {0:.4%} ({1} / {2})'.format(precision_mapq10plus, df[correct & mapq10plus & v_filter].shape[0], df[mapq10plus & v_filter].shape[0]))
        fdr_mapq10plus = 1 - precision_mapq10plus
        print ('fdr_mapq10plus       = {0:.4%}'.format(fdr_mapq10plus))
    except:
        print ('Warning: no read is 10+ mapq')
        fdr_mapq10plus = 1
    try:
        precision_mapq20plus = df[correct & mapq20plus & v_filter].shape[0] / df[mapq20plus & v_filter].shape[0]
        print ('precision_mapq20plus = {0:.4%} ({1} / {2})'.format(precision_mapq20plus, df[correct & mapq20plus & v_filter].shape[0], df[mapq20plus & v_filter].shape[0]))
        fdr_mapq20plus = 1 - precision_mapq20plus
        print ('fdr_mapq20plus       = {0:.4%}'.format(fdr_mapq20plus))
    except:
        print ('Warning: no read is 20+ mapq')
        fdr_mapq20plus = 1
    try:
        precision_mapq30plus = df[correct & mapq30plus & v_filter].shape[0] / df[mapq30plus & v_filter].shape[0]
        print ('precision_mapq30plus = {0:.4%} ({1} / {2})'.format(precision_mapq30plus, df[correct & mapq30plus & v_filter].shape[0], df[mapq30plus & v_filter].shape[0]))
        fdr_mapq30plus = 1 - precision_mapq30plus
        print ('fdr_mapq30plus       = {0:.4%}'.format(fdr_mapq30plus))
    except:
        print ('Warning: no read is 30+ mapq')
        fdr_mapq30plus = 1
    try:
        precision_mapq40plus = df[correct & mapq40plus & v_filter].shape[0] / df[mapq40plus & v_filter].shape[0]
        print ('precision_mapq40plus = {0:.4%} ({1} / {2})'.format(precision_mapq40plus, df[correct & mapq40plus & v_filter].shape[0], df[mapq40plus & v_filter].shape[0]))
        fdr_mapq40plus = 1 - precision_mapq40plus
        print ('fdr_mapq40plus       = {0:.4%}'.format(fdr_mapq40plus))
    except:
        print ('Warning: no read is 40+ mapq')
        fdr_mapq40plus = 1
    
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
    var_reads_fn = args.var_reads
    var_sample_fn = args.var_sample
    personalized = args.personalized
    chrm = args.chrm
    write_wrt_correctness = args.write_wrt_correctness
    mapq_threshold = args.write_wrt_mapq
    step = args.step_size
    read_len = args.read_len
    fn_ref = args.debug_ref
    fn_hapA = args.debug_hapA
    fn_hapB = args.debug_hapB
    
    #: Global variables
    global MAIN_CHRM, ALT_CHRM, MAIN_HAP, ALT_HAP, MAIN_STRAND, ALT_STRAND
    MAIN_STRAND = 'A'
    ALT_STRAND = 'B'
    MAIN_HAP = 'hap' + MAIN_STRAND
    ALT_HAP = 'hap' + ALT_STRAND
    MAIN_CHRM = chrm + MAIN_STRAND
    ALT_CHRM = chrm + ALT_STRAND

    if mapq_threshold:
        write_wrt_mapq(sam_fn, int(mapq_threshold))
        exit()
    
    USE_PREV_IF_POSSIBLE = False
    PLOT_HIST = False
    if USE_PREV_IF_POSSIBLE and os.path.isfile(sam_fn + '-stats.pkl'):
        print ('Read stats from {0}-stats.pkl'.format(sam_fn))
        df = pd.read_pickle(sam_fn + '-stats.pkl')
        
        if PLOT_HIST:
            plot_hist(sam_fn, df)

        print_df_stats(df, threshold, 'all')
        exit()

    global COMPARE_SEQ, HIGHC, TOTALNEAR, REF_G, HAPA_G, HAPB_G, CALL_D_ALT, SIM_D_ALT, CALL_D_ORIG, SIM_D_ORIG
    if (fn_ref != None) and (fn_hapA != None) and (fn_hapB != None):
        COMPARE_SEQ = True
        HIGHC = 0
        TOTALNEAR = 0
        CALL_D_ALT = []
        SIM_D_ALT = []
        CALL_D_ORIG = []
        SIM_D_ORIG = []
        REF_G = read_genome(fn_ref)
        HAPA_G = read_genome(fn_hapA)
        HAPB_G = read_genome(fn_hapB)
    else:
        COMPARE_SEQ = False

    all_indexes = build_all_indexes(
        var_reads_fn=var_reads_fn,
        var_sample_fn=var_sample_fn,
        personalized=personalized
    )
    golden_dic = load_golden_dic(golden_fn, 1)
    results_df = analyze_diploid_indels(
        sam_fn=sam_fn,
        golden_dic=golden_dic,
        threshold=threshold,
        all_indexes=all_indexes,
        personalized=personalized,
        chrm=chrm,
        step=step,
        read_len=read_len,
        write_wrt_correctness=write_wrt_correctness
    )

    print_df_stats(results_df, threshold, 'all')

    if COMPARE_SEQ:
        print ('Num of near alignments =', TOTALNEAR)
        print ('Num of alns have higher score than golden =', HIGHC)

        if HIGHC > 0:
            print ('Avg Lev. dist of called ALT alignments =', sum(CALL_D_ALT)/len(CALL_D_ALT))
            print ('Avg Lev. dist of simulated ALT alignments =', sum(SIM_D_ALT)/len(SIM_D_ALT))
        
        print ('Num of near alignments with higher golden score', TOTALNEAR - HIGHC)
        if TOTALNEAR - HIGHC > 0:
            print ('Avg Lev. dist of called ORIG alignments =', sum(CALL_D_ORIG)/(TOTALNEAR-HIGHC))
            print ('Avg Lev. dist of simulated ORIG alignments =', sum(SIM_D_ORIG)/(TOTALNEAR-HIGHC))
