'''
Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse, math, sys, os
import pandas as pd
# from analyze_sam import SamInfo, parse_line, load_golden_dic, Summary
from lib_compare_sam import SamInfo, parse_line, load_golden_dic, Summary, print_df_stats
from build_erg import read_var, read_genome
import constants

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
        '-c', '--chrom',
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
        '--write_wrt_correctness',
        action='store_true',
        #default=None,
        help='(int) If set, writes two files recording correct/incorrect alignments respectively. The output files use target sam prefix [None].'
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

def build_index(var_list):
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
        if v.strand == constants.MAIN_STRAND:
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

def build_offset_index(var_list):
    # , step, MAIN_STRAND, ALT_STRAND):
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
        if v.strand == constants.MAIN_STRAND:
            main_pos = v.alt_pos
            alt_pos = v.ref_pos + v.cor_offset
            main_offset = -v.offset + v.cor_offset
            alt_offset = v.offset - v.cor_offset
            if v.ref_pos == tmp_v:
                if SHOW_DUP_WARN:
                    print ('Warning: duplicated variant', v.line)
            tmp_v = v.ref_pos
        elif v.strand == constants.ALT_STRAND:
            main_pos = v.ref_pos + v.cor_offset
            alt_pos = v.alt_pos
            main_offset = v.offset - v.cor_offset
            alt_offset = -v.offset + v.cor_offset
        else:
            print ('Error: unspecified strand', v.strand)
            exit()
        
        i_main = math.ceil(main_pos / constants.STEP)
        while i_main >= len(main_offset_index):
            main_offset_index.append(main_offset_index[len(main_offset_index) - 1])
        main_offset_index[i_main] = main_offset
        i_alt = math.ceil(alt_pos / constants.STEP)
        while i_alt >= len(alt_offset_index):
            alt_offset_index.append(alt_offset_index[len(alt_offset_index) - 1])
        alt_offset_index[i_alt] = alt_offset
    
    return main_offset_index, alt_offset_index

def build_offset_index_ref(var_list):
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
        idx = math.ceil(v.alt_pos / constants.STEP)
        if v.strand == constants.MAIN_STRAND:
            while idx >= len(alt1_offset_index):
                alt1_offset_index.append(alt1_offset_index[len(alt1_offset_index) - 1])
            alt1_offset_index[idx] = offset
        elif v.strand == constants.ALT_STRAND:
            while idx >= len(alt2_offset_index):
                alt2_offset_index.append(alt2_offset_index[len(alt2_offset_index) - 1])
            alt2_offset_index[idx] = offset
        else:
            print ('Error: unspecified strand', v.strand)
            exit()
    
    return alt1_offset_index, alt2_offset_index

# TODO
def print_aln_within_distance(name, reads_offsets, sample_offsets, info, g_info, threshold, read_len, COMPARE_SEQ):
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
                print ('ORIG1  (%10d) = %s' % (g_info.pos - reads_offsets[0], REF_G[g_info.pos - reads_offsets[0] : g_info.pos - reads_offsets[0] + 80]))
                if reads_offsets[0] != reads_offsets[1]:
                    print ('ORIG2  (%10d) = %s' % (g_info.pos - reads_offsets[1], REF_G[g_info.pos - reads_offsets[1] : g_info.pos - reads_offsets[1] + 80]))
                print ('PERSON (#%9d) = %s' % (g_info.pos, HAPA_G[g_info.pos : g_info.pos + 80]))
        return

def compare_sam_info(
        info,
        ginfo,
        threshold,
        sample_offsets=[0],
        reads_offsets=[0],
        ignore_chrom=False
    ):
    '''
    Inputs:
        info:
            info from alignment
        ginfo:
            info from simulation profile (golden)
        offset:
            positiontal offset
        ignore_chrom:
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
    if (ignore_chrom is False) and (info.chrom != ginfo.chrom):
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

def diploid_compare(
    info, 
    g_info,
    name, 
    threshold, 
    dip_flag, 
    COMPARE_SEQ,
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
            i_low = int(info.pos / constants.STEP)
            i_high = math.ceil(info.pos / constants.STEP)
            if i_low >= len(sample_main_offset_index):
                sample_offset_low = sample_main_offset_index[len(sample_main_offset_index) - 1]
            else:
                sample_offset_low = sample_main_offset_index[i_low]
            if i_high >= len(sample_main_offset_index):
                sample_offset_high = sample_main_offset_index[len(sample_main_offset_index) - 1]
            else:
                sample_offset_high = sample_main_offset_index[i_high]
            sample_offsets = [sample_offset_low, sample_offset_high]
    elif dip_flag in ['same_id', 'same_var', 'diff_id', 'diff_var']:
        i_low = int(info.pos / constants.STEP)
        i_high = math.ceil(info.pos / constants.STEP)
        if info.chrom == constants.MAIN_CHROM:
            if i_low >= len(sample_main_offset_index):
                sample_offset_low = sample_main_offset_index[len(sample_main_offset_index) - 1]
            else:
                sample_offset_low = sample_main_offset_index[i_low]
            if i_high >= len(sample_main_offset_index):
                sample_offset_high = sample_main_offset_index[len(sample_main_offset_index) - 1]
            else:
                sample_offset_high = sample_main_offset_index[i_high]
        elif info.chrom == constants.ALT_CHROM:
            if i_low >= len(sample_alt_offset_index):
                sample_offset_low = sample_alt_offset_index[len(sample_alt_offset_index) - 1]
            else:
                sample_offset_low = sample_alt_offset_index[i_low]
            if i_high >= len(sample_alt_offset_index):
                sample_offset_high = sample_alt_offset_index[len(sample_alt_offset_index) - 1]
            else:
                sample_offset_high = sample_alt_offset_index[i_high]
        else:
            print ('Error: invalid chrom', info.chrom, constants.MAIN_CHROM, constants.ALT_CHROM)
            exit()
        sample_offsets = [sample_offset_low, sample_offset_high]
    else:
        print ('Error: undistinguished dip_flag: %s' % dip_flag)
        return False
    
    i_low = int(g_info.pos / constants.STEP)
    i_high = math.ceil(g_info.pos / constants.STEP)
    reads_offsets = []
    #: check hapA
    if name.find(constants.MAIN_HAP) > 0:
        if i_low >= len(reads_main_offset_index):
            reads_offsets.append(reads_main_offset_index[len(reads_main_offset_index) - 1])
        else:
            reads_offsets.append(reads_main_offset_index[i_low])
        if i_high >= len(reads_main_offset_index):
            reads_offsets.append(reads_main_offset_index[len(reads_main_offset_index) - 1])
        else:
            reads_offsets.append(reads_main_offset_index[i_high])
    #: check hapB
    elif name.find(constants.ALT_HAP) > 0:
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
        ignore_chrom=True
    )

    if (dist < 0 or dist > threshold) and COMPARE_SEQ:
        print_aln_within_distance(
            name=name,
            reads_offsets=reads_offsets,
            sample_offsets=sample_offsets,
            info=info,
            g_info=g_info,
            threshold=1000,
            read_len=read_len,
            COMPARE_SEQ=COMPARE_SEQ
        )

    return dist

def count_overlapping_vars(
    name,
    info,
    g_info,
    main_index,
    alt_index
):
    '''
    For an alignment, count the number of overlapping variants.
    The count is based on simulated position 
    (look up golden dictionary).
    '''
    num_var = 0   
    for i in range(g_info.pos, g_info.pos + constants.READ_LEN):
        if g_info.chrom == constants.MAIN_CHROM or g_info.chrom == constants.CHROM or g_info.chrom == 'chr' + constants.MAIN_CHROM or g_info.chrom == 'chr' + constants.CHROM:
        # if g_info.chrom == constants.MAIN_CHROM or g_info.chrom == constants.CHROM:
        # if g_info.chrom == constants.MAIN_CHROM:
            if main_index.get(i) != None:
                num_var += 1
        elif g_info.chrom == constants.ALT_CHROM or g_info.chrom == 'chr' + constants.ALT_CHROM:
        # elif g_info.chrom == constants.ALT_CHROM:
            if alt_index.get(i) != None:
                num_var += 1
        else:
            print ('Error: unexpected chrom', info.chrom)
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
    var_reads_list = read_var(
        var_reads_fn,
        remove_conflict=True,
        remove_homo_alt=False
    )
    main_index, alt_index = build_index(var_reads_list)
    reads_main_offset_index, reads_alt_offset_index = build_offset_index_ref(var_reads_list)
    #: diploid personalized ref
    if personalized == 2:
        var_sample_list = read_var(
            var_sample_fn,
            remove_conflict=True,
            remove_homo_alt=False
        )
        sample_main_offset_index, sample_alt_offset_index = build_offset_index_ref(var_sample_list)
    #: standard ref seq
    elif personalized == 0:
        #: major allele reference with indels
        if var_sample_fn != None:
            var_sample_list = read_var(
                var_sample_fn,
                remove_conflict=True,
                remove_homo_alt=False
            )
            sample_main_offset_index, _ = build_offset_index_ref(var_sample_list)
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
    chrom,
    step,
    read_len,
    write_wrt_correctness,
    COMPARE_SEQ
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
        #name, info = parse_line(line, erg=True)
        # single-end
        name, info = parse_line(line, erg=True, mason2=True, score=COMPARE_SEQ)
        # paired-end
        # name, info = parse_line(line, erg=True, mason2=False, score=COMPARE_SEQ)
        #: headers
        if name == 'header':
            continue
        summary.add_one()

        # first segment: 0
        # second segment: 1
        # g_info = golden_dic[name][info.is_first_seg() ^ 1]
        g_info = golden_dic[name][info.is_first_seg()]

        #: counts the number of overlapping variants for all aligned reads
        num_var = count_overlapping_vars(
            name=name,
            info=info,
            g_info=g_info,
            main_index=main_index,
            alt_index=alt_index
        )

        if info.is_unaligned():
            dist = -3
            flag = 'unaligned'
        #: alignment against personalized genomes
        elif personalized == 2:
            #: aligned to incorrect haplotype
            name_chrom_mismatch = (
                (name.find(constants.MAIN_HAP) > 0 and info.chrom != constants.MAIN_CHROM) or 
                (name.find(constants.ALT_HAP) > 0 and info.chrom != constants.ALT_CHROM)
            )
            if name_chrom_mismatch:
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
                # g_info=golden_dic[name],
                g_info=g_info,
                name=name, 
                threshold=threshold, 
                dip_flag=flag, 
                reads_main_offset_index = reads_main_offset_index, 
                reads_alt_offset_index = reads_alt_offset_index,
                sample_main_offset_index = sample_main_offset_index,
                sample_alt_offset_index = sample_alt_offset_index,
                COMPARE_SEQ=COMPARE_SEQ
            )
        #: alignment against standard ref (and ERG)
        elif personalized == 0:
            flag = 'same_strand_ref'
            dist = diploid_compare(
                info=info, 
                # g_info=golden_dic[name],
                g_info=g_info,
                name=name, 
                threshold=threshold, 
                dip_flag=flag, 
                reads_main_offset_index = reads_main_offset_index, 
                reads_alt_offset_index = reads_alt_offset_index,
                sample_main_offset_index = sample_main_offset_index,
                sample_alt_offset_index = sample_alt_offset_index,
                COMPARE_SEQ=COMPARE_SEQ
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

    summary.show_summary(has_answer=True)
    sam_f.close()

    results_df = pd.DataFrame(results, columns=['name', 'dist', 'mapq', 'numvar', 'category'])
    results_df.to_pickle(sam_fn + '-stats.pkl')

    return results_df

if __name__ == '__main__':
    args = parse_args()
    sam_fn = args.sam
    golden_fn = args.golden
    threshold = args.threshold
    var_reads_fn = args.var_reads
    var_sample_fn = args.var_sample
    personalized = args.personalized
    chrom = args.chrom
    write_wrt_correctness = args.write_wrt_correctness
    step = args.step_size
    read_len = args.read_len
    fn_ref = args.debug_ref
    fn_hapA = args.debug_hapA
    fn_hapB = args.debug_hapB
    
    constants.set_chrom(chrom)
    constants.set_step(step)
    constants.set_read_len(read_len)

    #USE_PREV_IF_POSSIBLE = False
    USE_PREV_IF_POSSIBLE = True
    if USE_PREV_IF_POSSIBLE and os.path.isfile(sam_fn + '-stats.pkl') and write_wrt_correctness == None:
        print ('Read stats from {0}-stats.pkl'.format(sam_fn))
        df = pd.read_pickle(sam_fn + '-stats.pkl')
        
        print_df_stats(df, threshold, 'all')
        exit()

    # global COMPARE_SEQ, 
    global HIGHC, TOTALNEAR, REF_G, HAPA_G, HAPB_G, CALL_D_ALT, SIM_D_ALT, CALL_D_ORIG, SIM_D_ORIG
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
    golden_dic = load_golden_dic(golden_fn)
    results_df = analyze_diploid_indels(
        sam_fn=sam_fn,
        golden_dic=golden_dic,
        threshold=threshold,
        all_indexes=all_indexes,
        personalized=personalized,
        chrom=chrom,
        step=step,
        read_len=read_len,
        write_wrt_correctness=write_wrt_correctness,
        COMPARE_SEQ=COMPARE_SEQ
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
