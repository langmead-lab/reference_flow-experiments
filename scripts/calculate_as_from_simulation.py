'''
This script reads a target sam file and a golden sam file (from simulation)
and summarizes the number of reads with following characteristics:
    1. Alignment score (AS:i, or AS) relationship between alignment
       (target, aln) and simulation (sim)
        * aln > sim
        * aln = sim
        * aln < sim
    2. Carrying a sequencing error or not. This characteristic is 
       calculated using AS_sim, if it's 0 there are no errors added 
       during simulation
        * no seq error
        * has seq error(s)
    3. Correctness. This script doesn't calculated correctness by itself,
       but retrieves the information in the given .sam-stats.pkl file
        * correct
        * incorrect
    4. #TODO Overlapping a unique variant or not.
'''
import pandas as pd
import numpy as np
import argparse
import re
import math
import constants
from analyze_sam import SamInfo, parse_line
from build_erg import VarInfo, read_var
from analyze_diploid_indels import build_index

parser = argparse.ArgumentParser()
parser.add_argument(
    '-n', '--fn_sam',
    help='target sam file'
)
parser.add_argument(
    '-p', '--fn_sam_stats',
    help='target sam-stats.pkl file'
)
parser.add_argument(
    '-g', '--fn_gold',
    help='golden .fq.sam file'
)
parser.add_argument(
    '-vr', '--fn_var_sim',
    help='var file specifying unique variants \
        of the personalized genome'
)
parser.add_argument(
    '-vs', '--fn_var_aln',
    help='var file specifying unique variants \
        of the aligned genome'
)
parser.add_argument(
    '-read_len', '--read_len', type=int, default=100,
    help='read length (100)'
)
args = parser.parse_args()
read_len = args.read_len
fn_gold = args.fn_gold
fn_sam = args.fn_sam
fn_sam_stats = args.fn_sam_stats
fn_var_sim = args.fn_var_sim
fn_var_aln = args.fn_var_aln

def report_stats(df):
    assert sum(df['AS_sim'] > 0) == 0
    #: unaligned
    num_unaligned = sum(df['AS_aln'] == 1)
    num_unaligned_noerr = df[(df['AS_aln'] == 1) & (df['AS_sim'] == 0)].shape[0]
    num_unaligned_haserr = df[(df['AS_aln'] == 1) & (df['AS_sim'] < 0)].shape[0]
    #: AS -- align > simulation
    num_aln_gt_sim = sum(df['AS_aln'] > df['AS_sim']) - num_unaligned
    num_aln_gt_sim_noerr  = df[(df['AS_aln'] > df['AS_sim']) & (df['AS_aln'] != 1) & (df['AS_sim'] == 0)].shape[0]
    num_aln_gt_sim_noerr_true  = df[(df['AS_aln'] > df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == True)  & (df['AS_sim'] == 0)].shape[0]
    num_aln_gt_sim_noerr_false = df[(df['AS_aln'] > df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == False) & (df['AS_sim'] == 0)].shape[0]
    num_aln_gt_sim_haserr = df[(df['AS_aln'] > df['AS_sim']) & (df['AS_aln'] != 1) & (df['AS_sim'] < 0)].shape[0]
    num_aln_gt_sim_haserr_true  = df[(df['AS_aln'] > df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == True)  & (df['AS_sim'] < 0)].shape[0]
    num_aln_gt_sim_haserr_false = df[(df['AS_aln'] > df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == False) & (df['AS_sim'] < 0)].shape[0]
    #: AS -- align == simulation
    num_aln_eq_sim = sum(df['AS_aln'] == df['AS_sim'])
    num_aln_eq_sim_noerr = df[(df['AS_aln'] == df['AS_sim']) & (df['AS_aln'] != 1) & (df['AS_sim'] == 0)].shape[0]
    num_aln_eq_sim_noerr_true  = df[(df['AS_aln'] == df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == True)  & (df['AS_sim'] == 0)].shape[0]
    num_aln_eq_sim_noerr_false = df[(df['AS_aln'] == df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == False) & (df['AS_sim'] == 0)].shape[0]
    num_aln_eq_sim_haserr = df[(df['AS_aln'] == df['AS_sim']) & (df['AS_aln'] != 1) & (df['AS_sim'] < 0)].shape[0]
    num_aln_eq_sim_haserr_true  = df[(df['AS_aln'] == df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == True) & (df['AS_sim'] < 0)].shape[0]
    num_aln_eq_sim_haserr_false = df[(df['AS_aln'] == df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == False) & (df['AS_sim'] < 0)].shape[0]
    #: AS -- align < simulation
    num_aln_lt_sim = sum(df['AS_aln'] < df['AS_sim'])
    num_aln_lt_sim_noerr = df[(df['AS_aln'] < df['AS_sim']) & (df['AS_aln'] != 1) & (df['AS_sim'] == 0)].shape[0]
    num_aln_lt_sim_noerr_true  = df[(df['AS_aln'] < df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == True)  & (df['AS_sim'] == 0)].shape[0]
    num_aln_lt_sim_noerr_false = df[(df['AS_aln'] < df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == False) & (df['AS_sim'] == 0)].shape[0]
    num_aln_lt_sim_haserr = df[(df['AS_aln'] < df['AS_sim']) & (df['AS_aln'] != 1) & (df['AS_sim'] < 0)].shape[0]
    num_aln_lt_sim_haserr_true  = df[(df['AS_aln'] < df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == True) & (df['AS_sim'] < 0)].shape[0]
    num_aln_lt_sim_haserr_false = df[(df['AS_aln'] < df['AS_sim']) & (df['AS_aln'] != 1) & (df['correctness'] == False) & (df['AS_sim'] < 0)].shape[0]

    print ('Num. AS_aln >  AS_sim   = {0:8d} ({1:.4%})'.format(num_aln_gt_sim, num_aln_gt_sim/df.shape[0]))
    if num_aln_gt_sim > 0:
        print ('    Num. no  seq. error =       {0:8d} ({1:.4%})'.format(num_aln_gt_sim_noerr, num_aln_gt_sim_noerr/num_aln_gt_sim))
        if num_aln_gt_sim_noerr > 0:
            print ('        Num. correct    =               {0:8d} ({1:.4%})'.format(num_aln_gt_sim_noerr_true, num_aln_gt_sim_noerr_true/num_aln_gt_sim_noerr))
            print ('        Num. incorrect  =               {0:8d} ({1:.4%})'.format(num_aln_gt_sim_noerr_false, num_aln_gt_sim_noerr_false/num_aln_gt_sim_noerr))
        print ('    Num. has seq. error =       {0:8d} ({1:.4%})'.format(num_aln_gt_sim_haserr, num_aln_gt_sim_haserr/num_aln_gt_sim))
        if num_aln_gt_sim_haserr > 0:
            print ('        Num. correct    =               {0:8d} ({1:.4%})'.format(num_aln_gt_sim_haserr_true, num_aln_gt_sim_haserr_true/num_aln_gt_sim_haserr))
            print ('        Num. incorrect  =               {0:8d} ({1:.4%})'.format(num_aln_gt_sim_haserr_false, num_aln_gt_sim_haserr_false/num_aln_gt_sim_haserr))

    print ('Num. AS_aln == AS_sim   = {0:8d} ({1:.4%})'.format(num_aln_eq_sim, num_aln_eq_sim/df.shape[0]))
    if num_aln_eq_sim > 0:
        print ('    Num. no  seq. error =       {0:8d} ({1:.4%})'.format(num_aln_eq_sim_noerr, num_aln_eq_sim_noerr/num_aln_eq_sim))
        if num_aln_eq_sim_noerr > 0:
            print ('        Num. correct    =               {0:8d} ({1:.4%})'.format(num_aln_eq_sim_noerr_true, num_aln_eq_sim_noerr_true/num_aln_eq_sim_noerr))
            print ('        Num. incorrect  =               {0:8d} ({1:.4%})'.format(num_aln_eq_sim_noerr_false, num_aln_eq_sim_noerr_false/num_aln_eq_sim_noerr))
        print ('    Num. has seq. error =       {0:8d} ({1:.4%})'.format(num_aln_eq_sim_haserr, num_aln_eq_sim_haserr/num_aln_eq_sim))
        if num_aln_eq_sim_haserr > 0:
            print ('        Num. correct    =               {0:8d} ({1:.4%})'.format(num_aln_eq_sim_haserr_true, num_aln_eq_sim_haserr_true/num_aln_eq_sim_haserr))
            print ('        Num. incorrect  =               {0:8d} ({1:.4%})'.format(num_aln_eq_sim_haserr_false, num_aln_eq_sim_haserr_false/num_aln_eq_sim_haserr))

    print ('Num. AS_aln <  AS_sim   = {0:8d} ({1:.4%})'.format(num_aln_lt_sim, num_aln_lt_sim/df.shape[0]))
    if num_aln_lt_sim > 0:
        print ('    Num. no  seq. error =       {0:8d} ({1:.4%})'.format(num_aln_lt_sim_noerr, num_aln_lt_sim_noerr/num_aln_lt_sim))
        if num_aln_lt_sim_noerr > 0:
            print ('        Num. correct    =               {0:8d} ({1:.4%})'.format(num_aln_lt_sim_noerr_true, num_aln_lt_sim_noerr_true/num_aln_lt_sim_noerr))
            print ('        Num. incorrect  =               {0:8d} ({1:.4%})'.format(num_aln_lt_sim_noerr_false, num_aln_lt_sim_noerr_false/num_aln_lt_sim_noerr))
        print ('    Num. has seq. error =       {0:8d} ({1:.4%})'.format(num_aln_lt_sim_haserr, num_aln_lt_sim_haserr/num_aln_lt_sim))
        if num_aln_lt_sim_haserr > 0:
            print ('        Num. correct    =               {0:8d} ({1:.4%})'.format(num_aln_lt_sim_haserr_true, num_aln_lt_sim_haserr_true/num_aln_lt_sim_haserr))
            print ('        Num. incorrect  =               {0:8d} ({1:.4%})'.format(num_aln_lt_sim_haserr_false, num_aln_lt_sim_haserr_false/num_aln_lt_sim_haserr))

    print ('Num. unaligned          = {0:8d} ({1:.4%})'.format(num_unaligned, num_unaligned/df.shape[0]))
    if num_unaligned > 0:
        print ('    Num. no  seq. error =       {0:8d} ({1:.4%})'.format(num_unaligned_noerr, num_unaligned_noerr/num_unaligned))
        print ('    Num. has seq. error =       {0:8d} ({1:.4%})'.format(num_unaligned_haserr, num_unaligned_haserr/num_unaligned))

#: quoting = csv.QUOTE_NONE (3)
#: This is important! If not set, there would be issues in the QUAL field 
#: and resuls in incorrect reading
df = pd.read_csv(fn_gold, sep='\t', header=None, quoting=3)
df.columns = \
    ['QNAME', 'FLAG', 'RNAME_sim', 'POS_sim', 'MAPQ', 
    'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 
    'NM', 'MD', 'XE', 'XS', 'XI']

#: processes MD string
split_md = []
for i in range(df.shape[0]):
    split_md.append(re.findall('([0-9]+|[a-zA-Z]+)', df['MD'].iloc[i].split(':')[-1]))
df['sMD'] = split_md

#: processes CIGAR string
split_cigar = []
for i in range(df.shape[0]):
    split_cigar.append(re.findall('([0-9]+[a-zA-Z]+)', df['CIGAR'].iloc[i]))

reduced_cigar = []
for i in range(df.shape[0]):
    assert len(df['SEQ'].iloc[i]) == read_len
    assert len(df['QUAL'].iloc[i]) == read_len

    ci = split_cigar[i]
    for c in ci:
        if c.count('P') > 0:
            ci.remove(c)
    
    ci_flag = True
    while ci_flag:
        ci_flag = False
        for j, c in enumerate(ci):
            if j == len(ci)-1:
                break
            elif ci[j+1][-1] == c[-1]:
                ci_flag = True
                count = int(c[:-1]) + int(ci[j+1][:-1])
                ci[j] = str(count) + c[-1]
                del ci[j+1]

    #: checks if length of cigar (I+M) matches read length
    len_ci = 0
    for c in ci:
        assert c[-1] in ['I','D','M']
        if c[-1] != 'D':
            len_ci += int(c[:-1])
    assert len_ci == read_len
    reduced_cigar.append(ci)
df['rCIGAR'] = reduced_cigar

#: calculates AS:i from simulation
list_as = [None] * df.shape[0]
for i in range(df.shape[0]):
    #: reads with no error
    if int(df['XE'].iloc[i].split(':')[-1]) == 0:
        list_as[i] = 0
    else:
        score = 0
        vec_aln = []
        for c in df['rCIGAR'].iloc[i]:
            #: "0" = match/mismatch
            if c[-1] == 'M':
                for z in range(int(c[:-1])):
                    vec_aln.append(0)
            elif c[-1] in ['I', 'D']:
                len_indel = int(c[:-1])
                score -= (5 + 3 * len_indel)
                #: "1" = insertion
                if c[-1] == 'I':
                    for z in range(int(c[:-1])):
                        vec_aln.append(1)
                #: "2" = deletion
                else:
                    for z in range(int(c[:-1])):
                        vec_aln.append(2)
        current_pos = 0
        current_count_del = 0
        current_cigar_pos = current_pos + current_count_del
        for c in df['sMD'].iloc[i]:
            if c.isalpha():
                #: checks INS
                overlapping_ins = vec_aln[current_cigar_pos: current_cigar_pos + len(c)].count(1)
                if overlapping_ins > 0:
                    while vec_aln[current_cigar_pos: current_cigar_pos + overlapping_ins + len(c)].count(1) > overlapping_ins:
                        overlapping_ins = vec_aln[current_cigar_pos: current_cigar_pos + overlapping_ins + len(c)].count(1)
                current_pos += overlapping_ins
                for a in c:
                    if vec_aln[current_cigar_pos] == 2:
                        current_count_del += 1
                        # print ('D at {}'.format(current_pos))
                    else:
                        qual = df['QUAL'].iloc[i][current_pos]
                        qual = ord(qual) - 33
                        score -= (2 + math.floor(4 * min(qual, 40)/40))
                        current_pos += 1
                    current_cigar_pos = current_pos + current_count_del
            else:
                # current_pos -= overlapping_ins
                #: checks INS
                overlapping_ins = vec_aln[current_cigar_pos: current_cigar_pos + int(c) + 1].count(1)
                if overlapping_ins > 0:
                    while vec_aln[current_cigar_pos: current_cigar_pos + overlapping_ins + int(c) + 1].count(1) > overlapping_ins:
                        overlapping_ins = vec_aln[current_cigar_pos: current_cigar_pos + overlapping_ins + int(c) + 1].count(1)
                current_pos += int(c) + overlapping_ins
                # current_pos += int(c)
            current_cigar_pos = current_pos + current_count_del
        if current_pos != read_len:
            print ('{0} != {1} at {2}'.format(current_cigar_pos, read_len, i))
        assert current_pos == read_len
        list_as[i] = score
df['AS_sim'] = list_as

#: reads target sam file
f_sam = open(fn_sam, 'r')
list_name = []
list_as_sam = []
list_pos = []
list_rname = []
for line in f_sam:
    name, info = parse_line(line)
    #: headers
    if name == 'header':
        continue
    list_name.append(name)
    list_as_sam.append(info.score)
    list_rname.append(info.chrm)
    list_pos.append(info.pos)

df_sam = pd.DataFrame(columns=['QNAME', 'AS_aln'])
df_sam['QNAME'] = list_name
df_sam['RNAME_aln'] = list_rname
df_sam['POS_aln'] = list_pos
df_sam['AS_aln'] = list_as_sam

df_merge = df.merge(df_sam, on='QNAME')
# df_merge.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NM', 'MD', 'XE', 'XS', 'XI', 'sMD', 'rCIGAR', 'AS_sim', 'AS_aln'] 

df_results = pd.read_pickle(fn_sam_stats)
list_correctness = [1 if i >=0 and i <=10 else 0 for i in df_results['dist']]
df_results['correctness'] = list_correctness
df_results.columns = ['QNAME', 'dist', 'mapq', 'numvar', 'category', 'correctness']
df_results_tmp = df_results[['QNAME', 'correctness']]

df_merge2 = df_merge.merge(df_results_tmp, on='QNAME')
report_stats(df_merge2)


def count_overlapping_vars_pos(
    chrom,
    pos,
    main_index,
    alt_index,
    main_chrom,
    alt_chrom
):
    '''
    For an alignment, count the number of overlapping variants.
    '''
    num_var = 0   
    for i in range(pos, pos + constants.READ_LEN):
        if chrom == main_chrom:
            if main_index.get(i) != None:
                num_var += 1
        elif chrom == alt_chrom:
            if alt_index.get(i) != None:
                num_var += 1
        #: unaligned
        elif chrom == '*':
            num_var = -1
            return num_var
        else:
            print ('Error: unexpected chrm', chrom)
            exit()
    return num_var

constants.set_chrom('21')
constants.set_read_len(read_len)

list_var_sim = read_var(
    fn_var_sim,
    remove_conflict=True,
    remove_homo_alt=False
)
main_index_sim, alt_index_sim = build_index(list_var_sim)

list_var_aln = read_var(
    fn_var_aln,
    remove_conflict=True,
    remove_homo_alt=False
)
main_index_aln, _ = build_index(list_var_aln)

list_num_overlapping_var_sim = []
list_num_overlapping_var_aln = []
for i in range(df_merge2.shape[0]):
    num_ov_sim = count_overlapping_vars_pos(
        df_merge2['RNAME_sim'].iloc[i],
        df_merge2['POS_sim'].iloc[i], 
        main_index_sim,
        alt_index_sim,
        constants.MAIN_CHROM,
        constants.ALT_CHROM
    )
    list_num_overlapping_var_sim.append(num_ov_sim)
    num_ov_aln = count_overlapping_vars_pos(
        df_merge2['RNAME_aln'].iloc[i],
        df_merge2['POS_aln'].iloc[i], 
        main_index_aln,
        _,
        constants.CHROM,
        constants.CHROM
    )
    list_num_overlapping_var_aln.append(num_ov_aln)
df_merge2['NumOV_sim'] = list_num_overlapping_var_sim
df_merge2['NumOV_aln'] = list_num_overlapping_var_aln

from scipy.stats import pearsonr
pearsonr(df_merge2['NumOV_sim'], df_merge2['correctness'])
pearsonr(df_merge2['NumOV_aln'], df_merge2['correctness'])