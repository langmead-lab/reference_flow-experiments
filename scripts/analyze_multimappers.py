import pandas as pd
import argparse
from analyze_sam import SamInfo, parse_line, load_golden_dic
from analyze_diploid_indels import build_index, count_overlapping_vars, diploid_compare, build_all_indexes
from build_erg import read_var

def check_accuracy_by_dist(dist, threshold):
    if dist >= 0 and dist <= threshold:
        return True
    return False

parser = argparse.ArgumentParser()
parser.add_argument(
    '-stats', '--stats',
    help='target stats file (-stats.pkl)'
)
parser.add_argument(
    '-sam', '--sam',
    help='target sam file'
)
parser.add_argument(
    '-g', '--fn_golden',
    help='ground truth sam file'
)
parser.add_argument(
    '-v', '--var',
    help='.var file specifying variants included in the genome'
)
parser.add_argument(
    '-read_len', '--read_len', type=int, default=100,
    help='read length (100)'
)
parser.add_argument(
    '-t', '--threshold', type=int, default=10,
    help='threshold for comparison (10)'
)
args = parser.parse_args()
fn_stats = args.stats
fn_golden = args.fn_golden
fn_sam = args.sam
fn_var = args.var
read_len = args.read_len
threshold = args.threshold
df = pd.read_pickle(fn_stats)

chrm='21'
MAIN_STRAND = 'A'
ALT_STRAND = 'B'
MAIN_HAP = 'hap' + MAIN_STRAND
ALT_HAP = 'hap' + ALT_STRAND
MAIN_CHRM = chrm + MAIN_STRAND
ALT_CHRM = chrm + ALT_STRAND
COMPARE_SEQ = False

dict_correct = {}
dict_incorrect = {}

for i, name in enumerate(df['name']):
    #: ignores if any alignment for a read is correct
    if name in dict_correct:
        continue
    else:
        #: finds a correct alignment
        if check_accuracy_by_dist(df['dist'][i], threshold=threshold):
            dict_correct[name] = df.iloc[i]
            if name in dict_incorrect:
                dict_incorrect.pop(name, None)
        else:
            if name in dict_incorrect:
                dict_incorrect[name].append(df.iloc[i])
            else:
                dict_incorrect[name] = df.iloc[i]

num_correct = len(dict_correct)
num_incorrect = len(dict_incorrect)
num_total = len(df.groupby('name'))
print ('Number of correct reads = {0}'.format(num_correct))
print ('Number of incorrect reads = {0}'.format(num_incorrect))
print ('Number of total reads = {0}'.format(num_total))
assert num_total == (num_correct + num_incorrect)

print ('Sensitivity = {0:.4%} ( {1} / {2} )'.format(num_correct/num_total, num_correct, num_total))

#: reads variants and builds index
var_reads_list = read_var(fn_var, remove_conflict=True, remove_coexist=False)
all_indexes = build_all_indexes(
    var_reads_fn=fn_var,
    var_sample_fn=fn_var,
    personalized=2,
    step=1000,
    MAIN_STRAND=MAIN_STRAND,
    ALT_STRAND=ALT_STRAND
)
main_index, alt_index, reads_main_offset_index, reads_alt_offset_index, sample_main_offset_index, sample_alt_offset_index = all_indexes

dict_sam_incorrect = {}
dict_sam_correct = {}
f_sam = open(fn_sam, 'r')
for line in f_sam:
    name, info = parse_line(line, erg=True)
    #: headers
    if name == 'header':
        continue
    
    #: aligned incorrectly
    if name in dict_incorrect:
        if name in dict_sam_incorrect:
            dict_sam_incorrect[name].append(info)
        else:
            dict_sam_incorrect[name] = [info]
        continue
    
    #: one of the multi-mappers is correct
    if name in dict_sam_correct:
        score = dict_sam_correct[name][0].score
        if info.score >= score:
            dict_sam_correct[name].append(info)
    else:
        dict_sam_correct[name] = [info]

golden_dic = load_golden_dic(fn_golden, 1)

dict_overlapping_var_correct = {}
for key in list(dict_sam_correct.keys()):
    if len(dict_sam_correct[key]) == 1:
        continue
    var_flag = False
    tmp_info = []
    for i, info in enumerate(dict_sam_correct[key]):
        num_var = count_overlapping_vars(
            name = key,
            info = info,
            g_info = info,
            main_index = main_index,
            alt_index = alt_index,
            MAIN_CHRM = MAIN_CHRM,
            ALT_CHRM = ALT_CHRM,
            read_len = read_len)
        if num_var > 0:
            var_flag = True
        dist = diploid_compare(
            info = info, 
            g_info = info,
            name = key, 
            threshold = threshold, 
            dip_flag = 'same_id', 
            step = 1000,
            COMPARE_SEQ=COMPARE_SEQ,
            reads_main_offset_index = reads_main_offset_index, 
            reads_alt_offset_index = reads_alt_offset_index,
            sample_main_offset_index = sample_main_offset_index,
            sample_alt_offset_index = sample_alt_offset_index,
            MAIN_CHRM=MAIN_CHRM,
            ALT_CHRM=ALT_CHRM,
            MAIN_HAP=MAIN_HAP,
            ALT_HAP=ALT_HAP
        )
        if dist >= 0 and dist <= threshold:
            tmp = '{0}  {1:10d}  {2}  {3}  {4}  {5}  {6}'.format(key, info.pos, info.chrm, info.score, num_var, dist, '*')
        else:
            tmp = '{0}  {1:10d}  {2}  {3}  {4}  {5}  {6}'.format(key, info.pos, info.chrm, info.score, num_var, dist, '')
        tmp_info.append(tmp)
    if var_flag:
        for t in tmp_info:
            print (t)
        # print (golden_dic[key].pos)
        # input ()