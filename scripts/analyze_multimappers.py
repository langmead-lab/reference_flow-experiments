import pandas as pd
import argparse
from analyze_sam import SamInfo, parse_line, load_golden_dic
from analyze_diploid_indels import build_index, diploid_compare, build_all_indexes
from build_erg import read_var

def check_accuracy_by_dist(dist, threshold):
    if dist >= 0 and dist <= threshold:
        return True
    return False

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


parser = argparse.ArgumentParser()
parser.add_argument(
    '-stats', '--fn_stats',
    help='target stats file (-stats.pkl)'
)
parser.add_argument(
    '-sam', '--fn_sam',
    help='target sam file'
)
parser.add_argument(
    '-g', '--fn_golden',
    help='ground truth sam file'
)
# parser.add_argument(
#     '-v', '--var',
#     help='.var file specifying variants included in the genome'
# )
parser.add_argument(
    '-vr', '--fn_var_reads',
    help='.var file specifying variants included in the genome'
)
parser.add_argument(
    '-vs', '--fn_var_sample',
    help='.var file specifying variants included in the genome'
)
parser.add_argument(
    '-o', '--fn_out',
    help='name of output log file'
)
parser.add_argument(
    '-read_len', '--read_len', type=int, default=100,
    help='read length (100)'
)
parser.add_argument(
    '-step', '--step', type=int, default=1000,
    help='step size (1000)'
)
parser.add_argument(
    '-t', '--threshold', type=int, default=10,
    help='threshold for comparison (10)'
)
args = parser.parse_args()
fn_stats = args.fn_stats
fn_golden = args.fn_golden
fn_out = args.fn_out
fn_sam = args.fn_sam
fn_var_sample = args.fn_var_sample
fn_var_reads = args.fn_var_reads
# fn_var = args.var
read_len = args.read_len
step = args.step
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
# var_reads_list = read_var(fn_var, remove_conflict=True, remove_coexist=False)
all_indexes = build_all_indexes(
    var_reads_fn=fn_var_reads,
    var_sample_fn=fn_var_sample,
    personalized=2,
    step=step,
    MAIN_STRAND=MAIN_STRAND,
    ALT_STRAND=ALT_STRAND
)

# main_index, alt_index, reads_main_offset_index, reads_alt_offset_index, sample_main_offset_index, sample_alt_offset_index = all_indexes
_, _, reads_main_offset_index, reads_alt_offset_index, sample_main_offset_index, sample_alt_offset_index = all_indexes
var_reads_list = read_var(fn_var_sample, remove_conflict=True, remove_coexist=False)
main_index, alt_index = build_index(var_reads_list, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)


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
        dict_sam_correct[name].append(info)
        # if info.score >= score:
        #     dict_sam_correct[name].append(info)
    else:
        dict_sam_correct[name] = [info]

golden_dic = load_golden_dic(fn_golden, 1)

df_out = pd.DataFrame(
    columns=['name', 'pos', 'chrom', 'score', 'num_var', 'dist', 'correctness', 'tie_broken']
)
dict_overlapping_var_correct = {}
for key in list(dict_sam_correct.keys()):
    if len(dict_sam_correct[key]) == 1:
        continue
    var_flag = False
    df_tmp = pd.DataFrame(
        columns=['name', 'pos', 'chrom', 'score', 'num_var', 'dist', 'correctness', 'tie_broken']
    )
    for i, info in enumerate(dict_sam_correct[key]):
        #: counts the number of overlapping variants for each aligned region
        num_var = count_overlapping_vars(
            name = key,
            info = info,
            g_info = info,
            # g_info = golden_dic[key],
            main_index = main_index,
            alt_index = alt_index,
            MAIN_CHRM = MAIN_CHRM,
            ALT_CHRM = ALT_CHRM,
            read_len = read_len)
        if num_var > 0:
            var_flag = True
        dist = diploid_compare(
            info = info, 
            g_info = golden_dic[key],
            name = key, 
            threshold = threshold, 
            dip_flag = 'same_id', 
            step = step,
            COMPARE_SEQ = COMPARE_SEQ,
            reads_main_offset_index = reads_main_offset_index, 
            reads_alt_offset_index = reads_alt_offset_index,
            sample_main_offset_index = sample_main_offset_index,
            sample_alt_offset_index = sample_alt_offset_index,
            MAIN_CHRM = MAIN_CHRM,
            ALT_CHRM = ALT_CHRM,
            MAIN_HAP = MAIN_HAP,
            ALT_HAP = ALT_HAP
        )
        if dist >= 0 and dist <= threshold:
            comp = 1
        else:
            comp = 0
        tmp = pd.DataFrame([[key, info.pos, info.chrm, info.score, num_var, dist, comp, 0]], columns=['name','pos','chrom','score','num_var','dist','correctness', 'tie_broken'])
        df_tmp = df_tmp.append(tmp)
    if var_flag:
        df_tmp = df_tmp.sort_values('score', ascending=False)
        best_score = max(df_tmp['score'])
        if df_tmp['score'].iloc[1] < best_score and df_tmp['correctness'].iloc[0] == 1:
            df_tmp['tie_broken'].iloc[0] = 1
        df_out = df_out.append(df_tmp)
    df_out.to_csv(fn_out, sep='\t', index=None)