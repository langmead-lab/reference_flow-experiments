import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '-s', '--stats',
    help='target stats file (-stats.pkl)'
)
args = parser.parse_args()
fn_stats = args.stats
df = pd.read_pickle(fn_stats)

dict_correct = {}
dict_incorrect = {}

for i, name in enumerate(df['name']):
    #: ignores if any alignment for a read is correct
    if name in dict_correct:
        continue
    else:
        #: finds a correct alignment
        if df['dist'][i] >= 0 and df['dist'][i] <= 10:
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