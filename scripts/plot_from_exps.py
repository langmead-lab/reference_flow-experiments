import argparse
import pandas as pd
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s', '--samples',
        help='txt file recording the order of samples used in experiments'
    )
    parser.add_argument(
        '-l', '--log',
        help='log file for experiments'
    )
    parser.add_argument(
        '-p', '--ped',
        help='ped file recoding population and family info'
    )
    return parser.parse_args()

def read_true_pos(log_fn):
    list_tp = []
    with open(log_fn, 'r') as f:
        for line in f:
            if line.startswith('sensitivity_all'):
                tp = int(line.split()[3][1:])
                list_tp.append(tp)
    return list_tp

def read_samples(sample_fn):
    list_samples = []
    with open(sample_fn, 'r') as f:
        for line in f:
            list_samples.append(line.rstrip())
    return list_samples

def read_ped(ped_fn):
    '''
    Reads the ped file and returns a dict where keys are indivs and 
    values are corresponding populations
    '''
    ped_df = pd.read_csv(ped_fn, sep='\t')
    dd = ped_df[['Individual ID', 'Population']]

    popd = {}
    for i in range(dd.shape[0]):
        key = dd['Individual ID'][i]
        value = dd['Population'][i]
        popd[key] = value

    return popd

if __name__ == '__main__':
    args = parse_args()
    sample_fn = args.samples
    log_fn = args.log
    ped_fn = args.ped
    
    TP_OFFSET = 19719900 # num of TP in the first pass
    TP_FLOOR  = 19908912 # num of TP using purely hg19
    TP_CEIL   = 19929495 # num of TP using personalized genomes

    list_tp = read_true_pos(log_fn)
    list_samples = read_samples(sample_fn)
    data_df = pd.DataFrame(list_tp)
    data_df.columns = ['Num TP']
    data_df.index = list_samples
    dict_pop = read_ped(ped_fn)
    list_pop = [dict_pop[data_df.index[i]] for i in range(data_df.shape[0])]
    data_df['Population'] = list_pop
    
    data_df['Num TP'] += TP_OFFSET

    list_tp = [TP_FLOOR] + list(data_df['Num TP']) + [TP_CEIL]
    list_sample = ['hg19'] + list(data_df.index) + ['personalized']
