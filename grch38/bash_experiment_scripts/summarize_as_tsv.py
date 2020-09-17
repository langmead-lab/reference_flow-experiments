import argparse
import os
import pandas as pd

def build_dict_indiv_to_pop(fn):
    dict_indiv_to_pop = {}
    df = pd.read_csv(fn, sep='\t')
    list_sample = df['Individual ID']
    list_pop = df['Population']
    assert len(list_sample) == len(list_pop)
    assert len(list_sample) == len(set(list_sample))
    for i, s in enumerate(list_sample):
        dict_indiv_to_pop[s] = list_pop[i]
    return dict_indiv_to_pop

def build_dict_pop_to_spop(fn):
    dict_pop_to_spop = {}
    df = pd.read_csv(fn, sep='\t', header=None)
    list_spop = df[2]
    list_pop = df[0]
    for i, p in enumerate(list_pop):
        dict_pop_to_spop[p] = list_spop[i]
    return dict_pop_to_spop

def summarize_as_tsv(fn_family, fn_spop, dir_acc, output_tsv, extension, metric):
    dict_indiv_to_pop = build_dict_indiv_to_pop(fn_family)
    dict_pop_to_spop = build_dict_pop_to_spop(fn_spop)
    
    df = pd.DataFrame()
    list_fn = os.listdir(dir_acc)
    list_fn = [os.path.join(dir_acc, fn) for fn in list_fn]
    list_exp = []
    list_indiv = []
    list_method = []
    list_tp = []
    list_all = []
    list_sensitivity = []
    for fn in list_fn:
        # if fn.endswith('.acc'):
        if fn.endswith(extension):
            with open(fn, 'r') as f:
                for i, line in enumerate(f):
                    if i == 0:
                        list_tp.append(int(line))
                    elif i == 1:
                        list_all.append(int(line))
            # fn = fn[: fn.rfind('.acc')]
            fn = fn[: fn.rfind(extension)]
            # list_exp.append(os.path.basename(fn))
            bn = os.path.basename(fn).split('-')
            indiv = bn[0]
            method = '-'.join(bn[1:])
            list_indiv.append(indiv)
            list_method.append(method)
    # df['Experiments'] = list_exp
    # df['Inidvidual'] = list_indiv
    df['Individual'] = list_indiv
    df['Super Population'] = [dict_pop_to_spop[dict_indiv_to_pop[s]] for s in list_indiv]
    df['Method'] = list_method
    # df['True Positive'] = list_tp
    df[metric] = list_tp
    df['Num Reads'] = list_all
    df.to_csv(output_tsv, sep='\t', index=None, float_format = '%.4f')


def summarize_lowq(fn_family, fn_spop, input_tsv, output_tsv):
    dict_indiv_to_pop = build_dict_indiv_to_pop(fn_family)
    dict_pop_to_spop = build_dict_pop_to_spop(fn_spop)
    
    df = pd.read_csv(input_tsv, sep=',')
    list_indiv = list(df.iloc[:,0])
    df['Super Population'] = [dict_pop_to_spop[dict_indiv_to_pop[s]] for s in list_indiv]
    df.to_csv(output_tsv, index=None, sep='\t')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fn_family',
        help='1KGP family file')
    parser.add_argument(
        '--fn_spop',
        help='1KGP superpopulation file')
    parser.add_argument(
        '--dir_acc',
        help='Directory storing all .acc files. For "acc" and "unaligned" modes.')
    parser.add_argument(
        '--output_tsv',
        help='Output TSV filename')
    parser.add_argument(
        '--input_tsv',
        help='Input TSV filename. For "unaligned" mode.')
    parser.add_argument(
        '--mode', default="acc",
        help='Operating mode. "acc", "lowq", "unaligned", "num_incorrect" ["acc"]')
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_args()
    if args.mode == 'acc':
        summarize_as_tsv(args.fn_family, args.fn_spop, args.dir_acc, args.output_tsv, '.acc', 'True Positive')
    elif args.mode == 'unaligned':
        summarize_as_tsv(args.fn_family, args.fn_spop, args.dir_acc, args.output_tsv, '.unaligned', 'Unaligned')
    elif args.mode == 'num_incorrect':
        summarize_as_tsv(args.fn_family, args.fn_spop, args.dir_acc, args.output_tsv, '.num_incorrect', 'Num Incorrect')
    elif args.mode == 'lowq':
        summarize_lowq(args.fn_family, args.fn_spop, args.input_tsv, args.output_tsv)
    else:
        print ('unsupported mode', args.mode)
        exit()
