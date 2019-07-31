import argparse
import pandas as pd
from build_erg import VarInfo

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-a', '--fn_a',
        help='Var set a'
    )
    parser.add_argument(
        '-b', '--fn_b',
        help='Var set b'
    )
    parser.add_argument(
        '-l', '--list-varset',
        help='List of paths to target .var files. \
            If this is set, dont look at -a and -b [None]'
    )
    parser.add_argument(
        '-op', '--output-prefix',
        help='Output prefix for TSVs in list mode. \
            The outputs will be <op>_inter.tsv, \
            <op>_union.tsv and <op_>_jaccard.tsv [None]'
    )
    args = parser.parse_args()
    return args

def create_set_from_file(fn):
    f = open(fn, 'r')
    l = []
    for line in f:
        v = VarInfo(line)
        encoding = str(v.chrm) + '_' + str(v.ref_pos) + '_' + str(v.ref_allele) + '_' + str(v.alt_allele)
        l.append(encoding)
    return set(l)

def calc_similarity(set_a, set_b):
    intersection = set_a.intersection(set_b)
    union = set_a.union(set_b)
    return len(intersection), len(union), len(intersection) / len(union)

def duo_mode(fn_a, fn_b):
    set_a = create_set_from_file(fn_a)
    print ('Size of set "{0:20s}" = {1}'.format(fn_a, len(set_a)))
    set_b = create_set_from_file(fn_b)
    print ('Size of set "{0:20s}" = {1}'.format(fn_b, len(set_b)))
    
    len_i, len_u, jaccard = calc_similarity(set_a, set_b)

    print ('Intersection       = {0}'.format(len_i))
    print ('Union              = {0}'.format(len_u))
    print ('Jaccard similarity = {0:.4f}'.format(jaccard))
    print ()

def list_mode(
    fn_list,
    output_prefix
):
    f_list = open(fn_list, 'r')
    
    list_set = []
    list_fn = []
    for line in f_list:
        fn = line.rstrip()
        list_set.append(create_set_from_file(fn))
        list_fn.append(fn)
        print ('Size of set "{0:20s}" = {1}'.format(fn, len(list_set[-1])))
    list_fn_base =[]
    for i in list_fn:
        base = i.split('/')[-1]
        base = base[: base.rfind('.var')]
        for s in ['EUR', 'AMR', 'AFR', 'EAS', 'SAS']:
            if base.find(s) > 0:
                base = base[base.find(s): base.find(s)+3]
                break
        list_fn_base.append(base)

    df_inter = pd.DataFrame(columns = list_fn_base, index=list_fn_base)
    df_union = pd.DataFrame(columns = list_fn_base, index=list_fn_base)
    df_jaccard = pd.DataFrame(columns = list_fn_base, index=list_fn_base)
    for i in range(df_inter.shape[0]):
        for j in range(df_inter.shape[1]):
            if i <= j:
                df_inter.iloc[i,j] = calc_similarity(list_set[i], list_set[j])[0]
                df_union.iloc[i,j] = calc_similarity(list_set[i], list_set[j])[1]
                df_jaccard.iloc[i,j] = calc_similarity(list_set[i], list_set[j])[2]

    df_inter.to_csv(output_prefix + '_inter.tsv', sep='\t')
    df_union.to_csv(output_prefix + '_union.tsv', sep='\t')
    df_jaccard.to_csv(output_prefix + '_jaccard.tsv', sep='\t')

if __name__ == '__main__':
    args = parse_args()
    fn_a = args.fn_a
    fn_b = args.fn_b
    fn_list = args.list_varset
    output_prefix = args.output_prefix

    if args.list_varset:
        list_mode(fn_list, output_prefix)
    else:
        duo_mode(fn_a, fn_b)
