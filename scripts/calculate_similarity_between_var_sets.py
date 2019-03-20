import argparse
from build_erg import VarInfo

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-a', '--varset_a_fn',
        help='var set a'
    )
    parser.add_argument(
        '-b', '--varset_b_fn',
        help='var set b'
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
    print ('intersection       = {0}'.format(len(intersection)))
    print ('union              = {0}'.format(len(union)))
    print ('Jaccard similarity = {0:.4f}'.format(len(intersection) / len(union)))

if __name__ == '__main__':
    args = parse_args()
    varset_a_fn = args.varset_a_fn
    varset_b_fn = args.varset_b_fn

    set_a = create_set_from_file(varset_a_fn)
    print ('size of set "{0:20s}" = {1}'.format(varset_a_fn, len(set_a)))
    set_b = create_set_from_file(varset_b_fn)
    print ('size of set "{0:20s}" = {1}'.format(varset_b_fn, len(set_b)))

    calc_similarity(set_a, set_b)
    print ()