'''
Merge two SAM files based on xxx (probably AS)
'''
import argparse
import sys
from analyze_sam import SamInfo, parse_line

def merge_sam(args):
    sam1_fn = args.sam1
    sam2_fn = args.sam2
    sam1_f = open(sam1_fn, 'r')
    sam2_f = open(sam2_fn, 'r')

    sam2_dic = {}
    num_replacement = 0
    for line in sam2_f:
        name, info = parse_line(line, erg=True)
        if name == 'header':
            continue
        if info.is_unaligned():
            continue
        sam2_dic[name] = [info, line]
    for line in sam1_f:
        name, info = parse_line(line, erg=True)
        if name == 'header':
            print (line[:-1])
            continue
        if sam2_dic.get(name) == None:
            print (line[:-1])
        else:
            s2_info = sam2_dic[name]
            s2_score = s2_info[0].score
            if s2_score > info.score:
                num_replacement += 1
                print (s2_info[1][:-1])
            else:
                print (line[:-1])
    sys.stderr.write ('num replacement = %d\n' % num_replacement)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n1', '--sam1',
        help='sam file 1'
    )
    parser.add_argument(
        '-n2', '--sam2',
        help='sam file 2, expect this to be aligned to a partial genome'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    merge_sam(args)