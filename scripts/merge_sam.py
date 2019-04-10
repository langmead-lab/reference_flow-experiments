'''
Merge two SAM files based on xxx (probably AS)
'''
import argparse
import sys
from analyze_sam import SamInfo, parse_line

def merge_sam(args):
    sam1_fn = args.sam1
    sam2_fn = args.sam2
    id1 = args.id1
    id2 = args.id2
    sam1_f = open(sam1_fn, 'r')
    sam2_f = open(sam2_fn, 'r')
    
    #: prepare output files
    sam1_prefix = sam1_fn[: sam1_fn.find('.sam')]
    sam2_prefix = sam2_fn[: sam2_fn.find('.sam')]
    sam1_out_fn = sam1_prefix.split('/')[-1] + '-merged_' + id2 + '.sam'
    sam2_out_fn = sam2_prefix.split('/')[-1] + '-merged_' + id1 + '.sam'
    sam1_out_f = open(sam1_out_fn, 'w')
    sam2_out_f = open(sam2_out_fn, 'w')

    sam2_dic = {}
    num_replacement = 0
    for line in sam2_f:
        name, info = parse_line(line, erg=True)
        if name == 'header':
            sam2_out_f.write(line)
            continue
        if info.is_unaligned():
            continue
        sam2_dic[name] = [info, line]
    for line in sam1_f:
        name, info = parse_line(line, erg=True)
        if name == 'header':
            sam1_out_f.write(line)
            continue
        if sam2_dic.get(name) == None:
            sam1_out_f.write(line)
        else:
            s2_info = sam2_dic[name]
            s2_score = s2_info[0].score
            if s2_score > info.score:
                num_replacement += 1
                sam2_out_f.write(s2_info[1])
            else:
                sam1_out_f.write(line)
    sys.stderr.write ('num replacement = %d\n' % num_replacement)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n1', '--sam1',
        help='sam file 1'
    )
    parser.add_argument(
        '-id1', '--id1',
        help='id for sam file 1'
    )
    parser.add_argument(
        '-n2', '--sam2',
        help='sam file 2, expect this to be aligned to a partial genome'
    )
    parser.add_argument(
        '-id2', '--id2',
        help='id for sam file 2'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    merge_sam(args)