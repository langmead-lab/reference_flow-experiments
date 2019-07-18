'''
Merge two SAM files.

The selection rule is Alignment Score (AS) > Mapping Quality (MAPQ)

If a read is aligned with same AS and same MAPQ in the provided SAM files,
it is assigned based on user-provided tie-breaking rule:
    (1) always use first SAM record
    (2) always use second SAM record
    (3) randomly assign
'''
import argparse
import sys
import random
from analyze_sam import SamInfo, parse_line

def merge_sam(args):
    sam1_fn = args.sam1
    sam2_fn = args.sam2
    id1 = args.id1
    id2 = args.id2
    tie_break = args.tie_break
    assert tie_break in ['1', '2', 'r']

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

    ## WIP
    num_sameAS_higherQ = 0
    num_sameAS_sameAS = 0

    #: save aligned reads into a dictionary
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
            s2_mapq = s2_info[0].mapq

            if s2_score > info.score:
                num_replacement += 1
                sam2_out_f.write(s2_info[1])
            #: Same AS, MAPQ_2 is higher
            elif s2_score == info.score and s2_mapq > info.mapq:
                num_replacement += 1
                num_sameAS_higherQ += 1
                sam2_out_f.write(s2_info[1])
            #: Same AS, same MAPQ
            elif s2_score == info.score and s2_mapq == info.mapq:
                num_sameAS_sameAS += 1
                if tie_break == 'r':
                    if random.random() > 0.5:
                        num_replacement += 1
                        sam2_out_f.write(s2_info[1])
                    else:
                        sam1_out_f.write(line)
                elif tie_break == '1':
                    sam1_out_f.write(line)
                else:
                    num_replacement += 1
                    sam2_out_f.write(s2_info[1])
            else:
                sam1_out_f.write(line)
    sys.stderr.write ('num replacement = %d\n' % num_replacement)
    print ('higherQ', num_sameAS_higherQ)
    print ('sameASsameQ', num_sameAS_sameAS)

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
    parser.add_argument(
        '-tb', '--tie-break',
        default='r',
        help='Ties happen when two alignments have same AS and MAPQ. \
            Tie breaking rule: \
                "1": always use -n1 \
                "2": always use -n2 \
                "r": random selection [r]'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    merge_sam(args)
