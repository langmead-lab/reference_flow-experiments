'''
Merge a set of SAM files.

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
import os
from analyze_sam import SamInfo, parse_line

def compare_score_and_mapq(list_info):
    order = list(range(len(list_info)))
    # list_is_unaligned = [0 if info.is_unaligned() else 1 for info in list_info]
    # list_as = [info.score for info in list_info]
    # list_mapq = [info.mapq for info in list_info]
    list_is_unaligned = [info[0] != 1 for info in list_info]
    list_as = [info[0] for info in list_info]
    list_mapq = [info[1] for info in list_info]

    #: sort order: if_aligned > AS > MAPQ
    #: so perform sorting in reversed order
    list_mapq, list_as, list_is_unaligned, order = \
        zip(
            *sorted(zip(list_mapq, list_as, list_is_unaligned, order), reverse = True)
        )
    list_as, list_mapq, list_is_unaligned, order = \
        zip(
            *sorted(zip(list_as, list_mapq, list_is_unaligned, order), reverse = True)
        )
    list_is_unaligned, list_as, list_mapq, order = \
        zip(
            *sorted(zip(list_is_unaligned, list_as, list_mapq, order), reverse = True)
        )

    for i in range(1, len(list_info)):
        #: if not a tie
        if (list_is_unaligned[i] != list_is_unaligned[0]) or \
            (list_as[i] != list_as[0]) or \
            (list_mapq[i] != list_mapq[0]):
            order = order[:i]
            break
    return random.sample(order, 1)

def get_name_as_mapq(line):
    if line[0] == '@':
        return 'header', None
    line = line.split()
    name = line[0]
    mapq = int(line[4])
    score = 1 #: represents unmapped
    for i in line:
        if i.startswith('AS:i'):
            score = int(i.split(':')[-1])
    return name, (score, mapq)

def merge_sam_list(args):
    fn_sam = args.sam_list
    fn_ids = args.id_list
    tie_break = args.tie_break
    assert tie_break in ['1', '2', 'r']
    random.seed(args.rand_seed)
    if args.rand_seed:
        print ('Set random seed: {}'.format(args.rand_seed))
    fn_log = args.log
    if fn_log != None:
        f_log = open(fn_log, 'w')

    with open(fn_sam, 'r') as f:
        list_fn_sam = [] #: list of SAM names
        list_f_sam = [] #: list of opened SAM files
        for line in f:
            fn_sam = line.rstrip()
            list_fn_sam.append(fn_sam)
            list_f_sam.append(open(fn_sam, 'r'))
    with open(fn_ids, 'r') as f:
        list_ids = []
        for line in f:
            list_ids.append(line.rstrip())

    assert len(list_fn_sam) == len(list_ids)
    num_sam = len(list_fn_sam)

    #: prepare output files
    list_f_out = []
    for i in range(num_sam):
        fn_sam = list_fn_sam[i]
        sam_prefix = fn_sam[: fn_sam.find('.sam')]
        #fn_out = sam_prefix.split('/')[-1] + '+' + ('_').join(list_ids[0:i]+list_ids[i+1:]) + '.sam'
        # fn_out = sam_prefix + '+' + ('_').join(list_ids[0:i]+list_ids[i+1:]) + '.sam'
        fn_out = fn_sam + '+' + ('_').join(list_ids[0:i]+list_ids[i+1:]) + '.sam'
        if fn_log != None:
            f_log.write(fn_out + '\n')
        print (fn_out)
        list_f_out.append(open(fn_out, 'w'))
    print ('------')
    
    #: save aligned reads into a list of dictionaries
    list_dicts = []
    for i in range(num_sam):
        list_dicts.append({})
        sys.stderr.write('Loading ' + str(i) + '\n')
        for line in list_f_sam[i]:
            # name, info = parse_line(line, erg=True)
            name, info = get_name_as_mapq(line)
            if name == 'header':
                list_f_out[i].write(line)
                continue
            list_dicts[i][name] = info
            # list_dicts[i][name] = [info, line]
        sys.stderr.write('Completed ' + str(i) + '\n')

    #: check if lists have the same key (i.e. read names)
    list_set_dict = []
    list_read_names = []
    for i in range(num_sam):
        list_set_dict.append(set(list_dicts[i].keys()))
        list_read_names.append([])
        if i > 0:
            assert list_set_dict[0] == list_set_dict[i]
    
    #: close files
    for i in range(num_sam):
        list_f_sam[i].close()

    for read_name in list_dicts[0].keys():
        curr = [list_dicts[i][read_name] for i in range(num_sam)]
        selected_idx = compare_score_and_mapq(curr)
        assert len(selected_idx) == 1
        selected_idx = selected_idx[0]
        list_read_names[selected_idx].append(read_name)
        # list_f_out[selected_idx].write(list_dicts[selected_idx][read_name][1])

    for i in range(num_sam):
        list_read_names[i] = set(list_read_names[i])
        sys.stderr.write(str(list_fn_sam[i]) + '\n')
        sys.stderr.write(str(len(list_read_names[i])) + '\n')
        with open(list_fn_sam[i], 'r') as f_sam:
            for line in f_sam:
                name, info = parse_line(line)
                if name in list_read_names[i]:
                    list_f_out[i].write(line)

    #: close files
    for i in range(num_sam):
        # list_f_sam[i].close()
        list_f_out[i].close()

def merge_sam(args):
    sam1_fn = args.sam1
    sam2_fn = args.sam2
    id1 = args.id1
    id2 = args.id2
    tie_break = args.tie_break
    assert tie_break in ['1', '2', 'r']
    random.seed(args.rand_seed)
    if args.rand_seed:
        print ('Set random seed: {}'.format(args.rand_seed))
    fn_log = args.log

    sam1_f = open(sam1_fn, 'r')
    sam2_f = open(sam2_fn, 'r')
    
    #: prepare output files
    sam1_prefix = sam1_fn[: sam1_fn.find('.sam')]
    sam2_prefix = sam2_fn[: sam2_fn.find('.sam')]
    # sam1_out_fn = sam1_prefix.split('/')[-1] + '-merged_' + id2 + '.sam'
    # sam2_out_fn = sam2_prefix.split('/')[-1] + '-merged_' + id1 + '.sam'

    # sam1_out_fn = sam1_prefix + '-merged_' + id2 + '.sam'
    # sam2_out_fn = sam2_prefix + '-merged_' + id1 + '.sam'

    sam1_out_fn = sam1_fn + '-merged_' + id2 + '.sam'
    sam2_out_fn = sam2_fn + '-merged_' + id1 + '.sam'
    sam1_out_f = open(sam1_out_fn, 'w')
    sam2_out_f = open(sam2_out_fn, 'w')

    if fn_log != None:
        with open(fn_log, 'w') as f:
            f.write(os.path.abspath(sam1_out_fn) + '\n')
            f.write(os.path.abspath(sam2_out_fn))

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
    parser.add_argument(
        '-ns', '--sam-list',
        help='list of paths of SAM files'
    )
    parser.add_argument(
        '-ids', '--id-list',
        help='list of ids of files'
    )
    parser.add_argument(
        '-rs', '--rand-seed',
        help='random seed for controlled randomness [None]'
    )
    parser.add_argument(
        '-l', '--log',
        help='log file recording paths of merged files'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    #if args.list == 1:
    if args.sam_list and args.id_list:
        print ('Note: list mode enabled, outputs:')
        merge_sam_list(args)
    else:
        merge_sam(args)
