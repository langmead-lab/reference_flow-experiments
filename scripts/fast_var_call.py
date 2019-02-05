'''
Calls variants without correction
'''
import argparse
import re
from analyze_sam import SamInfo, parse_line

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n', '--sam',
        help='target sam file'
    )
    # parser.add_argument(
    #     '-g', '--golden',
    #     help='golden sam file'
    # )
    # parser.add_argument(
    #     '-v', '--var',
    #     help='the file specifying variants'
    # )
    # parser.add_argument(
    #     '-t', '--threshold', type=int,
    #     default=10,
    #     help='(int) max allowed distance for a correct mapping [10]'
    # )
    # parser.add_argument(
    #     '--read_len', type=int,
    #     default=100,
    #     help='(int) read length [100]'
    # )
    # parser.add_argument(
    #     '-p', '--personalized', type=int,
    #     default=0,
    #     help='(int) specify whether the ref seq(s) are standard (0) or personalized-diploid (2) sample [0]'
    # )
    # parser.add_argument(
    #     '--step_size', type=int,
    #     default=1000,
    #     help='(int) the step size for main/alt offset indexes [1000]'
    # )
    # parser.add_argument(
    #     '--write_wrt_correctness', type=int,
    #     default=None,
    #     help='(int) If set, writes two files recording correct/incorrect alignments respectively. The output files use target sam prefix [None].'
    # )
    # parser.add_argument(
    #     '--write_wrt_mapq', type=int,
    #     default=None,
    #     help='(int) If specified, writes two files recording alignments with mapq >= t and mapq < t. This argument is the threshold. The output files use target sam prefix [None].'
    # )
    args = parser.parse_args()
    return args

def parse_cigar_and_md(info):
    tmp_var = []
    tmp_offset = 0
    # print (info.cigar)
    cigar = re.findall(r'\d+[M,I,D]', info.cigar)
    tag_md = re.findall(r'\d+|[\^,A-Z]+', info.tag_md)
    # print (re.findall(r'[\^,A-Z]+', info.tag_md))
    md_length = 0
    for md_element in tag_md:
        if md_element.isdigit():
            md_length += int(md_element)
        elif md_element.startswith('^') == False:
            tmp_var.append([md_length, md_element, ''])
            md_length += len(md_element)
        else:
            # deletion
            seq = md_element[1:]
            tmp_var.append([md_length, seq, '-'])
    cigar_length = 0
    ins_length = 0
    for cigar_element in cigar:
        if cigar_element[-1] == 'I':
            ins_idx = -1
            for i, v in enumerate(tmp_var):
                if ins_idx != -1:
                    print ('before', v[0])
                    v[0] += int(cigar_element[: -1])
                    print ('after', v[0])
                if v[0] > cigar_length and ins_idx == -1:
                    print ('before', v[0])
                    ins_idx = i
                    v[0] += int(cigar_element[: -1])
                    print ('after', v[0])
            tmp_var.append([cigar_length, '-', 'INS'])
            cigar_length += int(cigar_element[: -1])
            ins_length += int(cigar_element[: -1])
        elif cigar_element[-1] == 'M':
            cigar_length += int(cigar_element[: -1])
    
    assert md_length + ins_length == len(info.seq)
    if info.is_rc() == False:
        print ('cigar =', cigar)
        print ('MD:Z =', info.tag_md)
        print (tag_md)
        info.print(chrm=False, offset=False, mapq=False)
        print (tmp_var)
        input ()
    
    print ()
    
    
    '''
    How to represent a deletion and a SNV in a MD string?
        (Guess) adds one base to the deletion and adds an insertion (hidden in MD)
    '''

def fast_var_call(argrs):
    sam_fn = args.sam
    f = open(sam_fn, 'r')
    for line in f:
        name, info = parse_line(line, by_score=0, md=True, cigar=True)
        if name == 'header':
            continue
        if info.is_unaligned():
            continue
        # skipped perfect alignments
        try:
            # 100 matches with no SNV
            if len(info.seq) == int(info.tag_md) and info.cigar == str(len(info.seq)) + 'M':
                # print ('perfect', info.tag_md)
                continue
        except:
            parse_cigar_and_md(info)    
    
    print ('finished')

if __name__ == '__main__':
    args = parse_args()
    fast_var_call(args)
