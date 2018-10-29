'''
Compare aligned sam with synthetic golden data
'''
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n', '--sam',
        help='aligned sam file'
    )
    parser.add_argument(
        '-g', '--golden',
        help='golden sam file'
    )
    parser.add_argument(
        '-t', '--threshold', type=int,
        default=10,
        help='max allowed distance for a correct mapping [10]'
    )
    parser.add_argument(
        '-s', '--seg', type=int,
        default=1,
        help='0: paired-end mode, 1: unpaired_1, 2: unpaired_2 [1]'
    )
    parser.add_argument(
        '--strat_by_as', type=int,
        default=0,
        help='0: normal mode; 1: evaluate with stratification by alignment score [0]'
    )
    args = parser.parse_args()
    return args

def parse_line(line, strat_by_as):
    if line[0] == '@':
        return False, False, 0
    line = line.split()
    name = line[0]
    flag = '{0:012b}'.format(int(line[1]))
    chrm = line[2]
    pos = int(line[3])
    if strat_by_as > 0:
        if flag[9] is '1': # unaligned
            return name, [chrm, pos, flag], 1
        score = line[11] # AS:i:xx, bt2 format
        if score.startswith('AS:') is False:
            print ('Error: incorrect AS information!')
            print (line)
            return
        score = int(score.split(':')[-1])
        return name, [chrm, pos, flag], score
    return name, [chrm, pos, flag], 0

def compare_sam_info(info, ginfo, threshold):
    if __debug__:
        print (info)
        print (ginfo)
    if info[0] != ginfo[0]:
        # diff chromosome
        if __debug__: print ("False: chr")
        return False
    if info[2][7] != ginfo[2][7]:
        # diff direction
        if __debug__: print ("False: direction (%s, %s)" % (info[2][7], ginfo[2][7]))
        return False
    if abs(info[1] - ginfo[1]) > threshold:
        if __debug__: print ("False: 2")
        return False
    if __debug__: print ("True")
    return True

def analyze_sam(args):
    sam = args.sam
    golden = args.golden
    threshold = args.threshold
    seg = args.seg
    strat_by_as = args.strat_by_as
    num_reads = 0
    num_correct_reads = 0
    num_incorrect_reads = 0
    g_dic = {}
    with open(golden, 'r') as gfile:
        for line in gfile:
            name, info, _ = parse_line(line, 0)
            if name is False:
                continue
            flag = info[2]
            if seg == 1:
                if name in g_dic and flag[5] is '1':
                    print ("Error: duplicated reads in golden!")
                    print (name)
                    print (info)
                    print (g_dic[name])
                    return
                if flag[5] is '1': # is first_seg
                    g_dic[name] = info
                    num_reads += 1
        print ('Size of database:', len(g_dic))

    with open(sam, 'r') as infile:
        if strat_by_as > 0:
            cor_strat = {}
            incor_strat = {}
        for line in infile:
            name, info, score = parse_line(line, strat_by_as)
            if name is False:
                continue
            if info[2][3] == '1':
                # neglect secondary alignments
                continue
            comp = compare_sam_info(info, g_dic[name], threshold)
            if __debug__:
                input()
            if comp:
                num_correct_reads += 1
            else:
                num_incorrect_reads += 1
            if strat_by_as > 0:
                if cor_strat.get(score) is None:
                    cor_strat[score] = 0
                if incor_strat.get(score) is None:
                    incor_strat[score] = 0
                if comp:
                    cor_strat[score] += 1
                else:
                    incor_strat[score] += 1
    
    if strat_by_as > 0:
        print (cor_strat)
        print (incor_strat)

    if num_reads != num_correct_reads + num_incorrect_reads:
        print ('Warning: number of reads do not match!')
        print ('Total reads = %d, correct reads = %d, incorrect reads = %d' % \
                (num_reads, num_correct_reads, num_incorrect_reads))
        # return 
    print ('Alignment accuracy = %f (%d/%d)' % \
        (100 * float(num_correct_reads) / num_reads, num_correct_reads, num_reads))

if __name__ == '__main__':
    args = parse_args()
    analyze_sam(args)
