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
        help='0: paired-end mode, 1: unpaired_1, 2: unpaired_2 [0]'
    )
    args = parser.parse_args()
    return args

def parse_line(line):
    if line[0] == '@':
        return False, False
    line = line.split()
    name = line[0]
    flag = '{0:012b}'.format(int(line[1]))
    chrm = line[2]
    pos = int(line[3])
    return name, [chrm, pos, flag]

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
    num_reads = 0
    num_correct_reads = 0
    num_incorrect_reads = 0
    g_dic = {}
    with open(golden, 'r') as gfile:
        for line in gfile:
            name, info = parse_line(line)
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
                if flag[5] is '1': # first_seg
#                    print (int(flag, 2))
#                    print ("[5] flag = %s" % flag)
#                    print (flag[5:])
                    g_dic[name] = info
                    num_reads += 1
        print (len(g_dic))

    with open(sam, 'r') as infile:
        for line in infile:
            name, info = parse_line(line)
            if name is False:
                continue
            if info[2][3] == '1':
                # neglect secondary alignments
                continue
            comp = compare_sam_info(info, g_dic[name], threshold)
            if comp:
                num_correct_reads += 1
            else:
                num_incorrect_reads += 1
                #print (info)
                #print (g_dic[name])
                #input()
    
    if num_reads != num_correct_reads + num_incorrect_reads:
        print ('Error: number of reads do not match!')
        print ('Total reads = %d, correct reads = %d, incorrect reads = %d' % \
                (num_reads, num_correct_reads, num_incorrect_reads))
        return 
    print ('Alignment accuracy = %f (%d/%d)' % \
        (100 * float(num_correct_reads) / num_reads, num_correct_reads, num_reads))

if __name__ == '__main__':
    args = parse_args()
    analyze_sam(args)
