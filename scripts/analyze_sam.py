'''
Last update: 2018/11/28 by Nae-Chyun Chen

Analyzes aligned sam with synthetic golden data
and measures sensitivity
'''
import argparse
import pickle
import os

class SamInfo:
    '''
    Records information of a sam line
    '''
    pos = 0
    chrm = ''
    flag = 0
    mapq = 0
    score = 0

    def __init__(self, pos, chrm, flag, mapq, score):
        self.pos = int(pos)
        self.chrm = chrm
        self.flag = int(flag)
        self.mapq = int(mapq)
        self.score = int(score)
    
    def print(self):
        print ('pos =', self.pos)
        print ('chrm =', self.chrm)
        print ('flag =', self.flag)
        print ('mapq =', self.mapq)
        print ('score =', self.score)

    def is_unaligned(self):
        if self.flag & 4: return True
        return False

    def is_rc(self):
        if self.flag & 16: return True
        return False

    def is_first_seq(self):
        if self.flag & 64: return True
        return False

    def is_secondary(self):
        if self.flag & 256: return True
        return False

    def update_score(self, raw_score):
        if raw_score.startswith('AS:') is False:
            print ('Error: incorrect AS information!')
            print (raw_score)
            return False
        self.score = int(raw_score.split(':')[-1])
        return True


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
        '--haploid', type=int,
        default=1,
        help='1: haploid input, 2: diploid input [1]'
    )
    parser.add_argument(
        '--by_score', type=int,
        default=0,
        help='0: normal mode; 1: evaluate by alignment score [0]'
    )
    parser.add_argument(
        '--by_mapq', type=int,
        default=0,
        help='0: normal mode; 1: evaluate by mapping quality [0]'
    )
    parser.add_argument(
        '--secondary', type=int,
        default=0,
        help='0: do not consider secondary alignments; 1: consider secondary alignment [0]'
    )
    args = parser.parse_args()
    return args

def parse_line(line, by_score):
    if line[0] == '@':
        return False, False
    line = line.split()
    name = line[0]
    flag = int(line[1])
    chrm = line[2]
    pos = int(line[3])
    mapq = int(line[4])
    info = SamInfo(pos, chrm, flag, mapq, 0)

    if by_score > 0:
        if info.is_unaligned():
            info.score = 1
            return name, info
        if info.update_score(line[11]):
            return name, info
        else:
            return False
    return name, info

def compare_sam_info(info, ginfo, threshold):
    if info.chrm != ginfo.chrm:
        # diff chromosome
        if __debug__: print ("False: chr, mapq =", info.mapq)
        return False
    if (info.is_rc() ^ ginfo.is_rc()) is True:
        # diff direction
        if __debug__: 
            print ("False: direction (%s, %s)" % (info.is_rc(), ginfo.is_rc()), \
                    "mapq =", info.mapq)
        return False
    if abs(info.pos - ginfo.pos) > threshold:
        if __debug__: print ("False: distance > threshold, mapq =", info.mapq)
        return False
    if __debug__: print ("True, mapq =", info.mapq)
    return True

def dump_golden_dic(filename, seg):
    g_dic = {}
    with open(filename, 'r') as gfile:
        for line in gfile:
            name, info = parse_line(line, 0)
            if name is False:
                continue
            if seg == 1:
                if name in g_dic and info.is_first_seq():
                    print ("Error: duplicated reads in golden!")
                    print (name)
                    info.print()
                    g_dic[name].print()
                    return
                if info.is_first_seq():
                    g_dic[name] = info
        print ("Size of database built:", len(g_dic))
    
    pkl_filename = filename + '.pkl'
    print ("Dump %s using pickle" % pkl_filename)
    f = open(pkl_filename, 'wb')
    pickle.dump(g_dic, f)
    f.close()
    return g_dic

def load_golden_dic(filename, seg):
    pkl_filename = filename + '.pkl'
    if os.path.isfile(pkl_filename):
        f = open(pkl_filename, 'rb')
        g_dic = pickle.load(f)
        f.close()
        print ('Size of database loaded:', len(g_dic))
        return g_dic
    else:
        return dump_golden_dic(filename, seg)

def analyze_sam(args):
    sam = args.sam
    golden = args.golden
    threshold = args.threshold
    seg = args.seg
    haploid = args.haploid
    by_score = args.by_score
    by_mapq = args.by_mapq
    secondary = args.secondary

    num_correct_reads = 0
    num_incorrect_reads = 0
    dic_secondary = {}
    # build golden dictionary
    # g_dic = dump_golden_dic(golden, seg)
    # load pre-built golden dictionary
    g_dic = load_golden_dic(golden, seg)
    num_syn_reads = len(g_dic)
    with open(sam, 'r') as infile:
        if by_score > 0:
            dic_as_cor = {}
            dic_as_incor = {}
        if by_mapq > 0:
            dic_q_cor = {}
            dic_q_incor = {}
        for line in infile:
            name, info = parse_line(line, by_score)
            if name is False:
                continue
            if info.is_secondary() and secondary == 0: # neglect secondary alignments
                continue
            comp = compare_sam_info(info, g_dic[name], threshold)
            if __debug__ & 0:
                if comp is not True:
                    print (info)
                    print (g_dic[name])
                    input()
            if comp:
                num_correct_reads += 1
            else:
                num_incorrect_reads += 1
            if by_score > 0:
                if dic_as_cor.get(info.score) is None:
                    dic_as_cor[info.score] = 0
                if dic_as_incor.get(info.score) is None:
                    dic_as_incor[info.score] = 0
                if comp:
                    dic_as_cor[info.score] += 1
                else:
                    dic_as_incor[info.score] += 1
            if by_mapq > 0:
                if dic_q_cor.get(info.mapq) is None:
                    dic_q_cor[info.mapq] = 0
                if dic_q_incor.get(info.mapq) is None:
                    dic_q_incor[info.mapq] = 0
                if comp:
                    dic_q_cor[info.mapq] += 1
                else:
                    dic_q_incor[info.mapq] += 1
            # TODO
            if secondary > 0:
                #if dic_secondary.get(name) is None:
                #    dic_secondary[name] = 0
                if comp:
                    if dic_secondary.get(name) != 1:
                        dic_secondary[name] = 1
                        if __debug__:
                            info.print()
                            g_dic[name].print()
                            print ("SEC: precision +=1")
                            input()
                    # elif dic_secondary.get(name) == 1:
                    #     print (info)
                    #     print (g_dic[name])
    
    num_reads = num_correct_reads + num_incorrect_reads

    if num_syn_reads != num_correct_reads + num_incorrect_reads:
        print ('Warning: number of reads do not match!')
        print ('Synthesized reads = %d, correct reads = %d, incorrect reads = %d' % \
                (num_syn_reads, num_correct_reads, num_incorrect_reads))
        # return 
    print ('Alignment sensitivity = %f (%d/%d)' % \
        (100 * float(num_correct_reads) / num_reads, num_correct_reads, num_reads))
    
    if by_score > 0:
        print (dic_as_cor)
        print (dic_as_incor)
    if by_mapq > 0:
        print (dic_q_cor)
        print (dic_q_incor)
    if secondary > 0:
        print (len(dic_secondary))

if __name__ == '__main__':
    args = parse_args()
    analyze_sam(args)
