'''
Compare aligned sam with synthetic golden data
'''
import argparse
import pickle

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
        '-ga', '--golden_hapa',
        help='golden hapa sam file'
    )
    parser.add_argument(
        '-gb', '--golden_hapb',
        help='golden hapb sam file'
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
        return False, False, 0
    line = line.split()
    name = line[0]
    flag = '{0:012b}'.format(int(line[1]))
    chrm = line[2]
    pos = int(line[3])
    mapq = int(line[4])
    if by_score > 0:
        if flag[9] is '1': # unaligned
            return name, [chrm, pos, flag], 1
        score = line[11] # AS:i:xx, bt2 format
        if score.startswith('AS:') is False:
            print ('Error: incorrect AS information!')
            print (line)
            return
        score = int(score.split(':')[-1])
        return name, [chrm, pos, flag], score
    return name, [chrm, pos, flag, mapq], 0

def compare_sam_info(info, ginfo, threshold):
    if info[0] != ginfo[0]:
        # diff chromosome
        if __debug__: print ("False: chr, mapq =", info[3])
        return False
    if info[2][7] != ginfo[2][7]:
        # diff direction
        if __debug__: 
            print ("False: direction (%s, %s)" % (info[2][7], ginfo[2][7]), \
                    "mapq =", info[3])
        return False
    if abs(info[1] - ginfo[1]) > threshold:
        if __debug__: print ("False: distance > threshold, mapq =", info[3])
        return False
    if __debug__: print ("True, mapq =", info[3])
    return True

def dump_golden_dic(filename, seg):
    g_dic = {}
    with open(filename, 'r') as gfile:
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
        print ("Size of database built:", len(g_dic))
    
    pkl_filename = filename + '.pkl'
    print ("Dump %s using pickle" % pkl_filename)
    f = open(pkl_filename, 'wb')
    pickle.dump(g_dic, f)
    f.close()
    return g_dic

def load_golden_dic(pkl_filename):
    f = open(pkl_filename, 'rb')
    g_dic = pickle.load(f)
    f.close()
    print ('Size of database loaded:', len(g_dic))
    return g_dic

def analyze_sam(args):
    sam = args.sam
    golden = args.golden
    ga = args.golden_hapa
    gb = args.golden_hapb
    threshold = args.threshold
    seg = args.seg
    by_score = args.by_score
    by_mapq = args.by_mapq
    secondary = args.secondary

    num_correct_reads_a = 0
    num_incorrect_reads_a = 0
    num_correct_reads_b = 0
    num_incorrect_reads_b = 0
    dic_secondary = {}
    # build golden dictionary
    #g_dic = dump_golden_dic(golden, seg)
    #ga_dic = dump_golden_dic(ga, seg)
    #gb_dic = dump_golden_dic(gb, seg)
    # load pre-built golden dictionary
    #g_dic = load_golden_dic(golden + '.pkl')
    ga_dic = load_golden_dic(ga + '.pkl')
    gb_dic = load_golden_dic(gb + '.pkl')
    num_syn_reads = len(ga_dic)+len(gb_dic)
    #if len(gb_dic) != num_syn_reads:
    #    print ("Error: Size of dics don't match (%s, %s)" % (len(ga_dic), len(gb_dic)))
    #    return

    with open(sam, 'r') as infile:
        if by_score > 0:
            dic_as_cor_a = {}
            dic_as_incor_a = {}
            dic_as_cor_b = {}
            dic_as_incor_b = {}
        if by_mapq > 0:
            dic_q_cor_a = {}
            dic_q_incor_a = {}
            dic_q_cor_b = {}
            dic_q_incor_b = {}
        for line in infile:
            name, info, score = parse_line(line, by_score)
            if name is False:
                continue
            if info[2][3] == '1' and secondary == 0: # neglect secondary alignments
                continue
            if name in ga_dic:
                passa = True
                compa = compare_sam_info(info, ga_dic[name], threshold)
            else:
                passa = False
            if name in gb_dic:
                passb = True
                compb = compare_sam_info(info, gb_dic[name], threshold)
            else:
                passb = False
            if __debug__ & 0:
                if comp is not True:
                    print (info)
                    print (ga_dic[name])
                    print (gb_dic[name])
                    input()
            if passa and compa:
                num_correct_reads_a += 1
            elif passa:
                num_incorrect_reads_a += 1
            if passb and compb:
                num_correct_reads_b += 1
            elif passb:
                num_incorrect_reads_b += 1


            if by_score > 0:
                if dic_as_cor_a.get(score) is None:
                    dic_as_cor_a[score] = 0
                if dic_as_incor_a.get(score) is None:
                    dic_as_incor_a[score] = 0
                if dic_as_cor_b.get(score) is None:
                    dic_as_cor_b[score] = 0
                if dic_as_incor_b.get(score) is None:
                    dic_as_incor_b[score] = 0

                if passa and compa:
                    dic_as_cor_a[score] += 1
                elif passa:
                    dic_as_incor_a[score] += 1
                if passb and compb:
                    dic_as_cor_b[score] += 1
                elif passb:
                    dic_as_incor_b[score] += 1
            if by_mapq > 0:
                mapq = int(info[3])
                if dic_q_cor_a.get(mapq) is None:
                    dic_q_cor_a[mapq] = 0
                if dic_q_incor_a.get(mapq) is None:
                    dic_q_incor_a[mapq] = 0

                if dic_q_cor_b.get(mapq) is None:
                    dic_q_cor_b[mapq] = 0
                if dic_q_incor_b.get(mapq) is None:
                    dic_q_incor_b[mapq] = 0

                if passa and compa:
                    dic_q_cor_a[mapq] += 1
                elif passa:
                    dic_q_incor_a[mapq] += 1

                if passb and compb:
                    dic_q_cor_b[mapq] += 1
                elif passb:
                    dic_q_incor_b[mapq] += 1

            # TODO
            if secondary > 0:
                #if dic_secondary.get(name) is None:
                #    dic_secondary[name] = 0
                if comp:
                    if dic_secondary.get(name) != 1:
                        dic_secondary[name] = 1
                        if __debug__:
                            print (info)
                            print (g_dic[name])
                            print ("SEC: precision +=1")
                            input()
                    # elif dic_secondary.get(name) == 1:
                    #     print (info)
                    #     print (g_dic[name])
    
    num_reads_a = num_correct_reads_a + num_incorrect_reads_a
    num_reads_b = num_correct_reads_b + num_incorrect_reads_b
    num_reads = num_reads_a + num_reads_b 
    #num_correct_reads_a + num_incorrect_reads_a + num_correct_reads_b + num_incorrect_reads_b

    if num_syn_reads != num_reads:
        print ('Warning: number of reads do not match!')
        print ('Synthesized reads = %d, correct reads = (%d, %d), incorrect reads = (%d, %d)' % \
                (num_syn_reads, num_correct_reads_a, num_correct_reads_b, num_incorrect_reads_a, num_incorrect_reads_b))
        # return 
    print ('Alignment sensitivity = %f (%d/%d)' % \
        (100 * float(num_correct_reads_a+num_correct_reads_b) / num_reads, (num_correct_reads_a+num_correct_reads_b), num_reads))
    print ('Hap A = %f (%d/%d)' % \
        (100 * float(num_correct_reads_a) / num_reads_a, num_correct_reads_a, num_reads_a))
    print ('Hap B = %f (%d/%d)' % \
        (100 * float(num_correct_reads_b) / num_reads_b, num_correct_reads_b, num_reads_b))
 
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
