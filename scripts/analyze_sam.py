'''
Analyzes aligned sam with synthetic golden data
and measures sensitivity
'''
import argparse
import pickle
import os

class Summary:
    '''
    Records the summary of analysis
    '''
    num_total = 0
    num_same_strand = 0
    num_same_id = 0
    num_same_var = 0
    num_diff_id = 0
    num_diff_var = 0
    num_unaligned = 0
    # If golden is specified:
    has_answer = True
    num_true_pos = 0
    num_false_ss = 0    
    num_false_sid = 0
    num_false_svar = 0
    num_false_id = 0
    num_false_var = 0

    def __init__(self, has_answer):
        if has_answer == False:
            self.num_true_pos = None
            self.num_false_ss = None   
            self.num_false_id = None
            self.num_false_var = None
            self.has_answer = False

    def add_one(self):
        self.num_total += 1

    def add_by_categories(self, flag, comp):
        if flag == 'unaligned':
            self.num_unaligned += 1
            return
        
        if comp and self.has_answer:
            self.num_true_pos += 1
        
        if flag == 'diff_id':
            self.num_diff_id += 1
            if comp == False and self.has_answer:
                self.num_false_id += 1
        elif flag == 'diff_var':
            self.num_diff_var += 1
            if comp == False and self.has_answer:
                self.num_false_var += 1
        elif flag == 'same_strand':
            self.num_same_strand += 1
            if comp == False and self.has_answer:
                self.num_false_ss += 1
        elif flag == 'same_id':
            self.num_same_id += 1
            if comp == False and self.has_answer:
                self.num_false_sid += 1
        elif flag == 'same_var':
            self.num_same_var += 1
            if comp == False and self.has_answer:
                self.num_false_svar += 1
    
    def show_summary(self, has_answer):
        frac_ss_total = 100 * float(self.num_same_strand) / self.num_total
        frac_id_total = 100 * float(self.num_diff_id) / self.num_total
        frac_v_total = 100 * float(self.num_diff_var) / self.num_total
        frac_sid_total = 100 * float(self.num_same_id) / self.num_total
        frac_sv_total = 100 * float(self.num_same_var) / self.num_total
        frac_v_total = 100 * float(self.num_diff_var) / self.num_total
        frac_un_total = 100 * float(self.num_unaligned) / self.num_total
        print ('\n------ Alignment categories ------')
        print ('Same haplotype: %d (%.2f%%)' % (self.num_same_strand, frac_ss_total))
        print ('Same haplotype, read has no variants: %d (%.2f%%)' % (self.num_same_id, frac_sid_total))
        print ('Same haplotype, read with variant(s): %d (%.2f%%)' % (self.num_same_var, frac_sv_total))
        print ('Diff haplotype, read has no variants: %d (%.2f%%)' % (self.num_diff_id, frac_id_total))
        print ('Diff haplotype, read with variant(s): %d (%.2f%%)' % (self.num_diff_var, frac_v_total))
        print ('Unaligned: %d (%.2f%%)' % (self.num_unaligned, frac_un_total))
        print ('Total:', self.num_total)

        if has_answer:
            print ('\n------ Alignment accuracy ------')
            frac_tp_total = 100 * float(self.num_true_pos) / self.num_total
            frac_fss_total = 100 * float(self.num_false_ss) / self.num_total
            if self.num_same_strand == 0:
                print ('Warning: num_same_strand = 0')
                frac_fss_ss = 0
            else:
                frac_fss_ss = 100 * float(self.num_false_ss) / self.num_same_strand
            frac_fsid_total = 100 * float(self.num_false_sid) / self.num_total
            if self.num_same_id == 0:
                print ('Warning: num_same_id = 0')
                frac_fsid_sid = 0
            else:
                frac_fsid_sid = 100 * float(self.num_false_sid) / self.num_same_id
            frac_fsvar_total = 100 * float(self.num_false_svar) / self.num_total
            if self.num_same_var == 0:
                print ('Warning: num_same_var = 0')
                frac_fsvar_svar = 0
            else:
                frac_fsvar_svar = 100 * float(self.num_false_svar) / self.num_same_var
            frac_fid_total = 100 * float(self.num_false_id) / self.num_total
            if self.num_diff_id == 0:
                print ('Warning: num_diff_id = 0')
                frac_fid_id = 0
            else:
                frac_fid_id = 100 * float(self.num_false_id) / self.num_diff_id
            frac_fv_total = 100 * float(self.num_false_var) / self.num_total
            if self.num_diff_var == 0:
                print ('Warning: num_diff_var = 0')
                frac_fv_v = 0
            else:
                frac_fv_v = 100 * float(self.num_false_var) / self.num_diff_var

            print ('True            : %d (total=%.2f%%)' % (self.num_true_pos, frac_tp_total))
            print ('Unaligned       : %d (total=%.2f%%)' % (self.num_unaligned, frac_un_total))
            print ('Same     F/TOTAL: %d / %d (total=%.2f%%, cat=%.2f%%)' % (self.num_false_ss, self.num_same_strand, frac_fss_total, frac_fss_ss))
            print ('Same-id  F/TOTAL: %d / %d (total=%.2f%%, cat=%.2f%%)' % (self.num_false_sid, self.num_same_id, frac_fsid_total, frac_fsid_sid))
            print ('Same-var F/TOTAL: %d / %d (total=%.2f%%, cat=%.2f%%)' % (self.num_false_svar, self.num_same_var, frac_fsvar_total, frac_fsvar_svar))
            print ('Diff-id  F/TOTAL: %d / %d (total=%.2f%%, cat=%.2f%%)' % (self.num_false_id, self.num_diff_id, frac_fid_total, frac_fid_id))
            print ('Diff-var F/TOTAL: %d / %d (total=%.2f%%, cat=%.2f%%)' % (self.num_false_var, self.num_diff_var, frac_fv_total, frac_fv_v))

class SamInfo:
    '''
    Records information of a sam line
    '''
    pos = 0
    chrm = ''
    flag = 0
    mapq = 0
    score = 0
    offset = 0
    seq = ''
    cigar = ''
    tag_md = ''

    def __init__(self, line, erg=False, md=False, cigar=False):
        self.flag = int(line[1])
        self.pos = int(line[3])
        self.mapq = int(line[4])
        self.chrm = line[2]
        self.seq = line[9]
        if cigar:
            self.cigar = line[5]
        if erg and self.chrm.find('erg') > 0:
            tmp = self.chrm.split('-')
            self.offset = int(tmp[2])
            self.chrm = tmp[0]
            self.pos += self.offset - 1
        if md and (self.is_unaligned() == False):
            check_md = False
            for tag in line[11:]:
                if tag[:5] == 'MD:Z:':
                    self.tag_md = tag[5:]
                    check_md = True
                    break
            assert check_md
        self.update_score(line[11])
    
    def print(self,
            pos=True, chrm=True,
            flag=True, mapq=True,
            score=True, offset=True):
        if pos:    print ('  pos =', self.pos)
        if chrm:   print ('  chrm =', self.chrm)
        if flag:   print ('  flag =', self.flag)
        if mapq:   print ('  mapq =', self.mapq)
        if score:  print ('  score =', self.score)
        if offset: print ('  offset =', self.offset)
        return

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
        if self.is_unaligned():
            self.score = 1
            return
        elif raw_score.startswith('AS:') is False:
            self.score = 1
            return
        self.score = int(raw_score.split(':')[-1])


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

def parse_line(line, erg=False, md=False, cigar=False, mason2=False):
    if line[0] == '@':
        return 'header', False
    line = line.split()
    if mason2:
        name = line[0][:line[0].find('/')]
    else:
        name = line[0]
    info = SamInfo(line, erg, md, cigar)
    return name, info

def compare_sam_info(info, ginfo, threshold, offset = [0], ignore_chrm=False):
    '''
    Inputs:
        info:
            info from alignment
        ginfo:
            info from simulation profile (golden)
        offset:
            positiontal offset
        ignore_chrm:
            set True to ignore alignment against different chromosomes
    
    Output:
        an INT representing alignment correctness
        if < 0:
            -1: unmatched chromosome
            -2: unmatched direction
        if >= 0:
            the difference in alignment position
            0 is a perfect match
    '''
    if (ignore_chrm is False) and (info.chrm != ginfo.chrm):
        #: diff chromosome
        if __debug__:
            print ("False: chr, mapq =", info.mapq)
        # return False
        return -1
    if (info.is_rc() ^ ginfo.is_rc()) is True:
        #: diff direction
        if __debug__: 
            print ("False: direction (%s, %s)" % (info.is_rc(), ginfo.is_rc()), "mapq =", info.mapq)
        # return False
        return -2
    # for off in offset:
    #     if abs(info.pos + off - ginfo.pos) <= threshold:
    #         if __debug__:
    #             print ("True, mapq =", info.mapq)
    #         return True
    # if __debug__: 
    #     print ("False: distance > threshold, mapq =", info.mapq)
    # return False
    return min([abs(info.pos + off - ginfo.pos) for off in offset])

def dump_golden_dic(filename, seg):
    g_dic = {}
    with open(filename, 'r') as gfile:
        for line in gfile:
            name, info = parse_line(line)
            if name is 'header':
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
    # haploid = args.haploid
    by_score = args.by_score
    by_mapq = args.by_mapq
    secondary = args.secondary

    num_correct_reads = 0
    num_incorrect_reads = 0
    dic_secondary = {}
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
            name, info = parse_line(line)
            if name is 'header':
                continue
            if info.is_secondary() and secondary == 0: # neglect secondary alignments
                continue
            dist = compare_sam_info(info, g_dic[name], threshold)
            if dist < 0 or dist > threshold:
                comp = False
            else:
                comp = True
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
