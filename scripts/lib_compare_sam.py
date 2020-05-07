'''
We put common functions for SAM files comparison in this script
'''
import pickle
import os
import pandas as pd

class Summary:
    '''
    Summary of analysis
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
    Information of a SAM line
    '''
    pos = 0
    chrom = ''
    flag = 0
    mapq = 0
    score = 0
    offset = 0
    seq = ''
    cigar = ''
    tag_md = ''

    def __init__(self, line, erg=False, md=False, cigar=False, score=True):
        self.flag = int(line[1])
        self.pos = int(line[3])
        self.mapq = int(line[4])
        self.chrom = line[2]
        self.seq = line[9]
        if cigar:
            self.cigar = line[5]
        if erg and self.chrom.find('erg') > 0:
            tmp = self.chrom.split('-')
            self.offset = int(tmp[2])
            self.chrom = tmp[0]
            self.pos += self.offset - 1
        if md and (self.is_unaligned() == False):
            check_md = False
            for tag in line[11:]:
                if tag[:5] == 'MD:Z:':
                    self.tag_md = tag[5:]
                    check_md = True
                    break
            assert check_md
        if score:
            self.update_score(line[11])
    
    def print(self,
            pos=True, chrom=True,
            flag=True, mapq=True,
            score=True, offset=True):
        if pos:    print ('  pos =', self.pos)
        if chrom:  print ('  chrom =', self.chrom)
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

    def is_first_seg(self):
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

def parse_line(line, erg=False, md=False, cigar=False, mason2=False, score=True):
    if line[0] == '@':
        return 'header', False
    line = line.split()
    if mason2:
        # name = line[0].split('/')[-2]
        name = line[0][:line[0].rfind('/')]
    else:
        name = line[0]
    info = SamInfo(line, erg, md, cigar, score)
    return name, info

def compare_sam_info(
    info,
    ginfo,
    threshold,
    offset = [0],
    ignore_chrom=False
):
    '''
    Compare the info object from a SAM line with the gold standard

    Args:
    - info:
        info from alignment
    - ginfo:
        info from simulation profile (golden)
    - offset:
        positiontal offset
    - ignore_chrom:
        set True to ignore alignment against different chromosomes
    
    Returns:
        an INT representing alignment correctness
        if < 0:
            -1: unmatched chromosome
            -2: unmatched direction
        if >= 0:
            the difference in alignment position
            0 is a perfect match
    '''
    if (ignore_chrom is False) and (info.chrom != ginfo.chrom):
        #: diff chromosome
        if __debug__:
            print ("False: chr, mapq =", info.mapq)
        return -1
    if (info.is_rc() ^ ginfo.is_rc()) is True:
        #: diff direction
        if __debug__: 
            print ("False: direction (%s, %s)" % (info.is_rc(), ginfo.is_rc()), "mapq =", info.mapq)
        return -2
    return min([abs(info.pos + off - ginfo.pos) for off in offset])

def dump_golden_dic(filename):
    g_dic = {}
    with open(filename, 'r') as gfile:
        for line in gfile:
            name, info = parse_line(line)
            if name is 'header':
                continue
            if name in g_dic:
                assert g_dic[name].is_first_seg() != info.is_first_seg()
                if info.is_first_seg():
                    g_dic[name][0] = info
                else:
                    g_dic[name][1] = info
            else:
                if info.is_first_seg():
                    g_dic[name] = [info, None]
                else:
                    g_dic[name] = [None, info]

        print ("Size of database built:", len(g_dic))
    
    pkl_filename = filename + '.pkl'
    print ("Dump %s using pickle" % pkl_filename)
    f = open(pkl_filename, 'wb')
    pickle.dump(g_dic, f)
    f.close()
    return g_dic

def load_golden_dic(filename):
    pkl_filename = filename + '.pkl'
    if os.path.isfile(pkl_filename):
        f = open(pkl_filename, 'rb')
        g_dic = pickle.load(f)
        f.close()
        print ('Size of database loaded:', len(g_dic))
        return g_dic
    else:
        return dump_golden_dic(filename)

def print_df_stats(df, threshold, var_opt):
    '''
    Print summarized results

    Returns:
    - x: a list of FDR values
    - y: a list of sensitivity values
    '''

    print ()
    print ('--- Stats ---')
    unaligned = (df['dist'] == -3)
    correct = (df['dist'] >= 0) & (df['dist'] <= threshold)
    aligned_incorrect = (df['dist'] == -1) | (df['dist'] == -2)

    #: categories
    cat_same_id = (df['category'] == 'same_id')
    cat_same_var = (df['category'] == 'same_var')
    cat_diff_id = (df['category'] == 'diff_id')
    cat_diff_var = (df['category'] == 'diff_var')

    sensitivity_all = df[correct].shape[0] / df.shape[0]
    print ('sensitivity_all      = {0:.4%} ({1} / {2})'.format(sensitivity_all, df[correct].shape[0], df.shape[0]))
    precision_all = df[correct].shape[0] / (df.shape[0] - df[unaligned].shape[0])
    print ('precision_all        = {0:.4%} ({1} / {2})'.format(precision_all, df[correct].shape[0], df.shape[0] - df[unaligned].shape[0]))
    fdr_all = 1 - precision_all
    print ('fdr_all              = {0:.4%}'.format(fdr_all))
    unaligned_rate = df[unaligned].shape[0] / df.shape[0]
    print ('unaligned            = {0:.4%} ({1} / {2})'.format(unaligned_rate, df[unaligned].shape[0], df.shape[0]))
    
    print ()
    try:
        sensitivity_same_id = df[correct & cat_same_id].shape[0] / df[cat_same_id].shape[0]
        print ('sensitivity_same_id  = {0:.4%} ({1} / {2})'.format(sensitivity_same_id, df[correct & cat_same_id].shape[0], df[cat_same_id].shape[0]))
        fnr_same_id = 1 - sensitivity_same_id
        print ('fnr_same_id          = {0:.4%} ({1} / {2})'.format(fnr_same_id, df[cat_same_id].shape[0] - df[correct & cat_same_id].shape[0], df[cat_same_id].shape[0]))
    except:
        print ('Warning: no element in "same_id"')
    try:
        sensitivity_same_var = df[correct & cat_same_var].shape[0] / df[cat_same_var].shape[0]
        print ('sensitivity_same_var = {0:.4%} ({1} / {2})'.format(sensitivity_same_var, df[correct & cat_same_var].shape[0], df[cat_same_var].shape[0]))
        fnr_same_var = 1 - sensitivity_same_var
        print ('fnr_same_var         = {0:.4%} ({1} / {2})'.format(fnr_same_var, df[cat_same_var].shape[0] - df[correct & cat_same_var].shape[0], df[cat_same_var].shape[0]))
    except:
        print ('Warning: no element in "same_var"')
    try:
        sensitivity_diff_id = df[correct & cat_diff_id].shape[0] / df[cat_diff_id].shape[0]
        print ('sensitivity_diff_id  = {0:.4%} ({1} / {2})'.format(sensitivity_diff_id, df[correct & cat_diff_id].shape[0], df[cat_diff_id].shape[0]))
        fnr_diff_id = 1 - sensitivity_diff_id
        print ('fnr_diff_id          = {0:.4%} ({1} / {2})'.format(fnr_diff_id, df[cat_diff_id].shape[0] - df[correct & cat_diff_id].shape[0], df[cat_diff_id].shape[0]))
    except:
        print ('Warning: no element in "diff_id"')
    try:
        sensitivity_diff_var = df[correct & cat_diff_var].shape[0] / df[cat_diff_var].shape[0]
        print ('sensitivity_diff_var = {0:.4%} ({1} / {2})'.format(sensitivity_diff_var, df[correct & cat_diff_var].shape[0], df[cat_diff_var].shape[0]))
        fnr_diff_var = 1 - sensitivity_diff_var
        print ('fnr_diff_var         = {0:.4%} ({1} / {2})'.format(fnr_diff_var, df[cat_diff_var].shape[0] - df[correct & cat_diff_var].shape[0], df[cat_diff_var].shape[0]))
    except:
        print ('Warning: no element in "diff_var"')

    #: number of overlapping variants
    print ()
    var_all = (df['numvar'] >= 0)
    var0 = (df['numvar'] == 0)
    var1 = (df['numvar'] == 1)
    var2 = (df['numvar'] == 2)
    var3plus = (df['numvar'] >= 3)
    
    if var_opt == 'all': 
        v_filter = var_all
    elif var_opt == '0':
        v_filter = var0
    elif var_opt == '1':
        v_filter = var1
    elif var_opt == '2':
        v_filter = var2
    elif var_opt == '3+':
        v_filter = var3plus

    try:
        sensitivity_var0 = df[correct & var0].shape[0] / df[var0].shape[0]
        print ('sensitivity_var0     = {0:.4%} ({1} / {2})'.format(sensitivity_var0, df[correct & var0].shape[0], df[var0].shape[0]))
    except:
        print ('Warning: no read has 0 variants')
    try:
        sensitivity_var1 = df[correct & var1].shape[0] / df[var1].shape[0]
        print ('sensitivity_var1     = {0:.4%} ({1} / {2})'.format(sensitivity_var1, df[correct & var1].shape[0], df[var1].shape[0]))
    except:
        print ('Warning: no read has 1 variants')
    try:
        sensitivity_var2 = df[correct & var2].shape[0] / df[var2].shape[0]
        print ('sensitivity_var2     = {0:.4%} ({1} / {2})'.format(sensitivity_var2, df[correct & var2].shape[0], df[var2].shape[0]))
    except:
        print ('Warning: no read has 2 variants')
    try:
        sensitivity_var3plus = df[correct & var3plus].shape[0] / df[var3plus].shape[0]
        print ('sensitivity_var3plus = {0:.4%} ({1} / {2})'.format(sensitivity_var3plus, df[correct & var3plus].shape[0], df[var3plus].shape[0]))
    except:
        print ('Warning: no read has 3+ variants')
    
    #: mapq
    print ()
    mapq5plus = (df['mapq'] >= 5)
    mapq10plus = (df['mapq'] >= 10)
    mapq20plus = (df['mapq'] >= 20)
    mapq30plus = (df['mapq'] >= 30)
    mapq40plus = (df['mapq'] >= 40)

    sensitivity_mapqall = df[correct & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapqall    = {0:.4%} ({1} / {2})'.format(sensitivity_mapqall, df[correct & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq5plus = df[correct & mapq5plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq5plus  = {0:.4%} ({1} / {2})'.format(sensitivity_mapq5plus, df[correct & mapq5plus & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq10plus = df[correct & mapq10plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq10plus = {0:.4%} ({1} / {2})'.format(sensitivity_mapq10plus, df[correct & mapq10plus & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq20plus = df[correct & mapq20plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq20plus = {0:.4%} ({1} / {2})'.format(sensitivity_mapq20plus, df[correct & mapq20plus & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq30plus = df[correct & mapq30plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq30plus = {0:.4%} ({1} / {2})'.format(sensitivity_mapq30plus, df[correct & mapq30plus & v_filter].shape[0], df[v_filter].shape[0]))
    sensitivity_mapq40plus = df[correct & mapq40plus & v_filter].shape[0] / df[v_filter].shape[0]
    print ('sensitivity_mapq40plus = {0:.4%} ({1} / {2})'.format(sensitivity_mapq40plus, df[correct & mapq40plus & v_filter].shape[0], df[v_filter].shape[0]))
    
    print ()
    try:
        precision_mapqall = df[correct & v_filter].shape[0] / df[v_filter].shape[0]
        print ('precision_mapqall    = {0:.4%} ({1} / {2})'.format(precision_mapqall, df[correct & v_filter].shape[0], df[v_filter].shape[0]))
        # fdr_mapqall = 1 - precision_mapqall
        # print ('fdr_mapqall          = {0:.4%}'.format(fdr_mapqall))
    except:
        print ('Warning: no read is 0+ mapq')
        fdr_mapqall = 1
    try:
        precision_mapq5plus = df[correct & mapq5plus & v_filter].shape[0] / df[mapq5plus & v_filter].shape[0]
        print ('precision_mapq5plus  = {0:.4%} ({1} / {2})'.format(precision_mapq5plus, df[correct & mapq5plus & v_filter].shape[0], df[mapq5plus & v_filter].shape[0]))
        # fdr_mapq5plus = 1 - precision_mapq5plus
        # print ('fdr_mapq5plus        = {0:.4%}'.format(fdr_mapq5plus))
    except:
        print ('Warning: no read is 5+ mapq')
        fdr_mapq5plus = 1
    try:
        precision_mapq10plus = df[correct & mapq10plus & v_filter].shape[0] / df[mapq10plus & v_filter].shape[0]
        print ('precision_mapq10plus = {0:.4%} ({1} / {2})'.format(precision_mapq10plus, df[correct & mapq10plus & v_filter].shape[0], df[mapq10plus & v_filter].shape[0]))
        # fdr_mapq10plus = 1 - precision_mapq10plus
        # print ('fdr_mapq10plus       = {0:.4%}'.format(fdr_mapq10plus))
    except:
        print ('Warning: no read is 10+ mapq')
        fdr_mapq10plus = 1
    try:
        precision_mapq20plus = df[correct & mapq20plus & v_filter].shape[0] / df[mapq20plus & v_filter].shape[0]
        print ('precision_mapq20plus = {0:.4%} ({1} / {2})'.format(precision_mapq20plus, df[correct & mapq20plus & v_filter].shape[0], df[mapq20plus & v_filter].shape[0]))
        # fdr_mapq20plus = 1 - precision_mapq20plus
        # print ('fdr_mapq20plus       = {0:.4%}'.format(fdr_mapq20plus))
    except:
        print ('Warning: no read is 20+ mapq')
        fdr_mapq20plus = 1
    try:
        precision_mapq30plus = df[correct & mapq30plus & v_filter].shape[0] / df[mapq30plus & v_filter].shape[0]
        print ('precision_mapq30plus = {0:.4%} ({1} / {2})'.format(precision_mapq30plus, df[correct & mapq30plus & v_filter].shape[0], df[mapq30plus & v_filter].shape[0]))
        # fdr_mapq30plus = 1 - precision_mapq30plus
        # print ('fdr_mapq30plus       = {0:.4%}'.format(fdr_mapq30plus))
    except:
        print ('Warning: no read is 30+ mapq')
        fdr_mapq30plus = 1
    try:
        precision_mapq40plus = df[correct & mapq40plus & v_filter].shape[0] / df[mapq40plus & v_filter].shape[0]
        print ('precision_mapq40plus = {0:.4%} ({1} / {2})'.format(precision_mapq40plus, df[correct & mapq40plus & v_filter].shape[0], df[mapq40plus & v_filter].shape[0]))
        # fdr_mapq40plus = 1 - precision_mapq40plus
        # print ('fdr_mapq40plus       = {0:.4%}'.format(fdr_mapq40plus))
    except:
        print ('Warning: no read is 40+ mapq')
        fdr_mapq40plus = 1
    
    x = [fdr_mapqall, fdr_mapq5plus, fdr_mapq10plus, fdr_mapq20plus, fdr_mapq30plus, fdr_mapq40plus]
    y = [sensitivity_mapqall, sensitivity_mapq5plus, sensitivity_mapq10plus, sensitivity_mapq20plus, sensitivity_mapq30plus, sensitivity_mapq40plus]

    return x, y

