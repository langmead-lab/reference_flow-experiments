'''
Last update: 2018/11/29 by Nae-Chyun Chen

Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''

from analyze_sam import SamInfo, parse_line, load_golden_dic, compare_sam_info

# inputs
read_len = 100
threshold = 10
#sam_fn = '/scratch/groups/blangme2/naechyun/relaxing/alignments/diploid/na12878-chr9-AB2AB-bt2.sam'
sam_fn = '/scratch/groups/blangme2/naechyun/relaxing/alignments/diploid/processed/na12878-chr9-lowAB2AB-bt2.sam'
diff_snp_fn = '/scratch/groups/blangme2/naechyun/relaxing/na12878/diff_AB.snp'
golden_fn = '/scratch/groups/blangme2/naechyun/relaxing/syn_reads/diploid/na12878-chr9-phase3-1M.fq.sam'

def build_diff_dic():
    # build dictionary
    diff_snp_dic = {}
    with open(diff_snp_fn, 'r') as diff_snp_f:
        for line in diff_snp_f:
            if line.startswith('pos'):
                continue
            line = line.split()
            pos = int(line[0])
            diff_snp_dic[pos] = line[1], line[2]
    return diff_snp_dic

def diploid_compare(info, g_info, threshold, dip_flag):
    # don't check the other haplotype
    if dip_flag in ['c', 'n']:
        return compare_sam_info(info, g_info, threshold)
    # check the other haplotype
    elif dip_flag is 'i':
        info.chrm = '9A'
        comp1 = compare_sam_info(info, g_info, threshold)
        info.chrm = '9B'
        comp2 = compare_sam_info(info, g_info, threshold)
#        input ('i: %s' % (comp1 | comp2))
        return comp1 | comp2
    else:
        print ('Error: undistinguished dip_flag: %s' % dip_flag)
        return False

def show_info(name, info, g_info, dip_flag):
    if 1: #__debug__:
        print (name)
        print (dip_flag)
        print ('info')
        info.print(flag=False, mapq=False, score=False)
        print ('golden')
        g_info.print(flag=False, mapq=False, score=False)
        input ()

def analyze_mutimapped_regions():
    diff_snp_dic = build_diff_dic()
    golden_dic = load_golden_dic(golden_fn, 1)
    with open(sam_fn, 'r') as sam_f:
        total_counts = 0
        correct_hap_counts = 0
        incorrect_hap_id_counts = 0
        incorrect_hap_nid_counts = 0
        unaligned_counts = 0
        num_tp = 0
        num_f_n = 0
        num_f_i = 0
        num_f_c = 0
        for line in sam_f:
            name, info = parse_line(line, 0)
            # headers
            if name is False:
#                print (line[:line.find('\\')])
                continue
#            if info.mapq < 10:
#                print (line[:line.find('\\')])
#            continue
            total_counts += 1
            # aligned to incorrect haplotype
            if info.is_unaligned():
                unaligned_counts += 1
                dip_flag = 'u'
                comp = False
            elif (name.find('hapA') > 0 and info.chrm != '9A') \
                or (name.find('hapB') > 0 and info.chrm != '9B'):
                identical_haplotypes = True
                for i in range(info.pos, info.pos + read_len):
                    if diff_snp_dic.get(i) != None:
                        identical_haplotypes = False
                        '''if __debug__:
                            print (incorrect_hap_nid_counts)
                            print ('name =', name)
                            info.print()
                            print (i, diff_snp_dic[i])
                            input ()'''
                        break
                # aligned to incorrect haplotype and two haps are NOT equal
                if identical_haplotypes is False:
                    comp = diploid_compare(info, golden_dic[name], threshold, 'n')
                    dip_flag = 'n'
                    incorrect_hap_nid_counts += 1
                    if comp is False: num_f_n += 1
                else:
                    dip_flag = 'i'
                    comp = diploid_compare(info, golden_dic[name], threshold, 'i')
                    incorrect_hap_id_counts += 1
                    if comp is False: num_f_i += 1
            else:
                dip_flag = 'c'
                comp = diploid_compare(info, golden_dic[name], threshold, 'c')
                correct_hap_counts += 1
                if comp is False:
                    num_f_c += 1
                    #if total_counts > 1284737:
                    #    show_info(name, info, golden_dic[name], 'c')
            if comp:
                num_tp += 1
            else:
                if __debug__:
                    print (dip_flag)
                    print ('info')
                    info.print(flag=False, mapq=False, score=False)
                    print ('golden')
                    golden_dic[name].print(flag=False, mapq=False, score=False)
                    input ('comp=%s\n' % comp)

    print ('Correct hap:', correct_hap_counts)
    print ('Incorrect hap (identical):', incorrect_hap_id_counts) 
    print ('Incorrect hap (with varaints):', incorrect_hap_nid_counts)
    print ('Unaligned:', unaligned_counts)
    print ('Total:', total_counts)
    print ('#TruePos:', num_tp)
    print ('#False_c:', num_f_c)
    print ('#False_i:', num_f_i)
    print ('#False_n:', num_f_n)

if __name__ == '__main__':
    analyze_mutimapped_regions()
