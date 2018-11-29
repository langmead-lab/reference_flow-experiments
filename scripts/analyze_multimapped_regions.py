'''
Last update: 2018/11/28 by Nae-Chyun Chen

Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''

from analyze_sam import SamInfo, parse_line

# inputs
read_len = 100
sam_fn = '/scratch/groups/blangme2/naechyun/relaxing/alignments/diploid/na12878-chr9-AB2AB-bt2.sam'
diff_snp_fn = '/scratch/groups/blangme2/naechyun/relaxing/na12878/diff_AB.snp'

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

def analyze_mutimapped_regions():
    diff_snp_dic = build_diff_dic()
    with open(sam_fn, 'r') as sam_f:
        total_counts = 0
        correct_hap_counts = 0
        incorrect_hap_id_counts = 0
        incorrect_hap_nid_counts = 0
        unaligned_counts = 0
        for line in sam_f:
            name, info = parse_line(line, 0)
            if name is False:
                continue
            total_counts += 1
            # aligned to incorrect haplotype
            if info.is_unaligned():
                unaligned_counts += 1
            elif (name.find('hapA') > 0 and info.chrm != '9A') \
                or (name.find('hapB') > 0 and info.chrm != '9B'):
                identical_haplotypes = True
                for i in range(info.pos, info.pos + read_len):
                    if diff_snp_dic.get(i) != None:
                        identical_haplotypes = False
                        if __debug__:
                            print (incorrect_hap_nid_counts)
                            print ('name =', name)
                            info.print()
                            print (i, diff_snp_dic[i])
                            input ()
                        break
                # aligned to incorrect haplotype and two haps
                # are NOT equal
                if identical_haplotypes is False:
                    incorrect_hap_nid_counts += 1
                else:
                    incorrect_hap_id_counts += 1
            else:
                correct_hap_counts += 1
    print ('Correct hap:', correct_hap_counts)
    print ('Incorrect hap (identical):', incorrect_hap_id_counts) 
    print ('Incorrect hap (with varaints):', incorrect_hap_nid_counts)
    print ('Unaligned:', unaligned_counts)
    print ('Total:', total_counts)
    
if __name__ == '__main__':
    analyze_mutimapped_regions()
