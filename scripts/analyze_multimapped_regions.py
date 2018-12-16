'''
Last update: 2018/12/16 by Nae-Chyun Chen

Compares a diploid-to-diploid sam and 
checks if the multi-mapped regions are 
identical in two personalized refs
'''
import argparse
from analyze_sam import SamInfo, parse_line, load_golden_dic, compare_sam_info
from build_erg import build_erg

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
        '-r', '--read_len', type=int,
        default=100,
        help='read length [100]'
    )
    parser.add_argument(
        '-d', '--diff',
        help='the snp file specifying the differences between two haplotypes'
    )
    parser.add_argument(
        '--var',
        default=None,
        help='the file specifying the variants'
    )
    args = parser.parse_args()
    return args

def build_var_dic(var_fn):
    '''
    Build a dictionary for the .var file
    
    var file format:
    STRAND CHRM TYPE REFPOS ALTPOS REF ALT OFFSET

    key:
        STRAND_CHRM_ALTPOS (str)
    value:
        TYPE, REFPOS, REF, ALT, OFFSET (list)
    '''
    var_dic = {}
    with open(var_fn, 'r') as var_f:
        for line in var_f:
            line = line.split()
            strand = line[0]
            chrm = line[1]
            vtype = line[2]
            ref_pos = line[3]
            alt_pos = line[4]
            ref_allele = line[5]
            alt_allele = line[6]
            offset = line[7]
            c_offset = line[8]
            v_key = strand + '_' + chrm + '_' + alt_pos
            var_dic[v_key] = vtype, ref_pos, ref_allele, alt_allele, offset, c_offset
#            input (var_dic)
    return var_dic

def build_offset_index(var_fn):
    main_index, alt_index = \
        build_erg('', var_fn, hap_mode=1, f_len=100, mode='index') 
    print (len(main_index))
    print (len(alt_index))
    input ()

def build_diff_dic(diff_snp_fn):
    # build dictionary
    diff_snp_dic = {}
    with open(diff_snp_fn, 'r') as diff_snp_f:
        for line in diff_snp_f:
            if line.startswith('pos'):
                continue
            line = line.split()
            pos = int(line[0])
            diff_snp_dic[pos] = line[1], line[2]
#            input (diff_snp_dic)
    return diff_snp_dic

def diploid_compare(info, g_info, threshold, dip_flag):
    # don't check the other strand
    if dip_flag in ['c', 'n']:
        return compare_sam_info(info, g_info, threshold)
    # check the other strand
    elif dip_flag == 'i':
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
        print (' ')
        print ('golden')
        g_info.print(flag=False, mapq=False, score=False)
        input ()

def analyze_mutimapped_regions(args):
    sam_fn = args.sam
    golden_fn = args.golden
    threshold = args.threshold
    read_len = args.read_len
    diff_snp_fn = args.diff
    var_fn = args.var

    #var_dic = build_var_dic(var_fn)
    build_offset_index(var_fn)
    diff_snp_dic = build_diff_dic(diff_snp_fn)
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
            if name == 'header':
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
            # TODO: this condition should be more generic
            elif (name.find('hapA') > 0 and info.chrm != '9A') \
                or (name.find('hapB') > 0 and info.chrm != '9B'):
                identical_haplotypes = True
                v_id_hap = True
                if name.find('hapB') > 0:
                    seq_source = 'A_9_'
                else:
                    seq_source = 'B_9_'
                for i in range(info.pos, info.pos + read_len):
                    pos_on_ref = i + info.offset
                    k = seq_source + str(pos_on_ref)
                    print ('show pos_on_ref')
                    input (k)
                    if var_dic.get(k) != None:
                        v_id_hap = False
                        print ('var_dic', k, var_dic.get(k))
                    if diff_snp_dic.get(i) != None:
                        identical_haplotypes = False
                        print ('diff_dic', i, diff_snp_dic.get(i))
                        break
                print ('compare var_dic and diff_snp_dic')
                print (v_id_hap, identical_haplotypes)
                input ()
                # aligned to incorrect haplotype and two haps are NOT equal
                if identical_haplotypes == False:
                    comp = diploid_compare(info, golden_dic[name], threshold, 'n')
                    dip_flag = 'n'
                    incorrect_hap_nid_counts += 1
                    if comp == False: num_f_n += 1
                else:
                    dip_flag = 'i'
                    comp = diploid_compare(info, golden_dic[name], threshold, 'i')
                    incorrect_hap_id_counts += 1
                    if comp == False: num_f_i += 1
            else:
                dip_flag = 'c'
                comp = diploid_compare(info, golden_dic[name], threshold, 'c')
                correct_hap_counts += 1
                if comp == False:
                    num_f_c += 1
            if comp:
                num_tp += 1
            else:
                if __debug__:
                    print ('dip_flag =', dip_flag)
                    print ('info')
                    info.print(flag=False, mapq=False, score=False)
                    print ()
                    print ('golden')
                    golden_dic[name].print(flag=False, mapq=False, score=False)
                    input ('comp=%s\n' % comp)

    print ('\n------ Alignment category distribution ------')
    print ('Correct hap [c]:', correct_hap_counts)
    print ('Incorrect hap (identical) [i]:', incorrect_hap_id_counts) 
    print ('Incorrect hap (with varaints) [n]:', incorrect_hap_nid_counts)
    print ('Unaligned [u]:', unaligned_counts)
    print ('Total:', total_counts)

    print ('\n------ Alignment accuracy ------')
    print ('True: %d (%.2f%%)' % (num_tp, 100*float(num_tp)/total_counts))
    print ('False_c: %d (%.2f%%)' % (num_f_c, 100*float(num_f_c)/total_counts))
    print ('False_i: %d (%.2f%%)' % (num_f_i, 100*float(num_f_i)/total_counts))
    print ('False_n: %d (%.2f%%)' % (num_f_n, 100*float(num_f_n)/total_counts))
    print ('False_u: %d (%.2f%%)' % (unaligned_counts, 100*float(unaligned_counts)/total_counts))

if __name__ == '__main__':
    args = parse_args()
    analyze_mutimapped_regions(args)
