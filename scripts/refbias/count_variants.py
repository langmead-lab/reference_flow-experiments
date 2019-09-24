import argparse
import pandas as pd
from utils import get_het_from_list_format

def count_single_hap_for_hets(fn_vcf, hap, threshold):
    '''
    hap [INT]:
        0 for hapA
        1 for hapB
    '''
    list_var = []
    list_size = []
    with open(fn_vcf, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            spl = line.split()
            ref = spl[3]
            assert ref.count(',') == 0
            alt = spl[4]
            assert alt.count(',') == 0
            gt = [int(i) for i in spl[-1].split('|')]
            assert gt[0] + gt[1] == 1

            list_var.append(int(spl[1]))
            if gt[0] and (hap == 0):
                list_size.append(abs(len(ref) - len(alt)) + 1)
            elif gt[1] and (hap == 1):
                list_size.append(abs(len(ref) - len(alt)) + 1)
            else:
                list_size.append(0)
    list_num_var = []
    for i, het in enumerate(list_var):
        var = 0
        for j in range(i - 1, -1, -1):
            if abs(list_var[j] - het) > threshold:
                break
            else:
                var += list_size[j]
        for j in range(i, len(list_var)):
            if abs(list_var[j] - het) > threshold:
                break
            else:
                var += list_size[j]
        list_num_var.append(var)

    return list_var, list_num_var

def count_size_of_variants_from_vcf(fn_vcf, indiv):
    #: assume only one chromosome in the vcf
    list_var = [] #: list of positions
    list_size = [] #: list of 2-element lists [sizeA, sizeB] for each variant
    # list_gt = []
    idx_indiv = None
    with open(fn_vcf, 'r') as f:
        for line in f:
            if line[0:2] == '##':
                continue
            elif line[0] == '#':
                header = line.split()
                for i, h in enumerate(header):
                    if h == indiv:
                        idx_indiv = i
                        break
                continue
            line = line.split()
            gt = line[idx_indiv].split('|')
            assert len(gt) == 2
            if gt[0] == '0' and gt[1] == '0':
                continue
            list_var.append(line[1])
            # list_gt.append(line[idx_indiv])
            ref_allele = line[3]
            #: ref allele should not have commas
            assert ref_allele.count(',') == 0
            alt_allele = line[4].split(',')
            size = [0, 0]
            for i, g in enumerate(gt):
                g = int(g)
                if g == 0:
                    size[i] = 0
                else:
                    size[i] = max(len(ref_allele), len(alt_allele[g - 1]))
            list_size.append(size)
    return list_var, list_size

def count_num_nearby_var_for_a_list_of_var(fn_vcf, indiv, threshold, list_target):
    list_var, list_size = count_size_of_variants_from_vcf(fn_vcf, indiv)

    list_var_near_target = [[0] * len(list_target), [0] * len(list_target)]
    for hap in [0, 1]:
        idx_target_start = 0
        for i_var, pos in enumerate(list_var):
            pos = int(pos)
            for i_target in range(idx_target_start, len(list_target)):
                if pos < list_target[i_target] - threshold:
                    break
                elif pos > list_target[i_target] + threshold:
                    idx_target_start = i_target + 1
                else:
                    list_var_near_target[hap][i_target] += list_size[i_var][hap]
            if idx_target_start == len(list_target):
                break
    #: take the max out of two haplotypes
    list_out = [max(list_var_near_target[0][i], list_var_near_target[1][i]) for i in range(len(list_target))]
    return list_out

def count_var(fn_vcf, indiv = '', target_het = None):
    # list_varA, list_num_var_A = count_single_hap_for_hets(fn_vcf, hap=0, threshold=100)
    # list_varB, list_num_var_B = count_single_hap_for_hets(fn_vcf, hap=1, threshold=100)
    assert list_varA == list_varB

    if target_het == None:
        return list_varA, list_num_var_A, list_num_var_B
    else:
        list_filtered_var = []
        list_filtered_num_var_A = []
        list_filtered_num_var_B = []
        for i, v in enumerate(list_varA):
            if v in target_het:
                list_filtered_var.append(v)
                list_filtered_num_var_A.append(list_num_var_A[i])
                list_filtered_num_var_B.append(list_num_var_B[i])
        return list_filtered_var, list_filtered_num_var_A, list_filtered_num_var_B

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file for chromosomes')
    parser.add_argument('-o', '--out', help = 'output file')
    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_output = args.out
    print('vcf', fn_vcf)
    print('output', fn_output)

    #count_single_hap(fn_vcf, '', 100, 'NA12878')
    count_num_nearby_var_for_a_list_of_var(
        fn_vcf = fn_vcf,
        indiv = 'NA12878',
        threshold = 100,
        list_target = [9411410, 9411500, 9411602],
        hap = 1)

    # list_var, list_numA, list_numB = count_var(fn_vcf)
    # if fn_out != None:
    #     df = pd.DataFrame(index = list_var)
    #     df['Num Var A'] = list_numA
    #     df['Num Var B'] = list_numB
    #     df.to_csv(fn_out, sep = '\t')

