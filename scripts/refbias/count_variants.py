import argparse
import pandas as pd
from utils import get_het_from_list_format

def count_single_hap(fn_vcf, hap, threshold):
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

def count_var(fn_vcf, fn_out, target_het=None):
    list_varA, list_num_var_A = count_single_hap(fn_vcf, hap=0, threshold=100)
    list_varB, list_num_var_B = count_single_hap(fn_vcf, hap=1, threshold=100)
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

    list_var, list_numA, list_numB = count_var(fn_vcf, fn_output)
    if fn_out != None:
        df = pd.DataFrame(index = list_var)
        df['Num Var A'] = list_numA
        df['Num Var B'] = list_numB
        df.to_csv(fn_out, sep = '\t')

