import argparse
from utils import get_het_from_list_format

def main(fn_vcf, fn_het, fn_out):
    f_vcf = open(fn_vcf, 'r')

    het_list = []
    ref_alt_list = []
    hap_var_list = []

    for line in f_vcf:
        if not line.startswith("#"):
            spl = line.split()
            # het_list.append(get_het_from_list_format(spl[1]))
            het_list.append(int(spl[1]))
            alleleA = []
            alleleB = []
            if ',' in spl[3]:
                alleleA = spl[3].split(',')
            else:
                alleleA.append(spl[3])
            if ',' in spl[4]:
                alleleB = spl[4].split(',')
            else:
                alleleB.append(spl[4])
            alleles = (alleleA, alleleB)
            ref_alt_list.append(alleles)
            choices = spl[-1].split("|")
            hap_var_list.append((choices[0], choices[1]))

    #print(ref_alt_list)

    f_het = open(fn_het, 'r')
    f_out = open(fn_out, 'w')
    het_sites = []
    
    hapA_var_list = []
    hapB_var_list = []

    starting = 0
    change_starting = True
    count_line = 0
    for line in f_het:
        hapA = 0
        hapB = 0
        if count_line > 0:
            # het_site = int(line.split()[0])+1
            het_site = get_het_from_list_format(line.split('\t')[1]) + 1
            het_sites.append(het_site)
            ran = range(het_site-100, het_site+101)
            for i in range(starting, len(het_list)):
                choiceA = int(hap_var_list[i][0])
                choiceB = int(hap_var_list[i][1])
                if het_list[i] in ran and het_list[i] != het_site:
                    if change_starting:
                        starting = i
                        change_starting = False

                    if choiceA >= 1:
                        ref = ref_alt_list[i][0][0]
                        alt = ref_alt_list[i][1][choiceA-1]

                        if len(ref) == len(alt):
                            hapA += 1
                        else:
                            hapA += abs(len(ref) - len(alt))

                    if choiceB >= 1:
                        ref = ref_alt_list[i][0][0]
                        alt = ref_alt_list[i][1][choiceB-1]
                        if len(ref) == len(alt):
                            hapB += 1
                        else:
                            hapB += abs(len(ref) - len(alt))
                elif het_list[i] > het_site:
                    break
            hapA_var_list.append(hapA)
            hapB_var_list.append(hapB)
            change_starting = True

        count_line += 1

    f_out.write("HET SITE\tHapA\tHapB\t# Variants\n")
    for i in range(len(het_sites)):
        f_out.write(str(het_sites[i]))
        f_out.write("\t")
        ref_hap = find_ref_hap(het_sites[i], fn_vcf)
        if ref_hap == 'hapA':
            f_out.write("Yes\t\t")
            if hapA_var_list[i] > 10:
                f_out.write(str(10))
            else:
                f_out.write(str(hapA_var_list[i]+1))
        else:
            f_out.write("\tYes\t")
            if hapB_var_list[i] > 10:
                f_out.write(str(10))
            else:
                f_out.write(str(hapB_var_list[i]+1))
        f_out.write("\n")

def find_ref_hap(het_site, fn_vcf):
    file_in = open(fn_vcf, 'r')
    for line in file_in:
        if line.startswith("#"):
            continue
        else:
            spl = line.split()
            # if int(spl[1]) == het_site + 1: # because 0-base to 1-base
            if int(spl[1]) == het_site:#+1 because 0-base to 1-base
                hap = spl[8].split("|")
                if hap[0] == '0':
                    return 'hapA'
                else:
                    return 'hapB'
    print("error, het site not found")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file for chromosomes')
    parser.add_argument('-s', '--sit', help='file with HET sites')
    parser.add_argument('-o', '--out', help = 'output file')
    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_sit = args.sit
    fn_output = args.out
    print('vcf', fn_vcf)
    print('sit', fn_sit)
    print('output', fn_output)
    main(fn_vcf, fn_sit, fn_output)
