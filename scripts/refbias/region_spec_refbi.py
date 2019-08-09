import argparse

def main(fn_inp, fn_vcf, fn_out):


    file = open(fn_inp, 'r')
    het_list = []
    ref_bi = []
    mapQ = []
    het_site = ""
    count_a = 0.0
    count_b = 0.0
    total_map_q = 0
    count_reads = 0
    for line in file:
        if 'HET SITE' in line:
            if het_list and het_site in het_list:
                if reference_hap == 'hapA':
                    ref_bi.append((count_a/(count_a + count_b)))
                else:
                    ref_bi.append(count_b/(count_a + count_b))
                mapQ.append(float(total_map_q / count_reads))
            het_site = int(line.strip().split()[2])
            reference_hap = find_ref_hap(het_site, fn_vcf)
            #print("het site = {}, Reference Hap = {}".format(het_site, reference_hap))
            het_list.append(het_site)
            count_a = 0.0
            count_b = 0.0
            total_map_q = 0
            #print("count reads: ", count_reads)
            count_reads = 0
        else:
            if int(line.strip().split()[3]) == 255:
                print("SOMETHING = 255")
            total_map_q += int(line.strip().split()[3])
            count_reads += 1
            if 'hapA' in line.strip().split()[0]:
                count_a += 1.0
            else:
                count_b += 1.0
    reference_bias = 0.0
    if reference_hap == 'hapA':
        reference_bias = count_a/(count_a + count_b)
    else:
        reference_bias = count_b / (count_a + count_b)
        print("reference bias: ", reference_bias)
    ref_bi.append(reference_bias)
    mapQ.append(float(total_map_q / count_reads))

    f_out = open(fn_out, 'w')
    string = "HET SITE" + "\t" + "REFERENCE BIAS" + "\t" + "Average MapQ"
    f_out.write(string)
    f_out.write("\n")
    for i in range(len(het_list)):
        f_out.write(str(het_list[i]))
        f_out.write("\t")
        f_out.write(str(ref_bi[i]))
        f_out.write("\t")
        f_out.write(str(mapQ[i]))
        f_out.write("\n")

def find_ref_hap(het_site, fn_vcf):
    file_in = open(fn_vcf, 'r')
    for line in file_in:
        if line.startswith("#"):
            continue
        else:
            spl = line.split()
            if int(spl[1]) == het_site+1:#+1 because 0-base to 1-base
                hap = spl[9].split("|")
                if hap[0] == '0':
                    return 'hapA'
                else:
                    return 'hapB'
    print("error, het site not found")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inp', help='input file with HET sites and overlapped reads')
    parser.add_argument('-v', '--vcf', help = 'vcf file')
    parser.add_argument('-o', '--out', help='output file')


    args = parser.parse_args()
    fn_inp = args.inp
    fn_vcf = args.vcf
    fn_out = args.out

    print("fn_inp: ", fn_inp)
    print("fn_vcf: ", fn_vcf)
    print("fn_out: ", fn_out)

    main(fn_inp, fn_vcf, fn_out)