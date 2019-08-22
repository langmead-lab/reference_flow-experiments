'''
Reads a VCF (all variants should be HETs) and seperates it by chromosome
'''
import argparse

def main(fn_vcf):
    file = open(fn_vcf)
    header = []
    chr_num = []
    for line in file:
        if line.startswith("#"):
            header.append(line)
        else:
            spl = line.split()
            chr = int(spl[0])
            if chr not in chr_num:
                chr_num.append(chr)
                file_name = str(chr) + '_het.vcf'
                f = open(file_name, 'w+')
                for head in header:
                    f.write(head)
            f.write(line)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    args = parser.parse_args()
    fn_vcf = args.vcf
    #fn_output = args.out
    print('vcf', fn_vcf)
   # print('output', fn_output)
    main(fn_vcf)#, fn_sam, fn_output)
