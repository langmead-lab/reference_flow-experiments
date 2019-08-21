from sys import exit
import argparse
import os.path
from os import path

def main(fn_vcf, list_fn_sam, fn_fasta, list_fn_name, fn_out):
    import liftover_sam
    ref_bi_list = []
    counter = 1

    # name_list = ['MAJ-no-indel', 'AFR-no-indel', 'AMR-no-indel', 'EAS-no-indel', 'EUR-no-indel', 'SAS-no-indel']
    # name_list = ['perA-no-indel', 'perB-no-indel']
    for i in range(len(list_fn_sam)):
        print("i:", i)
        vcf = fn_vcf
        sam = list_fn_sam[i]
        fasta = fn_fasta
        
        print("vcf: ", vcf)
        print("sam: ", sam)
        print("fasta: ", fasta)
        
        # output = name_list[i] + '_refbias.txt'
        output = list_fn_sam[i] + '_refbias.txt'
        liftover_sam.main(vcf, sam, fasta, output)
        #cigar_whole_genome_sam.main(vcf, sam, fasta, output)
        ref_bi_list.append(output)
        counter += 1
    print("ref_bi_list: ", ref_bi_list)
    merge(ref_bi_list, fn_out)
    
def merge(ref_bi_list, fn_out):
    #question, how should i handle the HET sites when merging (b/c they are different depending on offsets
    het_site = []
    ref_count = []
    alt_count = []
    gap_count = []
    other_count = []
    num_reads = []
    chr = 21
    for i in range(len(ref_bi_list)):
        count = 0
        f_in = open(ref_bi_list[i], 'r') 
        for line in f_in:
            if '#' in line:
                continue
            else:
                #print("count: ", count)
                #print("line: ", line)
                spl = line.split()
                if i == 0:
                    if line == '\n':
                        het_site.append([])
                        ref_count.append(0)
                        alt_count.append(0)
                        gap_count.append(0)
                        other_count.append(0)
                        num_reads.append(0)
                    else:
                        het_site.append([int(spl[1])])
                        ref_count.append(int(spl[3]))
                        alt_count.append(int(spl[4]))
                        gap_count.append(int(spl[5]))
                        other_count.append(int(spl[6]))
                        num_reads.append(int(spl[7]))
                else: 
                    if line == '\n':
                        count += 1
                        continue
                    if count-1 < len(het_site): 
                        het_site[count].append(int(spl[1]))
                        ref_count[count] = ref_count[count] + int(spl[3])
                        alt_count[count] = alt_count[count] + int(spl[4])
                        gap_count[count] = gap_count[count] + int(spl[5])
                        other_count[count] = other_count[count] + int(spl[6])
                        num_reads[count] = num_reads[count] + int(spl[7])

            count += 1

    f_out = open(fn_out, 'w')
    f_out.write("#CHR\tHET SITE\tREFERENCE BIAS\tREF COUNT\tALT COUNT\tGAP COUNT\tOTHER COUNT\t# READS")
    f_out.write("\n")
    for i in range(len(ref_count)): 
        f_out.write(str(chr))
        f_out.write("\t")
        f_out.write(str(het_site[i]))
        f_out.write("\t")
        add = ref_count[i] + alt_count[i] + gap_count[i] + other_count[i]
        if add == 0:
            f_out.write("n/a")
        else:
            ref_bi = float(ref_count[i])/float(ref_count[i] + alt_count[i] + gap_count[i] + other_count[i])
            f_out.write(str(ref_bi))
        f_out.write("\t")
        f_out.write(str(ref_count[i]))
        f_out.write("\t")
        f_out.write(str(alt_count[i]))
        f_out.write("\t")
        f_out.write(str(gap_count[i]))
        f_out.write("\t")
        f_out.write(str(other_count[i]))
        f_out.write("\t")
        f_out.write(str(num_reads[i]))
        f_out.write("\n")

def get_paths_from_list(fn_list):
    l = []
    with open(fn_list, 'r') as f:
        for line in f:
            l.append(line.rstrip())
    return l

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-s', '--sam', help='list of sam file(s)')
    parser.add_argument('-n', '--name', help='list of partial output name (s)')
    parser.add_argument('-f', '--fasta', help='reference fasta file')
    parser.add_argument('-o', '--out', help='output joined reference bias files')

    args = parser.parse_args()

    fn_vcf = args.vcf
    list_fn_sam = get_paths_from_list(args.sam)
    list_fn_name = get_paths_from_list(args.name)
    fn_fasta = args.fasta
    fn_out = args.out

    print("vcf file: {}".format(fn_vcf))
    print("List of sam files: {}".format(list_fn_sam))
    print("fasta file: {}".format(fn_fasta))
    print("output: ", fn_out)
    main(
        fn_vcf=fn_vcf,
        list_fn_sam=list_fn_sam,
        fn_fasta=fn_fasta,
        list_fn_name=list_fn_name,
        fn_out=fn_out
    )
