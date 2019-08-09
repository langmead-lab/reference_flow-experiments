import argparse
import time

def main(
    fn_vcf,
    fn_sam,
    fn_output
):
    #file = open('/net/langmead-bigmem-ib.bluecrab.cluster/storage/sheila/21_out.vcf', 'r')#open('major_21_thousand.vcf', 'r')#open('21_out.vcf', 'r')
    
    t = time.time()
    
    file = open(fn_vcf, 'r')#open('major_21_thousand.vcf', 'r')#open('21_out.vcf', 'r')`

    count_line = 0
    index_ref = 0
    index_pos = 0
    index_alt = 0
    index_cat = 0

    vcf = {}

    for line in file:
        if not line.startswith("##") and line.startswith("#"):
            categories = line.split()
            index_cat = count_line
            for i in range(len(categories)):
                if categories[i] == 'REF':
                    index_ref = i
                elif categories[i] == 'ALT':
                    index_alt = i
                elif categories[i] == 'POS':
                    index_pos = i
        elif count_line > index_cat and index_cat != 0:
            now = line.split()
            if (len(now[index_alt]) == 1 and len(now[index_ref]) == 1):
                vcf[int(now[index_pos])] = [now[index_ref], now[index_alt]]
            elif ',' in now[index_alt]:
                possib = now[index_alt].split(',')
                vcf[int(now[index_pos])] = [now[index_ref], ''.join(possib)]
        count_line += 1

        #ifopen count_line > 1500:
         #   break
    #print("vcf: ", vcf)

    #file = open('/net/langmead-bigmem-ib.bluecrab.cluster/storage/sheila/real-NA12878/SRR622457-bt2-grch37.sam', 'r')#open('/scratch/groups/blangme2/naechyun/relaxing/chr21/experiments/10M-ref.sam', 'r')#open('chr21-10M-h37maj.sam', 'r')#open('chr21-h37maj-partial.sam', 'r')#open('/scratch/groups/blangme2/naechyun/test/SRR622457-bt2-chr21.sam', 'r')
    #/scratch/groups/blangme2/naechyun/test/SRR622457-bt2-chr21.sam
    file = open(fn_sam, 'r')
    #open('/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/1000Genomes/SRR622457/SRR622457_-chr21.sam', 'r')

    #f = open('full_bowtie_ref_bi.txt', 'w')#open('10M_ref_reference_bias.txt', 'w')#open('reference_bias.txt', 'w')
    f = open(fn_output, 'w')
    f.write("HET SITE\tREFERENCE BIAS\tREF COUNT\tALT COUNT\t# READS")
    f.write("\n")
    count_line = 0
    sam_pos = []
    sam_reads = []
    
    
    for line in file:
        if not line.startswith('@'):
            spl = line.split()
            #sam[int(spl[3])] = spl[9]
            sam_pos.append(int(spl[3]))
            sam_reads.append(spl[9])
        count_line += 1

        #if count_line > 1500:
        #    break
    #print("sam:  ", sam)
    #ref_bias = []
    het_site_list = list(vcf.keys())
    #print("before sort: ", het_site_list[:100])
    #het_site_list = het_site_list.sort()
    #print("after sort: ", het_site_list)
    het_site_list.sort()
    #print("after sort: ", het_site_list[:100])

    #print("sam_pos: ", sam_pos[:100])
    #print("sam_reads: ", sam_reads[:10])

    #print (vcf.keys())
    #print("sam keys: ", sam.keys())
    for pos in het_site_list:
        #input(("pos: ", pos))
        ref_count = 0
        alt_count = 0
        reads_at_het = 0
        #input (('------vcf', pos))
        for i in range(len(sam_pos)): 
            align = sam_pos[i]   
            ran = range(align, align+len(sam_reads[i]))#should this be 101, idts
            if pos in ran:
                reads_at_het += 1
                try:
                    #print("align: ", align)
                    #print("sam[align]: ", sam[align])
                    allele = sam_reads[i][pos-align]
                    #print("allele: ", allele)
                    #print("    string: ", sam[align][0:pos-align+1])
                    if allele in vcf[pos][0]:
                        ref_count += 1
                    elif allele in vcf[pos][1]:
                        alt_count += 1
                except:
                    print("in here. pos is: ", pos, "   and align is: ", align)
                    print("vcf[pos]: ", vcf[pos])
                    print("sam[align]: ", sam_reads[i])
            else:
                if align > pos: 
                    break
        #print(ref_count, alt_count)
        f.write(str(pos))
        f.write("\t")
        if ref_count == 0 and alt_count == 0:
            f.write("N/A")
        else:
            f.write(str(ref_count/float(ref_count + alt_count)))
        f.write("\t")
        f.write(str(ref_count))
        f.write("\t")
        f.write(str(alt_count))
        f.write("\t")
        f.write(str(reads_at_het))
        f.write("\n")

    print("time for program to run: ", time.time() - t)
    f.close()
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_sam = args.sam
    fn_output = args.out
    print ('vcf', fn_vcf)
    print ('sam', fn_sam)
    print ('output', fn_output)
    main(fn_vcf, fn_sam, fn_output)
