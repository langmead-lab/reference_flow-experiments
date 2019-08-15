import argparse
import re
import pickle
import os.path
from os import path


def main(fn_vcf, fn_sam, fn_fasta, fn_output):
    ref_file = open(fn_fasta, 'r')
    reference = ''
    for line in ref_file:
        if not line.startswith('>'):
            reference += line.strip()


    #file = open('/net/langmead-bigmem-ib.bluecrab.cluster/storage/sheila/21_out.vcf', 'r')#open('major_21_thousand.vcf', 'r')#open('21_out.vcf', 'r')
    file = open(fn_vcf, 'r')
    count_line = 0
    index_ref = 0
    index_pos = 0
    index_alt = 0
    index_cat = 0

    vcf_file_name = fn_output+'.chr_vcf.pickle'
    sam_file_name = fn_output+'.chr_sam.pickle'

    f = open(fn_output, 'a+')

    count_het = 0
    if not path.exists(vcf_file_name):
        chr_vcf = {}
        for line in file:
            if line.startswith("##"):
                f.write(line.strip())
                f.write("\n")
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
                #print("line: ", line)
                if line == '\n':
                    if not chr_vcf:
                        chr_vcf[chr] = [[],[]]#this works since we are doing individual chromosomes
                    else:
                        chr_vcf[chr][0].append('N/A')
                        chr_vcf[chr][1].append(['N/A', 'N/A'])
                    continue
                chr = int(now[0])
                if not chr_vcf:
                    chr_vcf[chr] = [[], []]
                elif chr not in chr_vcf.keys():
                    chr_vcf[chr] = [[], []]

                if (len(now[index_alt]) == 1 and len(now[index_ref]) == 1):
                    chr_vcf[chr][0].append(int(now[index_pos]) - 1)
                    chr_vcf[chr][1].append([now[index_ref], now[index_alt]])
                elif ',' in now[index_alt]:
                    possib = now[index_alt].split(',')
                    chr_vcf[chr][0].append(int(now[index_pos]) - 1)
                    chr_vcf[chr][1].append([now[index_ref], ''.join(possib)])
            count_line += 1
        pickle.dump(chr_vcf, open(vcf_file_name, 'wb'))
        print ('Dump to {}'.format(vcf_file_name))
    else:
        print ('Load from {}'.format(vcf_file_name))
        chr_vcf = pickle.load(open(vcf_file_name, "rb"))
    #print("chr_vcf: ", chr_vcf[1]) 


    #file = open('/net/langmead-bigmem-ib.bluecrab.cluster/storage/sheila/real-NA12878/SRR622457-bt2-grch37.sam', 'r')#open('/scratch/groups/blangme2/naechyun/relaxing/chr21/experiments/10M-ref
    #/scratch/groups/blangme2/naechyun/test/SRR622457-bt2-chr21.sam
    file = open(fn_sam, 'r')
    #open('/net/langmead-bigmem-ib.bluecrab.cluster/storage/naechyun/1000Genomes/SRR622457/SRR622457_-chr21.sam', 'r')

    #f = open(fn_output, 'a+')
    f.write("#CHR\tHET SITE\tREFERENCE BIAS\tREF COUNT\tALT COUNT\tGAP COUNT\tOTHER COUNT\t# READS")
    f.write("\n")
    count_line = 0
    chr_sam = {}

    if not path.exists(sam_file_name):
        for line in file:
            if not line.startswith('@'):
                spl = line.split() 
                tag = int(spl[1])
                if (tag & 4):
                    continue
                chr = int(spl[2][0:2])
                cigar = spl[5]
                start_pos = int(spl[3]) - 1 #TODO
                sequence = spl[9]
                mod_sequence = ''
                #: ignores unmapped reads
                #: tag 4: unaligned
                if tag & 4:
                    continue
                #print(str(len(sequence)) + 'M')
                if not cigar == (str(len(sequence))+'M') and not (tag & 4):
                    #ref = reference[start_pos:start_pos + 102]
                    change = 0
                    start = start_pos
                    count_del = 0
                    count_ins = 0
                    for num1, idm in re.findall('(\d+)([IDMS])', cigar):
                        # print(start_pos)
                        if idm == 'M':
                            mod_sequence += sequence[change:change + int(num1)]
                        elif idm == 'D':
                            count_del += int(num1)
                            for i in range(int(num1)):
                                mod_sequence += '-'
                        elif idm == 'I':
                            count_ins += int(num1)
                        elif idm == 'S':
                            count_ins += int(num1)
                        else:
                            print ('error: unexpected cigar letter', cigar)
                            exit ()

                        if idm != 'D':
                            change += int(num1)
                    #print("tag: {0}, cigar: {1}".format(tag, cigar))
                    ref = reference[start_pos:start_pos + len(sequence) + count_del - count_ins]
                    try:
                        assert len(ref) == len(mod_sequence)
                    except:
                        if len(ref) == 0:
                            continue
                        else:
                            print("read name: ", spl[0])
                            print ("ref      seq", len(ref))
                            print ("modified seq", mod_sequence)
                            print ("start_pos", start_pos)
                            print ("count_ins", count_ins)
                            print ("count_del", count_del)
                            print ("cigar", cigar)
                            print ("flag", tag) 
                            exit()
                    #print("ref: ", ref)
                    #print("mod: ", mod_sequence)
                    
                    # input()
                else:
                    mod_sequence = sequence
                if not chr_sam:
                    chr_sam[chr] = [[], [], [], []]
                elif chr not in chr_sam.keys():
                    chr_sam[chr] = [[], [], [], []]
                chr_sam[chr][0].append(start_pos)#position
                chr_sam[chr][1].append(mod_sequence)#sequence
                chr_sam[chr][2].append(tag)
                chr_sam[chr][3].append(cigar)
            count_line += 1
        pickle.dump(chr_sam, open(sam_file_name, 'wb'))
        print ('Dump to {}'.format(sam_file_name))
    else:
        print ('Load from {}'.format(sam_file_name))
        chr_sam = pickle.load(open(sam_file_name, 'rb'))

    #pickle.dump(chr_vcf, open(fn_output+'.chr_vcf.pickle', 'w'))
    #pickle.dump(chr_sam, open(fn_output+'.chr_sam.pickle', 'w'))

    #print("chr_sam: ", chr_sam[1])

    chr_list = list(chr_vcf.keys())
    chr_list.sort()

    for chr in chr_list:
        het_site_list = chr_vcf[chr][0]
        options = chr_vcf[chr][1]
        sam_pos = chr_sam[chr][0]
        sam_reads = chr_sam[chr][1]

        have_started = False
        starting_point = 0
        count_pos = 0

        total_ref_count = 0
        total_alt_count = 0
        total_gap_count = 0
        total_other_count = 0

        for pos in het_site_list:
            ref_count = 0
            alt_count = 0
            gap_count = 0
            other_count = 0
            reads_at_het = 0
            if pos == 'N/A':
                f.write('\n')
                count_pos += 1
                continue
            for i in range(starting_point, len(sam_pos)):
                align = sam_pos[i]
                ran = range(align, align + len(sam_reads[i]))  # should this be 101, idts
                if pos in ran:
                    #print("here, it overlaps")
                    if not have_started:
                        have_started = True
                        starting_point = i
                    reads_at_het += 1
                    
                    '''#: to comment
                    allele = sam_reads[i][pos - align]
                    #if (allele not in options[count_pos][0]) and (sam_reads[i].count('N') <= 10):
                    # if (allele not in [options[count_pos][0], options[count_pos][1]]):
                    if (allele not in [options[count_pos][0]]):
                        ref_a = options[count_pos][0]
                        alt_a = options[count_pos][1]
                        print ("ref = {0}, alt = {1} @ {2}".format(ref_a, alt_a, pos+1))
                        print ('##########'+options[count_pos][0]+'###########')
                        print (reference[pos-10 : pos+12])
                        tmp_reads = '..........'+sam_reads[i]+'..........'
                        print (sam_reads[i][pos-align-10 : pos-align+12])
                        print (tmp_reads[pos-align : pos-align+22])
                        print ('##########'+allele+'###########')
                        print (chr_sam[chr][2][i], chr_sam[chr][1][i])
                        print ("cigar", chr_sam[chr][3][i])
                        input ()
                    if allele in options[count_pos][0]:
                        ref_count += 1
                    elif allele in options[count_pos][1]:
                        alt_count += 1
                    elif allele == '-':
                        gap_count += 1
                    else:
                        other_count += 1
                    '''
                    try:
                        allele = sam_reads[i][pos - align]
                        if allele in options[count_pos][0]:
                            ref_count += 1
                            total_ref_count += 1
                        elif allele in options[count_pos][1]:
                            alt_count += 1
                            total_alt_count += 1
                        elif allele == '-':
                            gap_count += 1
                            total_gap_count += 1
                        else:
                            other_count += 1
                            total_other_count += 1
                            if (pos in range(9560476-2000, 9560476+2000)):# or (pos in range(11083563-1000, 11083563+1000)):
                                print ('vcf at', pos)
                                print (sam_reads[i])
                                print (align)
                                print (chr_sam[chr][2][i])
                                print (chr_sam[chr][3][i])
                                print ('##########'+options[count_pos][0]+'###########')
                                print (reference[pos-10 : pos+12])
                                tmp_reads = '..........'+sam_reads[i]+'..........'
                                print (tmp_reads[pos-align : pos-align+22])
                                print ('##########'+allele+'###########')
                    except:
                        print("in here. pos is: ", pos, "   and align is: ", align)
                        #print("vcf[pos]: ", vcf[pos])
                        print("sam[align]: ", sam_reads[i])
                else:
                    if align > pos:
                        break
                    # print(ref_count, alt_count)
            have_started = False
            f.write(str(chr))
            f.write("\t")
            f.write(str(pos))
            f.write("\t")
            if ref_count == 0 and alt_count == 0:
                f.write("N/A")
            else:
                f.write(str(ref_count / float(ref_count + alt_count + gap_count + other_count)))
            f.write("\t")
            f.write(str(ref_count))
            f.write("\t")
            f.write(str(alt_count))
            f.write("\t")
            f.write(str(gap_count))
            f.write("\t")
            f.write(str(other_count))
            f.write("\t")
            f.write(str(reads_at_het))
            f.write("\n")
            count_pos += 1


    #ref_string = 'Total Ref Count:' + '/t' + str(total_ref_count)
    #alt_string = 'Total Alt Count:' + '/t' + str(total_alt_count)
    #gap_string = 'Total Gap Count:' + '/t' + str(total_gap_count)
    #other_string = 'Total Other Count:' + '/t' + str(total_other_count)

    #f.write(ref_string)
    #f.write("\n")
    #f.write(alt_string)
    #f.write("\n")
    #f.write(gap_string)
    #f.write("\n")
    #f.write(other_string)
    #f.write("\n")
    f.close()
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-f', '--fasta', help='reference fasta file')
    parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_sam = args.sam
    fn_fasta = args.fasta
    fn_output = args.out
    print('vcf', fn_vcf)
    print('sam', fn_sam)
    print('fasta', fn_fasta)
    print('output', fn_output)
    main(fn_vcf, fn_sam, fn_fasta, fn_output)
