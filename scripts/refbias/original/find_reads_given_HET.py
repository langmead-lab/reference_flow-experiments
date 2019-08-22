import re
import os.path
import argparse
#import pandas as pd

def main(fn_sam, fn_fil, fn_out):#, fn_merge):
    
    #f_out = open(fn_out, 'w')
    for i in range(len(fn_sam)):
        element = fn_sam[i]
        #if os.path.exists(fn_out[i]):
        #    continue
        f_out = open(fn_out[i], 'w')
        print(element)
        sam_file = open(element, 'r')
        sam_name = []
        sam_reads = []
        sam_offset = []
        sam_mapQ = []
        for line in sam_file:
            if not line.startswith('@'):
                spl = line.split()
                start_pos = int(spl[3])-1
                tag =int(spl[1])
                if (tag & 4):
                    continue
                chr = int(spl[2])
                cigar = spl[5]
                #start_pos = int(spl[3]) - 1 #T
                sequence = spl[9]
                mod_sequence = ''
                if tag & 4:
                    continue  
            
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
                else:
                    mod_sequence = sequence

                sam_offset.append(start_pos)
                sam_reads.append(mod_sequence)
                sam_name.append(spl[0])
                sam_mapQ.append(spl[4])
            #print("just finished a sam file")
        starting = 0

        file = open(fn_fil, 'r')
        count_line = 0
        for line in file:
            have_started = False
            if count_line >= 2:
                het_pos = int(line.strip())
                toWrite = "HET SITE " + str(het_pos)
                f_out.write(toWrite)
                f_out.write("\n")
                for i in range(starting, len(sam_name)):
                    align = sam_offset[i]
                    ran = range(align, align + len(sam_reads[i]))
                    if het_pos in ran:
                        if not have_started:
                            have_started = True
                            starting = i
                        toWrite = sam_name[i] + "\t" + str(align) + "\t" + sam_reads[i] + "\t" + sam_mapQ[i] + "\n"
                        f_out.write(toWrite)
                    else:
                        if align > het_pos:
                            #print("i broke")
                            break
                    #print("starting: ", starting)
            count_line += 1
        f_out.close()

    #merge(fn_out, fn_merge)

def merge(fn_output, fn_merge):
    het_list = {}
    lines = []
    recent_site = -1
    num_file = 0
    for name in fn_output:
        out_file = open(name, 'r')
        for line in out_file: 
            if 'HET SITE' in line:
                recent_site = int(line.strip().split()[2])
                if num_file == 0:
                    het_list[recent_site] = []
            else:
                het_list[recent_site].append(line)
        num_file += 1
    f_out = open(fn_merge, 'w')
    for item in list(het_list.keys()):
        toWrite = 'HET SITE ' + str(item) + "\n"
        f_out.write(toWrite)
        for element in het_list[item]:
            f_out.write(element)
    f_out.close()
    print(list(het_list.keys()))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-s', '--sam', action='store', dest='sam_list',
                                    type=str, nargs='*', default=['item1', 'item2', 'item3'], help='list of sam file names')
    parser.add_argument('-f', '--fil', help='file with HET sites')
    #parser.add_argument('-o', '--out', help='output file')
    parser.add_argument('-o', '--out', action='store', dest='out_list', type=str, nargs='*', default=['item1', 'item2', 'item3'], help='list of output file names')
    #parser.add_argument('-l', '--l', help = 'merged output file name')
    args = parser.parse_args()
    
    fn_sam = args.sam_list
    fn_fil = args.fil
    fn_out = args.out_list
    #fn_merge = args.l

    print("fn_sam: ", fn_sam)
    print("fn_fil: ", fn_fil)
    print("fn_out: ", fn_out)
    #print("fn_merge: ", fn_merge)
    main(fn_sam, fn_fil, fn_out)#, fn_merge)
