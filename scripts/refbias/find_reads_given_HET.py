import argparse
#import pandas as pd

def main(fn_sam, fn_fil, fn_out):

    f_out = open(fn_out, 'w')


    sam_file = open(fn_sam, 'r')
    sam_name = []
    sam_reads = []
    sam_offset = []

    for line in sam_file:
        if not line.startswith('@'):
            spl = line.split()
            # sam[int(spl[3])] = spl[9]
            sam_name.append(spl[0])
            sam_offset.append(int(spl[3]))
            sam_reads.append(spl[9])

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
                    toWrite = sam_name[i] + "\t" + str(align) + "\t" + sam_reads[i] + "\n"
                    f_out.write(toWrite)
                else:
                    if align > het_pos:
                        print("i broke")
                        break
                print("starting: ", starting)
        count_line += 1
    f_out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sam', help='sam file')
    parser.add_argument('-f', '--fil', help='file with HET sites')
    parser.add_argument('-o', '--out', help='output file')


    args = parser.parse_args()
    fn_sam = args.sam
    fn_fil = args.fil
    fn_out = args.out

    print("fn_sam: ", fn_sam)
    print("fn_fil: ", fn_fil)
    print("fn_out: ", fn_out)

    main(fn_sam, fn_fil, fn_out)
