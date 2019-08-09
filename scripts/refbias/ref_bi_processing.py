import argparse

def main(fn_ref, fn_high_thresh, fn_low_thresh, fn_output):

    THRESHOLD_one = int(fn_high_thresh)#15
    THRESHOLD_two = int(fn_low_thresh)#3

    #file = open('2_ref_bi_no_n_het_sites.txt', 'r')
    file = open(fn_ref, 'r')
    f_out = open(fn_output, 'w')
    lines = []
    count = 0
    hets = {}
    #print("THRESHOLD: ", THRESHOLD_one)
    for line in file:
        lines.append(line)
        if count > 0:
            spl = line.split()
            if int(spl[4]) >= THRESHOLD_one:
                hets[int(spl[0])] = spl[1]
        count += 1

    #print("hets: ", hets)

    sum_ref_bi = 0
    for ref_bi in hets.values():
        sum_ref_bi += float(ref_bi)
    f_out.write("Number of reads: ")
    f_out.write(str(len(list(hets.keys()))))
    f_out.write("\n")
    f_out.write("Average reference bias: ")
    f_out.write(str(sum_ref_bi/len(list(hets.keys()))))
    f_out.write("\n")
    #print("Number of reads: ", len(list(hets.keys())))
    #print("Average reference bias: ", sum_ref_bi/len(list(hets.keys())), "\n")

    count = 0
    start = -1
    end = -1
    index_start = -1
    index_end = -1
    locs = []
    final = []
    for line in lines:
        if count > 0:
            spl = line.split()
            if int(spl[4]) <= THRESHOLD_two and start == -1:
                start = int(spl[0])
                index_start = count

            if int(spl[4]) <= THRESHOLD_two and start != -1:
                end = int(spl[0])
                index_end = count

            if int(spl[4]) > THRESHOLD_two and start != -1:
                locs.append((start, end, index_start, index_end))
                start = -1
                end = -1
        count += 1

    for pair in locs:
        if pair[1] - pair[0] >= 2000:
            final.append(pair)

    #print("final: ", final)
    count = 1
    for pair in final:
        s = pair[2]
        e = pair[3]
        string = "REGION " + str(count) + ": " + str(pair[1] - pair[0])  + " bp long"
        f_out.write(string)
        f_out.write("\n")
        #print("REGION ", count, ": ", pair[1] - pair[0], "bp long")
        for i in range(s, e+1):
            f_out.write(lines[i])
            #print(lines[i])
            f_out.write("\n")
        #print("\n")
        count += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref', help='reference bias text file')
    parser.add_argument('-g', '--good', help='good coverage threshold')
    parser.add_argument('-l', '--low', help='low coverage threshold')
    parser.add_argument('-o', '--out', help='output file')
    args = parser.parse_args()
    fn_ref = args.ref
    fn_high_thresh = args.good
    fn_low_thresh = args.low
    fn_output = args.out
    print('ref', fn_ref)
    print('high', fn_high_thresh)
    print('low', fn_low_thresh)
    print('output', fn_output)
    main(fn_ref, fn_high_thresh, fn_low_thresh, fn_output)