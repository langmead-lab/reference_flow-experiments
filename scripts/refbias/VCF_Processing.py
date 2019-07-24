import argparse
#import time #('time: ', 19.679501056671143)

def main(fn_vcf, fn_output):
    #t = time.time()
    file = open(fn_vcf, 'r')
    f = open(fn_output, 'a+')

    index = -1
    toAdd = []
    should = True
    for line in file:
        if should and line.startswith("##"):
            should = should
            #f.write(line)
        else:
            if line.startswith("#"):
                categories = line.split()
                for i in range(len(categories)):
                    if categories[i] == 'NA12878':
                        index = i
                toAdd = [categories[i] for i in range(8)]
                toAdd.append(categories[index])
            else:
                spl = line.split()
                if should and int(spl[0]) == 1:
                    f.write('\t'.join(toAdd))
                    should = False
                ref = spl[3]
                alt = spl[4]
                goAhead = True
                if len(ref) > 1 and ',' in ref:
                    for element in ref.split(","):
                        if len(element) > 1:
                            goAhead = False
                elif len(alt) > 1 and ',' in alt:
                    for element in alt.split(","):
                        if len(element) > 1:
                            goAhead = False
                elif len(alt) > 1 or len(ref) > 1:
                    goAhead = False

                if goAhead:
                    s = set(spl[index].split("|"))
                    if len(s) > 1 and '0' in s:
                        toAdd = [spl[i] for i in range(8)]
                        toAdd.append(spl[index])
                        f.write("\t".join(toAdd))
                        f.write("\n")
        #print("time: ", time.time() - t)
    f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file for chromosomes')
    parser.add_argument('-o', '--out', help='output vcf file with HET sites')
    args = parser.parse_args()
    fn_vcf = args.vcf
    fn_output = args.out
    print('vcf', fn_vcf)
    print('output', fn_output)
    main(fn_vcf, fn_output)
