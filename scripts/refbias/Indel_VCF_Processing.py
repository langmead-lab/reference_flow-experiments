import argparse

def main(fn_liftover, fn_normal, fn_output):
    file = open(fn_liftover, 'r')#open('offset.vcf', 'r')#open("21_major.recode.vcf")

    maj_pos = []
    offset = []
    ref_alleles = []
    alt_alleles = []
    for line in file:
        if '#' not in line:
            spl = line.split()
            pos = int(spl[1])
            ref = spl[3]
            alt = spl[4]
            if len(ref) > len(alt) or len(alt) > len(ref):
                offset.append(len(alt) - len(ref))
                maj_pos.append(pos)
                ref_alleles.append(ref)
                alt_alleles.append(alt)

    #print("maj_pos: ", maj_pos)
    #print("offset: ", offset)
    #print("ref_alleles: ", ref_alleles)
    #print("alt_alleles: ", alt_alleles)

    sites_to_ignore = []
    file = open(fn_normal, 'r')#open('test_ref_het.vcf', 'r')#open("21_thousand_out.vcf")
    out = open(fn_output, 'w')#open('test_het_mejor.vcf', 'w')#open("major_21_thousand.vcf", "w")
    for line in file:
        if "#" not in line:
            spl = line.split()
            location = int(spl[1])
            #print("location: ", location)
            dist = 0
            if location in maj_pos:
                #print("here")
                index = maj_pos.index(location)
                #print("index: ", index)
                #print("ref[index]: ", ref_alleles[index])
                #print("alt[index]: ", alt_alleles[index])
                if len(ref_alleles[index]) > len(alt_alleles[index]):#if it is a deletion for maj:
                    for i in range(0, len(ref)):
                        sites_to_ignore.append(location + i)
                    continue
            for i in range(len(maj_pos)):
                pos = maj_pos[i]
                if pos < location:
                    dist += offset[i]
            #print(location, "-->", location + dist)
            new_line = [str(location + dist) if i == 1 else spl[i] for i in range(len(spl))]
            string = '\t'.join(new_line)
            #print(string)
            out.write(string)
            out.write("\n")
        #else:
        #    out.write(line)
    print("het sites ignored: ", sites_to_ignore)
    print("# het sites ignored: ", len(sites_to_ignore))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--liftover', help='vcf file from grch37 to maj')
    parser.add_argument('-n', '--normal', help='vcf with normal het site file')
    parser.add_argument('-o', '--out', help='output vcf file with lifted het sites')
    args = parser.parse_args()
    fn_liftover = args.liftover
    fn_normal = args.normal
    fn_output = args.out
    print('vcf', fn_liftover)
    print('sam', fn_normal)
    print('output', fn_output)
    main(fn_liftover, fn_normal, fn_output)
