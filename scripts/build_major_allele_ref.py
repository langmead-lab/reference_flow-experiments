import sys

vcf_fname = "/scratch/groups/blangme2/naechyun/relaxing/na12878/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
indiv = "NA12878"

def get_info(line):
    vcf_d = {}
    vcf_d['chrom'] = line[0]
    vcf_d['pos'] = line[1]
    vcf_d['rsid'] = line[2]
    #vcf_d['ref'] = line[3]
    vcf_d['alt'] = line[4]
    #vcf_d['qual'] = line[5]
    vcf_d['filt'] = line[6]
    vcf_d['info'] = line[7]
    vcf_d['af'] = 0
    vcf_d['type'] = ''
    tags = line[7].split(';')
    for f in tags:
        if f.startswith('VT=SNP'): # type: snp
            vcf_d['type'] = 'single'
        if f.startswith('VT=INDEL'): # type: indel
            vcf_d['type'] = 'indel'
        if f.startswith('AF='): # overall allele frequency
            vcf_d['af'] = f
    if vcf_d['type'] == '':
        sys.stderr.write("WARNING: No variant type found\n")
        str_tags = ','.join(tags)
        sys.stderr.write(str_tags)
    #input (vcf_d)
    return vcf_d

with open(vcf_fname, 'r') as vcf_f:
    sys.stderr.write("WARNING: Multiple ALTs are so far not dealt with!\n")
    for line in vcf_f:
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            line = line.split('\t')
            idx_format = 0
            for i, e in enumerate(line):
                if e == 'FORMAT':
                    idx_format = i
                if e == indiv:
                    index = i
                    break
            continue
        elif len(line) < 1:
            break
        line = line.split('\t')
        # "indiv" genotypes don't count
        line[index] = '-'
        #info = line[7]
        genotypes = line[idx_format + 1 : ]
        # remove '\' in the last column
        genotypes[-1] = genotypes[-1][: genotypes[-1].find('\\')]
        cnt_total = len(genotypes) * 2
        cnt_var = 0
        for e in genotypes:
            if e == '1|0' or e == '0|1':
                cnt_var += 1
            elif e == '1|1':
                cnt_var += 2
        alt_is_major = False
        if cnt_var > 0.5 * cnt_total:
            alt_is_major = True
            vcf_d = get_info(line)
            #print (vcf_d['af'])
            #print ("[%s] %d / %d (%f): %r" % \
                    (line[1], cnt_var, cnt_total, float(cnt_var) / cnt_total, alt_is_major))
            if vcf_d['type'] is 'single':
                print ("%s\t%s\t%s\t%s\t%s" % \
                    (vcf_d['rsid'], vcf_d['type'], vcf_d['chrom'], vcf_d['pos'], vcf_d['alt']))
            #input()
