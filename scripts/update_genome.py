'''
Parse a vcf file and lists all the snvs and indels for a sample
'''

import sys
import argparse

# Update a genome with a set of SNPs and write to a new file


def get_mutation_type(info):
    '''
        Return the value of the VT attribute from the info field
    '''

    attrs = info.split(';')
    for a in attrs:
        if a[:3] == 'VT=':
            return a[3:]


def update_genome(indiv, seq, label, vcf, chrom, out_prefix, indels, var_only): 
    '''
    ##fileformat=VCFv4.1
    '''
    hapA = list(seq[:])
    hapB = list(seq[:])
    if var_only == 0:
        if indiv != None:
            fA = open(out_prefix + '_hapA.fa', 'w')
            split_label = label.split()
            label_A = split_label[0] + 'A ' + ' '.join(split_label[1:]) + '\n'
            fA.write(label_A)
            fB = open(out_prefix + '_hapB.fa', 'w')
            label_B = split_label[0] + 'B ' + ' '.join(split_label[1:]) + '\n'
            fB.write(label_B)
        else:
            fA = open(out_prefix + '.fa', 'w')
            fA.write(label)

    f_var = open(out_prefix + '.var', 'w')

    '''
    Format of a .var file:
    hap(A/B) chrm var_type ref_pos hap_pos offset
    '''
    if vcf != None:
        f = open(vcf, 'r')
    else:
        f = sys.stdin
    labels = None
    line_id = 0
    offsetA = 0
    offsetB = 0
    headA = 0
    headB = 0
    for line in f:
        # Skip header lines
        if line[0] == '#' and line[1] == '#':
            continue

        #: if "indiv" is set, select corresponding columns
        if not labels and indiv != None:
            labels = line.rstrip().split('\t')
            # print ('labels', labels)
            col = None
            for i in range(9, len(labels)):
                if labels[i] == indiv:
                    col = i
            if not col:
                print('Error! Couldn\'t find individual %s in VCF' % indiv)
                exit()
        else:
            row = line.rstrip().split('\t')
            type = get_mutation_type(row[7])

            # if type == 'SNP' or (indels and type == 'INDEL'):
            #: supports tri-allelic
            if type == 'SNP' or (indels and type in ['INDEL', 'SNP,INDEL']):
                # chrom = row[0]
                if row[0] != chrom:
                    continue
                loc = int(row[1])

                orig = row[3]
                alts = row[4].split(',')

                if indiv != None:
                    alleleA = int(row[col][0])
                    alleleB = int(row[col][2])
                else:
                    alleleA = 1
                    alleleB = 0
                    #: ignore multiallelic loci
                    # if len(alts) > 1:
                    #     continue
                    
                if alleleA > 0:
                    if len(orig) != len(alts[alleleA-1]):
                        type = 'INDEL'
                    else:
                        type = 'SNP'
                    if indels:
                        #: ignores conflicts
                        if loc >= headA:
                            new_offsetA = add_alt(hapA, loc-1, orig, alts[alleleA-1], offsetA)
                    else:
                        new_offsetA = 0
                        hapA[loc+offsetA-1] = alts[alleleA-1]
                    
                    if indels and (loc < headA):
                        pass
                    else:
                        f_var.write(
                            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                            ('A', chrom, type, str(loc), str(loc+offsetA), orig, alts[alleleA-1], str(new_offsetA), str(offsetB) )
                        )
                        offsetA = new_offsetA
                        headA = loc + len(orig)
                
                if alleleB > 0 and indiv != None:
                    if len(orig) != len(alts[alleleB-1]):
                        type = 'INDEL'
                    else:
                        type = 'SNP'
                    if indels:
                        #: ignores conflicts
                        if loc >= headB:
                            new_offsetB = add_alt(hapB, loc-1, orig, alts[alleleB-1], offsetB)
                    else:
                        new_offsetB = 0
                        hapB[loc+offsetB-1] = alts[alleleB-1]
                    
                    if indels and (loc < headB):
                        pass
                    else:
                        f_var.write(
                            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                            ('B', chrom, type, str(loc), str(loc+offsetB), orig, alts[alleleB-1], str(new_offsetB), str(offsetA) )
                        )
                        offsetB = new_offsetB
                        headB = loc + len(orig)

                line_id += 1

    f_var.close()
    
    if var_only == 0:
        for i in range(0, len(hapA), 60):
            fA.write(''.join(hapA[i:i+60])  + '\n')
        fA.close()

        if indiv != None:
            for i in range(0, len(hapB), 60):
                fB.write(''.join(hapB[i:i+60])  + '\n') 
            fB.close()
    f.close()

def add_alt(genome, loc, orig, alt, offset):
    '''
    loc here is the index for str
    i.e., loc = vcf_loc -1
    '''
    loc += offset

    if len(orig) == 1 and len(alt) == 1:
        # SNP
        # if genome[loc] != orig:
        #    print (loc, genome[loc], alt)
        genome[loc] = alt
    elif len(orig) > len(alt):
        # Deletion
        for i in range(len(alt)):
            genome[loc+i] = alt[i]
        del genome[loc+len(alt):loc+len(orig)]
        offset -= (len(orig) - len(alt))
    elif len(alt) > len(orig):
        # Insertion
        for i in range(len(orig)):
            genome[loc+i] = alt[i]
        genome[loc+len(orig):loc+len(orig)] = list(alt[len(orig):])
        offset += len(alt) - len(orig)
    else:
        # len(orig)=len(alt)>1 : Weird combination of SNP/In/Del
        for i in range(len(alt)):
            genome[loc+i] = alt[i]

    return offset

def read_chrom(ref, chrom):
    with open(ref, 'r') as f_in:
        label = None
        seq = ''
        for line in f_in:
            if line[0] == '>':
                if label:
                    return label, seq
                curr_chrom = line[1:].split(' ')[0]
                curr_chrom = curr_chrom.rstrip()
                if (curr_chrom == chrom) or (curr_chrom == 'chr' + chrom):
                    label = line
            elif label:
                seq += line.rstrip()

        return label, seq

if __name__ == '__main__':

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-r', '--ref', type=str, required=True, help='Path to fasta file containing reference genome'
    )
    parser.add_argument(
        '-v', '--vcf', type=str, help="Path to VCF file containing mutation information"
    )
    parser.add_argument(
        '-c', '--chrom', type=str, required=True, help="Chromosome to process"
    )
    parser.add_argument(
        '-op', '--out-prefix', type=str, required=True, help="Path to output prefix"
    )
    parser.add_argument(
        '-s', '--name', type=str, help="Name of individual in VCF to process; leave blank to allow all variants [None]"
    )
    parser.add_argument(
        '-i', '--include-indels', type=int, default=0, help="Set 1 to extract both SNPs and INDELs [0]."
    )
    parser.add_argument(
        '--var-only', type=int, default=0, help="Set 1 to report .var file only (no .fa) [0]."
    )


    args = parser.parse_args(sys.argv[1:])

    if args.name == None:
        print ('Note: no individual specified, all variants in chrom %s are included' % args.chrom)

    label, genome = read_chrom(args.ref, args.chrom)
    # update_genome(args.name, genome, label, args.vcf, args.out_prefix, args.include_indels)
    update_genome(
        indiv=args.name,
        seq=genome,
        label=label,
        vcf=args.vcf,
        chrom=args.chrom,
        out_prefix=args.out_prefix,
        indels=args.include_indels,
        var_only=args.var_only
    )
