#! /usr/bin/env python2.7
'''
Last update: 2018/12/11 by Nae-Chyun Chen

(Original script is from FORGe-experiments by Jacob Pritt)
Parse a vcf file and lists all the snps and
indels for a sample
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


def update_genome(indiv, seq, label, vcf, out_prefix, indels=None): 
    hapA = list(seq[:])
    hapB = list(seq[:])
    fA = open(out_prefix + '_hapA.fa', 'w')
    fA.write(label)
    fB = open(out_prefix + '_hapB.fa', 'w')
    fB.write(label)

    f_var = open(out_prefix + '.var', 'w')

    '''
    Format of a .var file:
    hap(A/B) chrm var_type ref_pos hap_pos offset
    '''

    with open(vcf, 'r') as f:
        labels = None

        line_id = 0

        offset = 0
        offsetA = 0
        offsetB = 0
        for line in f:
            # Skip header lines
            if line[0] == '#' and line[1] == '#':
                continue

            if not labels:
                labels = line.rstrip().split('\t')

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
                if type == 'SNP' or (indels and type == 'INDEL'):
                    chrom = row[0]
                    name = row[2]
                    loc = int(row[1])

                    orig = row[3]
                    alts = row[4].split(',')

                    alleleA = int(row[col][0])
                    alleleB = int(row[col][2])

                    if alleleA > 0:
                        #f_var.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('A', chrom, type, str(loc), str(loc+offsetA), orig, alts[alleleA-1], str(offsetA) ))
                        f_var.write(
                            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                            ('A', chrom, type, str(loc), str(loc+offsetA), orig, alts[alleleA-1], str(offsetA), str(offsetB) )
                        )
                        if indels:
                            offsetA = add_alt(hapA, loc-1, orig, alts[alleleA-1], offsetA)
                        else:
                            hapA[loc+offsetA-1] = alts[alleleA-1]
                    if alleleB > 0:
                        #f_var.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('B', chrom, type, str(loc), str(loc+offsetB), orig, alts[alleleB-1], str(offsetB) ))
                        f_var.write(
                            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                            ('B', chrom, type, str(loc), str(loc+offsetB), orig, alts[alleleB-1], str(offsetB), str(offsetA) )
                        )
                        if indels:
                            offsetB = add_alt(hapB, loc-1, orig, alts[alleleB-1], offsetB)
                        else:
                            hapB[loc+offsetB-1] = alts[alleleB-1]

                    line_id += 1

    f_var.close()
    for i in range(0, len(hapA), 60):
        fA.write(''.join(hapA[i:i+60])  + '\n')
    for i in range(0, len(hapB), 60):
        fB.write(''.join(hapB[i:i+60])  + '\n') 

    fA.close()
    fB.close()

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
                if curr_chrom == chrom:
                    label = line
            elif label:
                seq += line.rstrip()

        return label, seq

if __name__ == '__main__':

    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ref', type=str, required=True, help='Path to fasta file containing reference genome')
    parser.add_argument("--vcf", type=str, required=True, help="Path to VCF file containing mutation information")
    parser.add_argument("--chrom", type=str, required=True, help="Chromosome to process")
    parser.add_argument("--out-prefix", type=str, required=True, help="Path to output prefix")
    parser.add_argument("--name", type=str, help="Name of individual in VCF to process")
    parser.add_argument("--include-indels", type=str, help="If present, extract both SNPs and INDELs.")

    args = parser.parse_args(sys.argv[1:])

    label,genome = read_chrom(args.ref, args.chrom)
    update_genome(args.name, genome, label, args.vcf, args.out_prefix, args.include_indels)

