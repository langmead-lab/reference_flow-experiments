'''
Updates a genome with a set of SNPs (and INDELS) and
write to a new file.
'''

import sys
import argparse
import random

def get_mutation_type(info):
    '''
    Returns the value of the VT attribute from the info field
    '''

    attrs = info.split(';')
    for a in attrs:
        if a[:3] == 'VT=':
            return a[3:]

def get_allele_freq(info, num_haps):
    '''
    Returns allele frequency for a variant.
    Not using the "AF" attribute because it's calculated
    based on the entire 1KG population.
    '''
    attrs = info.split(';')
    for a in attrs:
        if a[:3] == 'AC=':
            try:
                count = int(a[3:])
            #: when there are multiple alleles,
            #: use the highest frequency
            except:
                a = a[3:].split(',')
                inta = [int(i) for i in a]
                count = max(inta)
    return float(count) / num_haps

def update_genome(
    indiv,
    seq,
    label,
    vcf,
    chrom,
    out_prefix,
    indels,
    var_only,
    is_stochastic,
    block_size
):
    '''
    ##fileformat=VCFv4.1
    '''
    hapA = list(seq[:])
    hapB = list(seq[:])
    # if var_only == 0:
    if not var_only:
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
    f_vcf = open(out_prefix + '.vcf', 'w')

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
    num_haps = 0

    current_block_pos = 0
    rr = 0 # initial number for random.random()

    for line in f:
        #: Skip header lines
        if line[0] == '#' and line[1] == '#':
            f_vcf.write(line)
            continue
        if line[0] == '#':
            f_vcf.write(line)
            labels = line.rstrip().split('\t')
            num_haps = 2 * (len(labels) - 9)
            #: if "indiv" is set, select corresponding columns
            if indiv != None:
                col = None
                for i in range(9, len(labels)):
                    if labels[i] == indiv:
                        col = i
                if not col:
                    print('Error! Couldn\'t find individual %s in VCF' % indiv)
                    exit()
            continue
        row = line.rstrip().split('\t')
        type = get_mutation_type(row[7])

        if row[0] != chrom:
            continue
        loc = int(row[1])

        if is_stochastic:
            freq = get_allele_freq(row[7], num_haps)
            #: only updates the random number when exceeding current block
            if loc >= current_block_pos + block_size:
                # print ('--update block--')
                # print ('prev rr = {0}, block_pos = {1}'.format(rr, current_block_pos))
                rr = random.random()
                current_block_pos = int(loc / block_size) * block_size
                # print ('updt rr = {0}, block_pos = {1}'.format(rr, current_block_pos))

            if rr > freq:
                continue
            # print ('selected, rr = {}'.format(rr), row[:2], freq)

        #: supports tri-allelic
        if type == 'SNP' or (indels and type in ['INDEL', 'SNP,INDEL']):
            # if row[0] != chrom:
            #     continue

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
                flag_skip = False
                if indels:
                    #: ignores conflicts or overlapped variants
                    #: but accepts overlapped INS
                    if loc == headA-1 and (len(orig) < len(alts[alleleA-1])):
                        print ('Warning: overlapped INS at {} for hapA'.format(loc))
                        new_offsetA = add_alt(hapA, loc-1, orig, alts[alleleA-1], offsetA, True)
                    elif loc >= headA:
                        new_offsetA = add_alt(hapA, loc-1, orig, alts[alleleA-1], offsetA, False)
                    else:
                        flag_skip = True
                        print ('Warning: conflict at {} for hapA'.format(loc))
                else:
                    new_offsetA = 0
                    hapA[loc+offsetA-1] = alts[alleleA-1]
                
                if not flag_skip:
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
                flag_skip = False
                if indels:
                    #: ignores conflicts or overlapped variants
                    #: but accepts overlapped INS
                    if loc == headB-1 and (len(orig) < len(alts[alleleB-1])):
                        print ('Warning: overlapped INS at {} for hapB'.format(loc))
                        new_offsetB = add_alt(hapB, loc-1, orig, alts[alleleB-1], offsetB, True)
                    elif loc >= headB:
                        new_offsetB = add_alt(hapB, loc-1, orig, alts[alleleB-1], offsetB, False)
                    else:
                        flag_skip = True
                        print ('Warning: conflict at {} for hapB'.format(loc))
                else:
                    new_offsetB = 0
                    hapB[loc+offsetB-1] = alts[alleleB-1]
                
                if not flag_skip:
                    f_var.write(
                        '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                        ('B', chrom, type, str(loc), str(loc+offsetB), orig, alts[alleleB-1], str(new_offsetB), str(offsetA) )
                    )
                    offsetB = new_offsetB
                    headB = loc + len(orig)

            if (alleleA > 0) or \
                (alleleB > 0 and indiv != None):
                f_vcf.write(line)
            line_id += 1
    
    f_vcf.close()
    f_var.close()
    
    # if var_only == 0:
    if not var_only:
        for i in range(0, len(hapA), 60):
            fA.write(''.join(hapA[i:i+60])  + '\n')
        fA.close()

        if indiv != None:
            for i in range(0, len(hapB), 60):
                fB.write(''.join(hapB[i:i+60])  + '\n') 
            fB.close()
    f.close()

def add_alt(genome, loc, orig, alt, offset, overlap_ins):
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
        #: Insertion
        
        # if overlap_ins:
        #     for i in range(len(orig)):
        #         if genome[loc+i] != alt[i]:
        #             print ('Warning: genome ({0}) differs from ALT ({1}) at {2}'.format(genome[loc+i], alt[i], loc))
        
        #: don't replace if overlap_ins is True
        if not overlap_ins:
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
        if label == None:
            print ('Error: no matched chromosome. Label = {}'.format(chrom))
            exit ()
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
        '-i', '--include-indels', action='store_true', help="Set to extract both SNPs and INDELs"
    )
    # parser.add_argument(
    #     '-i', '--include-indels', type=int, default=0, help="Set 1 to extract both SNPs and INDELs [0]"
    # )
    # parser.add_argument(
    #     '-S', '--stochastic', type=int, default=0, help="Set 1 to enable stochastic flipping [0]"
    # )
    parser.add_argument(
        '-S', '--stochastic', action='store_true', help="Set to enable stochastic flipping"
    )
    parser.add_argument(
        '-rs', '--rand-seed', help="random seed for controlled randomness [None]"
    )
    # parser.add_argument(
    #     '--var-only', type=int, default=0, help="Set 1 to report .var file only (no .fa) [0]"
    # )
    parser.add_argument(
        '--var-only', action='store_true', help="Set to report .var file only (no .fa output)"
    )
    parser.add_argument(
        '-b', '--block-size', type=int, default=1, help="Size of block for stochastic update [1]"
    )


    args = parser.parse_args(sys.argv[1:])

    if args.name == None:
        print ('Note: no individual specified, all variants in chrom %s are included' % args.chrom)
    # if args.stochastic == 1:
    if args.stochastic:
        print ('Note: stochastic update is enabled')
        if args.rand_seed:
            random.seed(args.rand_seed)
            print ('Set random seed: {}'.format(args.rand_seed))

    label, genome = read_chrom(args.ref, args.chrom)
    update_genome(
        indiv = args.name,
        seq = genome,
        label = label,
        vcf = args.vcf,
        chrom = args.chrom,
        out_prefix = args.out_prefix,
        indels = args.include_indels,
        var_only = args.var_only,
        is_stochastic = args.stochastic,
        block_size = args.block_size
    )
