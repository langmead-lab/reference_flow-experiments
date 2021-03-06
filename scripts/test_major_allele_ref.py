'''
Test the accuracy of a major allele reference
'''
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m', '--ma_ref',
        help='the major allele ref file (fasta) to be tested'
    )
    parser.add_argument(
        '-r', '--ref',
        help='standard ref seq (fasta)'
    )
    parser.add_argument(
        '-c', '--check_points',
        default=None,
        help='sampled check points'
    )
    args = parser.parse_args()
    return args

def read_genome(fn, gatk=False):
    '''
    Reads genome and supports more than one chromosomes

    Inputs
        fn: fasta file name
        gatk: set True to use the second item of the header as chromosome name
    Output:
        genome in a list
    '''
    genome = {}
    f = open(fn, 'r')
    #: seq[0] is empty to fit vcf coordinate (1-based)
    seq = '^'
    name = None
    for line in f:
        if line.startswith('>'):
            if len(seq) > 1:
                genome[name] = seq
                seq = '^'
            #: update name
            if gatk:
                name = line.split()[1]
                name = name[: name.find(':')]
            else:
                name = line.split()[0][1:]
            continue
        seq += line[: line.find('\\')]
    genome[name] = seq
    return genome

def test_major_allele_ref(fn_ma_ref, fn_ref, fn_check_points):
    genome_ref = read_genome(fn_ref)
    genome_ma_ref = read_genome(fn_ma_ref)
    #genome_ma_ref = read_genome(fn_ma_ref, gatk=True)
    check_length = True
    '''
    #: the code below works if chromosome names are '1', '2', ...
    #: but doesn't work for 'chr1', 'chr2', ...
    for i in range(1,23):
        i = str(i)
        if len(genome_ref[i]) != len(genome_ma_ref[i]):
            check_length = False
            print ('Error: length mismatches')
            print ('ref', i, len(genome_ref[i]))
            print ('ma ', i, len(genome_ma_ref[i]))
    if check_length:
        print ('Length check: pass')
    '''

    #: supports stdin for pipelined commands
    if fn_check_points == None:
        f_cp = sys.stdin
    else:
        f_cp = open(fn_check_points, 'r')

    num_pass = 0
    num_fail = 0
    for line in f_cp:
        line = line.split()
        [chrom, vid, pos, ref, alt, af] = line[:]
        pos = int(pos)
        #: case-insentivie comparison, all chars are casted to uppercase
        genome_allele = genome_ref[chrom][pos: pos+len(ref)].upper()
        ma_genome_allele = genome_ma_ref[chrom][pos: pos+len(alt)].upper()
        vcf_ref_allele = ref.upper()
        vcf_alt_allele = alt.upper()
        if genome_allele == vcf_ref_allele and ma_genome_allele == vcf_alt_allele:
        #if (genome_ref[chrom][pos: pos+len(ref)].upper() == ref.upper()) and (genome_ma_ref[chrom][pos: pos+len(alt)].upper() == alt.upper()):
            num_pass += 1
            #print ('pass')
        else:
            num_fail += 1
            if __debug__:
                print ('false')
                # if genome_ref[chrom][pos : pos+len(ref)].upper() != ref.upper():
                if genome_allele != vcf_ref_allele:
                    print ('error in ref')
                    print ('ref:', genome_allele)
                    print ('vcf:', vcf_ref_allele)
                    # print ('ref:', genome_ref[chrom][pos : pos+len(ref)])
                    # print ('vcf:', ref)
                #if genome_ma_ref[chrom][pos : pos+len(alt)].upper() != alt.upper():
                if ma_genome_allele != vcf_alt_allele:
                    print ('error in ma_ref')
                    print ('ma:', ma_genome_allele)
                    print ('vcf:', vcf_alt_allele)
                    # print ('ma:', genome_ma_ref[chrom][pos : pos+len(alt)])
                    # print ('vcf:', alt)
                print (line)
                print ()

    print ('check', fn_check_points)
    print ('  pass:', num_pass)
    print ('  fail:', num_fail)

if __name__ == '__main__':
    args = parse_args()
    fn_ma_ref = args.ma_ref
    fn_ref = args.ref
    fn_check_points = args.check_points
    test_major_allele_ref(fn_ma_ref, fn_ref, fn_check_points)
