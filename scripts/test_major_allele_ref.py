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
        help='sampled check points'
    )
    args = parser.parse_args()
    return args

def read_genome(fn):
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
            name = line.split()[0][1:]
            continue
        seq += line[: line.find('\\')]
    genome[name] = seq
    return genome

def test_major_allele_ref(fn_ma_ref, fn_ref, fn_check_points):
    ref = read_genome(fn_ref)
    ma_ref = read_genome(fn_ma_ref)
    f_cp = open(fn_check_points, 'r')
    for line in f_cp:
        line = line.split()
        [chrom, vid, pos, ref, alt, af] = line[:]
        pos = int(pos)
        if ref[chrom][pos] == ref and ma_ref[chrom][pos] == alt:
            print ('pass')
        else:
            input (line)
        
    # it = [*ref]
    # for i in it:
    #     print (len(ref[i]))
    # print (ref)

if __name__ == '__main__':
    args = parse_args()
    fn_ma_ref = args.ma_ref
    fn_ref = args.ref
    fn_check_points = args.check_points
    test_major_allele_ref(fn_ma_ref, fn_ref, fn_check_points)
