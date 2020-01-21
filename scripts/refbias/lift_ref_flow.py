from sys import exit
import argparse
import pandas as pd
#from os import path

# def main(fn_vcf, list_fn_sam, fn_fasta, list_fn_name, fn_out):
def main(fn_vcf, list_fn_sam, fn_fasta, fn_out):
    import liftover_sam
    ref_bi_list = []
    counter = 1

    for i in range(len(list_fn_sam)):
        print("i:", i)
        vcf = fn_vcf
        sam = list_fn_sam[i]
        fasta = fn_fasta
        
        print("vcf: ", vcf)
        print("sam: ", sam)
        print("fasta: ", fasta)
        
        output = list_fn_sam[i] + '_refbias.txt'
        liftover_sam.main(vcf, sam, fasta, output)
        #cigar_whole_genome_sam.main(vcf, sam, fasta, output)
        ref_bi_list.append(output)
        counter += 1
    print("ref_bi_list: ", ref_bi_list)
    merge(ref_bi_list, fn_out)
    
def merge(ref_bi_list, fn_out):
    list_df = []
    #: read TSVs from sub files
    for i, bias_file in enumerate(ref_bi_list):
        sub_df = pd.read_csv(bias_file, sep = '\t', comment = '#')
        list_df.append(sub_df)
        if i > 0:
            #: check columns are matched
            assert ((sub_df.columns) == (list_df[0].columns)).all()
            #: check shapes are matched
            assert sub_df.shape == list_df[0].shape
            #: check HET sites are matched
            assert (sub_df[['CHR', 'HET_SITE']] == list_df[0][['CHR', 'HET_SITE']]).all().all()

    df_out = pd.DataFrame(columns = list_df[0].columns)
    df_out[['CHR', 'HET_SITE']] = list_df[0][['CHR', 'HET_SITE']]
    df_out[['REF_COUNT', 'ALT_COUNT', 'GAP_COUNT', 'OTHER_COUNT', 'NUM_READS', 'SUM_MAPQ']] = list_df[0][['REF_COUNT', 'ALT_COUNT', 'GAP_COUNT', 'OTHER_COUNT', 'NUM_READS', 'SUM_MAPQ']]
    for i, sub_df in enumerate(list_df):
        if i > 0:
            df_out[['REF_COUNT', 'ALT_COUNT', 'GAP_COUNT', 'OTHER_COUNT', 'NUM_READS', 'SUM_MAPQ']] += sub_df[['REF_COUNT', 'ALT_COUNT', 'GAP_COUNT', 'OTHER_COUNT', 'NUM_READS', 'SUM_MAPQ']]
    df_out['REFERENCE_BIAS'] = df_out['REF_COUNT'] / (df_out['REF_COUNT'] + df_out['ALT_COUNT'])
    # df_out['REFERENCE_BIAS'] = df_out['REF_COUNT'] / df_out['NUM_READS']
    df_out.to_csv(fn_out, sep = '\t', index = None, float_format = '%.4f')
    return 

def get_paths_from_list(fn_list):
    l = []
    with open(fn_list, 'r') as f:
        for line in f:
            l.append(line.rstrip())
    return l

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', help='vcf file')
    parser.add_argument('-s', '--sam', help='list of sam file(s)')
    # parser.add_argument('-n', '--name', help='list of partial output name (s)')
    parser.add_argument('-f', '--fasta', help='reference fasta file')
    parser.add_argument('-o', '--out', help='output joined reference bias files')

    args = parser.parse_args()

    fn_vcf = args.vcf
    list_fn_sam = get_paths_from_list(args.sam)
    # list_fn_name = get_paths_from_list(args.name)
    fn_fasta = args.fasta
    fn_out = args.out

    print("vcf file: {}".format(fn_vcf))
    print("List of sam files: {}".format(list_fn_sam))
    print("fasta file: {}".format(fn_fasta))
    print("output: ", fn_out)
    main(
        fn_vcf=fn_vcf,
        list_fn_sam=list_fn_sam,
        fn_fasta=fn_fasta,
        # list_fn_name=list_fn_name,
        fn_out=fn_out
    )
