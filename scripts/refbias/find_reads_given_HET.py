import re
import os.path
import argparse
import pandas as pd
import random
from utils import get_paths_from_list
from count_variants import count_var

def main(fn_sam, fn_vcf, fn_het, fn_out, fn_merge, list_range, sample_rate):
    df_het = pd.read_csv(fn_het, sep='\t')
    filter = False
    for r in list_range:
        print (r)
        filter = filter | df_het['REFERENCE BIAS'].between(r[0], r[1])
    df_pass = df_het[filter]
    list_het_raw = df_pass['HET SITE']
    list_het = []
    for het_raw in list_het_raw:
        het_raw = het_raw.rstrip()
        het = [int(i) for i in het_raw.strip('[]').split(',')]
        for h in het:
            assert h == het[0]
        list_het.append(het[0])
    df_pass['HET SITE'] = list_het

    # f_vcf = open(fn_vcf, 'r')
    dict_het_gt = {}
    with open(fn_vcf, 'r') as f_vcf:
        for line in f_vcf:
            if line.startswith('#'):
                continue
            line = line.split()
            if int(line[1]) in list_het:
                dict_het_gt[int(line[1])] = line[-1]
    assert len(list_het) == len(dict_het_gt.keys())

    #: skips this part if all the output files exist
    all_out_exist = True
    for fn in fn_out:
        if not os.path.exists(fn):
            all_out_exist = False
    if all_out_exist:
        merge(fn_out, fn_merge, fn_vcf, dict_het_gt, df_pass)
        return

    set_het_sampled = set()
    for i_sam in range(len(fn_sam)):
        element = fn_sam[i_sam]
        #if os.path.exists(fn_out[i_sam]):
        #    continue
        f_out = open(fn_out[i_sam], 'w')
        print(element)
        sam_file = open(element, 'r')
        sam_name = []
        sam_reads = []
        sam_offset = []
        sam_mapQ = []
        sam_line = []
        for line in sam_file:
            if not line.startswith('@'):
                spl = line.split()
                start_pos = int(spl[3])
                # start_pos = int(spl[3])-1
                tag =int(spl[1])
                if (tag & 4):
                    continue
                # chr = int(spl[2])
                cigar = spl[5]
                #start_pos = int(spl[3]) - 1 #T
                sequence = spl[9]
                mod_sequence = ''
                if tag & 4:
                    continue  
            
                if not cigar == (str(len(sequence))+'M') and not (tag & 4):
                    #ref = reference[start_pos:start_pos + 102]
                    change = 0
                    start = start_pos
                    count_del = 0
                    count_ins = 0
                    for num1, idm in re.findall('(\d+)([IDMS])', cigar):
                        # print(start_pos)
                        if idm == 'M':
                            mod_sequence += sequence[change:change + int(num1)]
                        elif idm == 'D':
                            count_del += int(num1)
                            for mod_i in range(int(num1)):
                                mod_sequence += '-'
                        elif idm == 'I':
                            count_ins += int(num1)
                        elif idm == 'S':
                            count_ins += int(num1)
                        else:
                            print ('error: unexpected cigar letter', cigar)
                            exit ()

                        if idm != 'D':
                            change += int(num1)
                else:
                    mod_sequence = sequence

                sam_offset.append(start_pos)
                sam_reads.append(mod_sequence)
                sam_name.append(spl[0])
                sam_mapQ.append(spl[4])
                sam_line.append(line)
            #print("just finished a sam file")
        starting = 0

        have_started = False
        for het_pos in list_het:
            if sample_rate != 1:
                if i_sam == 0:
                    if random.random() > sample_rate:
                        continue
                    else:
                        set_het_sampled.add(het_pos)
                elif het_pos not in set_het_sampled:
                    continue

            toWrite = "HET SITE " + str(het_pos)
            f_out.write(toWrite)
            f_out.write("\n")
            for i in range(starting, len(sam_name)):
                align = sam_offset[i]
                ran = range(align, align + len(sam_reads[i]))
                if het_pos in ran:
                    if not have_started:
                        have_started = True
                        starting = i
                    toWrite = sam_line[i]
                    # toWrite = sam_name[i] + "\t" + str(align) + "\t" + sam_reads[i] + "\t" + sam_mapQ[i] + "\n"
                    f_out.write(toWrite)
                else:
                    if align > het_pos:
                        #print("i broke")
                        break
                #print("starting: ", starting)
        f_out.close()

    merge(fn_out, fn_merge, fn_vcf, dict_het_gt, df_pass)

def merge(fn_output, fn_merge, fn_vcf, dict_het_gt, df):
    het_list = {}
    lines = []
    recent_site = -1
    num_file = 0
    for name in fn_output:
        out_file = open(name, 'r')
        for line in out_file: 
            if 'HET SITE' in line:
                recent_site = int(line.strip().split()[2])
                if num_file == 0:
                    het_list[recent_site] = []
            else:
                het_list[recent_site].append(line)
        num_file += 1
    f_out = open(fn_merge, 'w')
    list_read_bias = []
    list_ref_read_bias = []
    list_avg_mapq = []
    for item in list(het_list.keys()):
        toWrite = 'HET SITE ' + str(item) + "\n"
        f_out.write(toWrite)
        count_hapA = 0
        count_hapB = 0
        sum_mapq = 0
        for element in het_list[item]:
            f_out.write(element)
            assert element.count('hapA') + element.count('hapB') == 1
            if element.count('hapA') > 0:
                count_hapA += 1
            else:
                count_hapB += 1
            sum_mapq += int(element.rstrip().split('\t')[4])

        ref_read_bias = calc_read_bias(het_list[item], dict_het_gt[item])
        list_ref_read_bias.append(ref_read_bias)
        read_bias = count_hapA / (count_hapA + count_hapB)
        list_read_bias.append(read_bias)
        list_avg_mapq.append(sum_mapq / len(het_list[item]))
    f_out.close()

    print ('{} HET sites to merge...'.format(len(het_list)))

    list_var, list_numA, list_numB = \
        count_var(fn_vcf, '', target_het = het_list)
    assert set(het_list) == set(list_var)
    list_num_var = [max(list_numA[i], list_numB[i]) for i in range(len(list_numA))]

    df = df.loc[df['HET SITE'].isin(list_var)]
    # df.assign(num_var = list_num_var, read_bias = list_read_bias, ref_read_bias = list_ref_read_bias, avg_mapq = list_avg_mapq)
    df['NUM VAR'] = list_num_var

    df['READ BIAS'] = list_read_bias
    df['REF READ BIAS'] = list_ref_read_bias
    df['AVG MAPQ'] = list_avg_mapq
    df.to_csv(fn_merge + '.tsv', sep='\t')

def calc_read_bias(list_reads, gt):
    gt = [int(i) for i in gt.split('|')]
    assert gt[0] + gt[1] > 0
    assert gt[0] * gt[1] == 0
    if gt[0] != 0:
        ref = 'hapB'
        alt = 'hapA'
    else:
        ref = 'hapA'
        alt = 'hapB'

    count_ref = 0
    count_alt = 0
    for r in list_reads:
        assert r.count(ref) + r.count(alt) == 1
        if r.count(ref) > 0:
            count_ref += 1
        else:
            count_alt += 1
    return count_ref / (count_ref + count_alt)

def range_to_lists(r):
    r = r.split(',')
    list_r = []
    for pairs in r:
        pairs = pairs.split('-')
        assert len(pairs) == 2
        list_r.append([float(p) for p in pairs])
    return list_r

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sam', help='list of sam file(s)')
    parser.add_argument('-v', '--vcf', help='vcf file specifying HETs')
    parser.add_argument('-f', '--het', help='file with HET sites')
    # parser.add_argument('-o', '--out', help='list of output file name(s)')
    parser.add_argument('-id', '--id', help='list of id(s)')
    parser.add_argument('-m', '--merge', help = 'merged output file name')
    parser.add_argument(
        '-r', '--range', default='0-1',
        help='ref bias range, values are separated by commas, \
            e.g. 0-0.2,0.8-1 removes HETs with bias from 0.2 to 0.8 [0-1]')
    #parser.add_argument('-v', '--val', help='ref bias range specified by user')
    parser.add_argument('--sample', type=float, default=1.0, help='sampling rate for HET sites [1.0]')
    args = parser.parse_args()
    
    list_fn_sam = get_paths_from_list(args.sam)
    fn_vcf = args.vcf
    fn_het = args.het
    fn_merge = args.merge
    list_fn_out = get_paths_from_list(
        args.id,
        prefix = fn_merge + '_'
        # suffix = '_' + args.range + '.reads'
        # auto_prefix = list_fn_sam
    )
    list_range = range_to_lists(args.range)
    sample_rate = args.sample

    random.seed(0)
    
    print("fn_sam:   ", list_fn_sam)
    print("fn_vcf:   ", fn_vcf)
    print("fn_het:   ", fn_het)
    print("fn_out:   ", list_fn_out)
    print("fn_merge: ", fn_merge)
    #print("var     : ", val)
    print("range    : ", list_range)
    print("sample   : ", sample_rate)

    # main(list_fn_sam, fn_het, list_fn_out, fn_merge, list_range)
    main(
        fn_sam = list_fn_sam,
        fn_vcf = fn_vcf,
        fn_het = fn_het,
        fn_out = list_fn_out,
        fn_merge = fn_merge,
        list_range = list_range,
        sample_rate = sample_rate
        )
