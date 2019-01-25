'''
Compares two sam files and looks into the difference.
'''
import argparse
from analyze_sam import SamInfo, parse_line
from analyze_diploid_indels import build_index, check_var_in_region
from build_erg import read_var

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-rc', '--ref_correct',
        help='correct ref sam file'
    )
    parser.add_argument(
        '-ric', '--ref_incorrect',
        help='incorrect ref sam file'
    )
    parser.add_argument(
        '-pc', '--per_correct',
        help='correct per sam file'
    )
    parser.add_argument(
        '-pic', '--per_incorrect',
        help='incorrect per sam file'
    )
    parser.add_argument(
        '-v', '--var',
        help='the file specifying variants'
    )
    parser.add_argument(
        '--read_len', type=int,
        default=100,
        help='read length [100]'
    )
    args = parser.parse_args()
    global MAIN_HAP, ALT_HAP, MAIN_STRAND, ALT_STRAND, MAIN_CHRM, ALT_CHRM, REF_CHRM
    MAIN_HAP = 'hapA'
    ALT_HAP = 'hapB'
    MAIN_STRAND = 'A'
    ALT_STRAND = 'B'
    MAIN_CHRM = '9A'
    ALT_CHRM = '9B'
    REF_CHRM = '9'
    global READ_LEN
    READ_LEN = args.read_len
    return args

def check_var(dict, names, ref_index, main_index, alt_index, per):
    list_names = list(names)
    list_num_var = []
    for ll in list_names:
        var = dict[ll]
        if per:
            num_vars = check_var_in_region(var, main_index=main_index, alt_index=alt_index, MAIN_CHRM=MAIN_CHRM, ALT_CHRM=ALT_CHRM, READ_LEN=READ_LEN)
        elif per == 0:
            num_vars = check_var_in_region(var, main_index=ref_index, alt_index={}, MAIN_CHRM=REF_CHRM, ALT_CHRM='', READ_LEN=READ_LEN)
        else:
            print ('Error: unspecified per parameter', per)
        list_num_var.append(num_vars)
    sum_num_var = float(sum(list_num_var))
    len_list = len(list_num_var)
    avg_num_var = sum_num_var / len_list
    max_num_var = max(list_num_var)
    print ('  Avg num of vars =', avg_num_var)
    print ('  Max num of vars =', max_num_var)
    print ()

def compare_sam(args):
    ref_c_fn = args.ref_correct
    ref_ic_fn = args.ref_incorrect
    per_c_fn = args.per_correct
    per_ic_fn = args.per_incorrect
    var_fn = args.var

    ref_c_f = open(ref_c_fn, 'r')
    ref_ic_f = open(ref_ic_fn, 'r')
    per_c_f = open(per_c_fn, 'r')
    per_ic_f = open(per_ic_fn, 'r')
    
    var_list = read_var(var_fn, remove_conflict=True, remove_coexist=True)
    main_index, alt_index = build_index(var_list, per=2, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    ref_index, _ = build_index(var_list, per=0, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)

    dict_per_c = {}
    dict_per_ic = {}
    dict_per_u = {}
    for line in per_c_f:
        name, info = parse_line(line, by_score=0)
        if info.is_unaligned():
            dict_per_u[name] = info
        else:
            dict_per_c[name] = info
    for line in per_ic_f:
        name, info = parse_line(line, by_score=0)
        if info.is_unaligned():
            dict_per_u[name] = info
        else:
            dict_per_ic[name] = info
    dict_ref_c = {}
    dict_ref_ic = {}
    dict_ref_u = {}
    for line in ref_c_f:
        name, info = parse_line(line, by_score=0)
        if info.is_unaligned():
            dict_ref_u[name] = info
        else:
            dict_ref_c[name] = info
    for line in ref_ic_f:
        name, info = parse_line(line, by_score=0)
        if info.is_unaligned():
            dict_ref_u[name] = info
        else:
            dict_ref_ic[name] = info
    # item-based intersection:
    # intersection = dict(dict1.items() & dict2.items())
    # key-based intersection:
    # dict1.keys() & dict2.keys()
    dict_pc_rc = dict_per_c.keys() & dict_ref_c.keys()
    dict_pc_ric = dict_per_c.keys() & dict_ref_ic.keys()
    dict_pc_ru = dict_per_c.keys() & dict_ref_u.keys()
    dict_pu_rc = dict_per_u.keys() & dict_ref_c.keys()
    dict_pu_ru = dict_per_u.keys() & dict_ref_u.keys()
    dict_pu_ric = dict_per_u.keys() & dict_ref_ic.keys()
    dict_pic_rc = dict_per_ic.keys() & dict_ref_c.keys()
    dict_pic_ru = dict_per_ic.keys() & dict_ref_u.keys()
    dict_pic_ric = dict_per_ic.keys() & dict_ref_ic.keys()

    num_total = len(dict_pc_rc) + len(dict_pc_ric) + len(dict_pc_ru) + len(dict_pu_rc) + len(dict_pu_ru) + len(dict_pu_ric) + len(dict_pic_rc) + len(dict_pic_ru) + len(dict_pic_ric)

    print ('Num of ref-c and per-c =', len(dict_pc_rc))
    print ('  Looking dict_pc_rc in ref...')
    check_var(dict_ref_c, dict_pc_rc, ref_index, main_index, alt_index, per=0)
    print ('  Looking dict_pc_rc in per...')
    check_var(dict_per_c, dict_pc_rc, ref_index, main_index, alt_index, per=2)

    print ('Num of ref-c and per-ic =', len(dict_pic_rc))
    print ('  Looking dict_pic_rc in ref...')
    check_var(dict_ref_c, dict_pic_rc, ref_index, main_index, alt_index, per=0)
    print ('  Looking dict_pic_rc in per...')
    check_var(dict_per_ic, dict_pic_rc, ref_index, main_index, alt_index, per=2)
    
    print ('Num of ref-c and per-u =', len(dict_pu_rc))
    print ('  Looking dict_pu_rc in ref...')
    check_var(dict_ref_c, dict_pu_rc, ref_index, main_index, alt_index, per=0)

    print ('Num of ref-ic and per-c =', len(dict_pc_ric))
    print ('  Looking dict_pc_ric in ref...')
    check_var(dict_ref_ic, dict_pc_ric, ref_index, main_index, alt_index, per=0)
    print ('  Looking dict_pc_ric in per...')
    check_var(dict_per_c, dict_pc_ric, ref_index, main_index, alt_index, per=2)

    print ('Num of ref-ic and per-ic =', len(dict_pic_ric))
    print ('  Looking dict_pic_ric in ref...')
    check_var(dict_ref_ic, dict_pic_ric, ref_index, main_index, alt_index, per=0)
    print ('  Looking dict_pic_ric in per...')
    check_var(dict_per_ic, dict_pic_ric, ref_index, main_index, alt_index, per=2)

    print ('Num of ref-ic and per-u =', len(dict_pu_ric))
    print ('  Looking dict_pu_ric in ref...')
    check_var(dict_ref_ic, dict_pu_ric, ref_index, main_index, alt_index, per=0)

    print ('Num of ref-u and per-c =', len(dict_pc_ru))
    print ('  Looking dict_pc_ru in per...')
    check_var(dict_per_c, dict_pc_ru, ref_index, main_index, alt_index, per=2)

    print ('Num of ref-u and per-ic =', len(dict_pic_ru))
    print ('  Looking dict_pic_ru in per...')
    check_var(dict_per_ic, dict_pic_ru, ref_index, main_index, alt_index, per=2)

    print ('Num of ref-u and per-u =', len(dict_pu_ru))

    print ('\nNum total =', num_total)

if __name__ == '__main__':
    args = parse_args()
    compare_sam(args)