'''
This script reads a var file and removes variant records 
which had occurred before.

If there's a conflict, the script exits and shows error report.
'''
import argparse
from build_erg import VarInfo


MAIN_STRAND = 'A'
ALT_STRAND = 'B'

def remove_redundant_vars(var_fn, out_fn):
    SHOW_REMOVED = False
    removed_count = 0
    var_f = open(var_fn, 'r')
    out_f = open(out_fn, 'w')
    del_pos = {MAIN_STRAND:[], ALT_STRAND:[]}
    del_allele = {MAIN_STRAND:[], ALT_STRAND:[]}
    var_list = {MAIN_STRAND:[], ALT_STRAND:[]}
    for line in var_f:
        v = VarInfo(line)
        if v.is_del():
            vd = del_pos[v.strand]
            if len(vd) > 0 and vd[len(vd) - 1] < v.ref_pos:
                del_pos[v.strand] = []
                del_allele[v.strand] = []
            for i in range(len(v.ref_allele)):
                del_pos[v.strand].append(v.ref_pos + i)
                del_allele[v.strand].append(v.ref_allele[i])
            var_list[v.strand].append(v)
            out_f.write(v.line + '\n')
            continue
        if v.ref_pos in del_pos[v.strand]:
            d = del_allele[v.strand][del_pos[v.strand].index(v.ref_pos)]
            if d != v.ref_allele:
                print ('Error: conflict at', v.ref_pos, d)
                input (v.line)
                exit()
            removed_count += 1
            if SHOW_REMOVED:
                print ('"%s"' % v.strand, v.ref_pos, d)
        else:
            var_list[v.strand].append(v)
            out_f.write(v.line + '\n')
    
    print (removed_count, 'redundant vars are removed')
    print ('Num of variants =', len(var_list[MAIN_STRAND]) + len(var_list[ALT_STRAND]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--var',
        help='original var file'
    )
    parser.add_argument(
        '-o', '--out_var',
        help='filtered var file'
    )
    args = parser.parse_args()
    var_fn = args.var
    out_fn = args.out_var
    remove_redundant_vars(var_fn, out_fn)
