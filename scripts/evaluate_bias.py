from build_erg import read_var
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    '-v', '--fn_var',
    help='target .var file'
)
args = parser.parse_args()
fn_var = args.fn_var
MAIN_STRAND = 'A'
ALT_STRAND = 'B'

if __name__ == '__main__':
    #: removes INDELS, (homozygous) ALT #TODO tri-allelic variants
    var_list = read_var(fn_var, remove_conflict=True, remove_homo_alt=True, remove_indel=True, remove_tri_allelic=True, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('removes CONFLICT, INDEL, ALT and TRI', len(var_list))
    for i in var_list:
        print (i.line)
    # var_list = read_var(fn_var, remove_conflict=True, remove_homo_alt=True, remove_indel=True, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('removes CONFLICT, INDEL and ALT', len(var_list))
    # var_list = read_var(fn_var, remove_conflict=True, remove_homo_alt=False, remove_indel=True, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('removes CONFLICT and INDEL', len(var_list))
    # var_list = read_var(fn_var, remove_conflict=True, remove_homo_alt=False, remove_indel=False, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('removes CONFLICT', len(var_list))
    # var_list = read_var(fn_var, remove_conflict=False, remove_homo_alt=False, remove_indel=False, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('all', len(var_list))