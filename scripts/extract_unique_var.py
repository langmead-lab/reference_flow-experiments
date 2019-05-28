'''
This script reads two .var files and outputs unique variants
with regard to each .var file
'''
import argparse
from build_erg import VarInfo

parser = argparse.ArgumentParser()
parser.add_argument(
    '-v1', '--fn_var1',
    help='var file 1'
)
parser.add_argument(
    '-v2', '--fn_var2',
    help='var file 2'
)
parser.add_argument(
    '-i', '--invert', type=int, default=0,
    help='set 1 to invert haplotype A -> B, B -> A (0)'
)
parser.add_argument(
    '-u', '--unique', type=int, default=0,
    help='set 1 to output unique var files (0)'
)
args = parser.parse_args()
fn_var1 = args.fn_var1
fn_var2 = args.fn_var2
invert = args.invert
unique = args.unique

def build_var_dict(f, dict_rev=None):
    dict_var = {}
    for line in f:
        v = VarInfo(line)
        if dict_rev:
            v.strand = dict_rev[v.strand]
        dict_key = v.strand + '_' + v.chrm + '_' + str(v.ref_pos)
        if dict_var.get(dict_key):
            #: conflict of a snp and indel -> keep indel
            if dict_var[dict_key].vtype == 'INDEL' and v.vtype == 'SNP':
                pass
            elif dict_var[dict_key].vtype == 'SNP' and v.vtype == 'INDEL':
                dict_var[dict_key] = v
            else:
                print ('ERROR: dictionary conflict with \n{0}'.format(v.line))
                print (dict_var.get(dict_key).line)
        else:
            dict_var[dict_key] = v
    return dict_var

#: TODO two combinations of haplotypes should both be reported (AB and BA)
def read_file_and_get_unique_var(f, dict_var1):
    #: complete var2 dictionary
    dict_var2 = {}
    dict_var1_unique = dict_var1
    dict_var2_unique = {}
    for line in f:
        v = VarInfo(line)
        dict_key = v.strand + '_' + v.chrm + '_' + str(v.ref_pos)
        if dict_var2.get(dict_key):
            #: conflict of a snp and indel -> keep indel
            if dict_var2[dict_key].vtype == 'INDEL' and v.vtype == 'SNP':
                pass
            elif dict_var2[dict_key].vtype == 'SNP' and v.vtype == 'INDEL':
                dict_var2[dict_key] = v
            else:
                print ('ERROR: dictionary conflict with \n{0}'.format(v.line))
                print (dict_var2.get(dict_key).line)
            # print ('ERROR: dictionary conflict with {}'.format(v.line))
        else:
            dict_var2[dict_key] = v

        #: not unique
        if dict_var1.get(dict_key):
            # print ('not unique @ {}'.format(dict_key))
            del dict_var1_unique[dict_key]
        else:
            dict_var2_unique[dict_key] = v
    return dict_var2, dict_var1_unique, dict_var2_unique

dict_rev = {'A':'B', 'B':'A'}

f_var1 = open(fn_var1, 'r')
if invert == 1:
    print ('*** Invert is on, haplotype A -> B, B -> A ***')
    dict_var1 = build_var_dict(f_var1, dict_rev=dict_rev)
else:
    dict_var1 = build_var_dict(f_var1, dict_rev=None)
f_var1.close()
print ('Length of dict_var1 = {}'.format(len(dict_var1)))

f_var2 = open(fn_var2, 'r')
dict_var2, dict_var1_unique, dict_var2_unique = read_file_and_get_unique_var(f_var2, dict_var1)
f_var2.close()
print ('Length of dict_var2 = {}'.format(len(dict_var2)))
print ('Length of dict_var1_unique = {}'.format(len(dict_var1_unique)))
print ('Length of dict_var2_unique = {}'.format(len(dict_var2_unique)))

if unique:
    fn_unique1 = fn_var1.split('.')[0] + '_' + fn_var2.split('.')[0] + 'unique.var'
    fn_unique2 = fn_var2.split('.')[0] + '_' + fn_var1.split('.')[0] + 'unique.var'
    f_unique1 = open(fn_unique1, 'w')
    f_unique2 = open(fn_unique2, 'w')
    for i in dict_var1_unique.keys():
        f_unique1.write(dict_var1_unique[i].line + '\n')
    for i in dict_var2_unique.keys():
        f_unique2.write(dict_var2_unique[i].line + '\n')