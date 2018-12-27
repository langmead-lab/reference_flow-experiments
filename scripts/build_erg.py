'''
Build enhanced reference genome using a 
pre-identified haplotype and a var file
'''
import argparse

# In var file, the strand occurs
# earlier is the main_strand
MAIN_STRAND = 'A'
ALT_STRAND = 'B'

class VarInfo():
    '''
    Records information of a var line
    '''
    strand = ''
    chrm = ''
    vtype = ''
    ref_pos = 0
    alt_pos = 0
    ref_allele = ''
    alt_allele = ''
    offset = 0
    cor_offset = 0 # the offset of the other strand
    line = ''

    def __init__(self, line):
        self.line = line[: line.find('\\')]
        line = line.split()
        self.strand = line[0]
        self.chrm = line[1]
        self.vtype = line[2]
        self.ref_pos = int(line[3])
        self.alt_pos = int(line[4])
        self.ref_allele = line[5]
        self.alt_allele = line[6]
        self.offset = int(line[7])
        self.cor_offset = int(line[8])
    
    def is_indel(self):
        return (self.vtype == 'INDEL')

    def is_ins(self):
        if self.vtype == 'SNP':
            return False
        if len(self.ref_allele) == len(self.alt_allele):
            print ('Error: incorrectly specified indel')
            print (self.line)
            exit()
        if len(self.ref_allele) > len(self.alt_allele):
            return False
        return True

    def is_del(self):
        if self.vtype == 'SNP':
            return False
        if len(self.ref_allele) == len(self.alt_allele):
            print ('Error: incorrectly specified indel')
            print (self.line)
            exit()
        if len(self.ref_allele) < len(self.alt_allele):
            return False
        return True

    def samevar(self, a):
        if type(a) != type(self):
            return False
        if a.vtype == self.vtype and \
            a.ref_pos == self.ref_pos and \
            a.alt_allele == self.alt_allele:
            return True
        return False


def read_genome(fn, print_main=False):
    with open(fn, 'r') as f:
        # seq[0] is empty to fit vcf coordinate (starting from 1)
        seq = '^'
        for line in f:
            line = line[: line.find('\\')]
            if print_main:
                print (line)
            if line.startswith('>') == False:
                seq += line
    return seq

def write_erg(var_list, main_genome, f_len, test_genome, ref_genome):
    '''
    Write one ERG seq
    '''
    erg = ''
    # offset_init = 0
    # assume same chromosom
    chrm = var_list[0].chrm
    
    # coordinate on REF
    v = var_list[0]
    if v.strand == MAIN_STRAND:
        offset_start = v.cor_offset
    else:
        offset_start = v.offset
    ref_start_pos = v.ref_pos - f_len

    v = var_list[len(var_list) - 1]
    ins_len = len(v.alt_allele) - len(v.ref_allele)
    if v.strand == MAIN_STRAND:
        offset_end = v.cor_offset
    else:
        offset_end = v.offset + ins_len
    ref_end_pos = v.ref_pos + f_len + len(v.alt_allele)
    erg_ref = ref_genome[ref_start_pos : ref_end_pos]
    
    nxt_start_pos = ref_start_pos
    for v in var_list:
        if v.strand == ALT_STRAND:
            # v: ref hapB
            ref_pos = v.ref_pos
            erg += ref_genome[nxt_start_pos : ref_pos]
            erg += v.alt_allele
            nxt_start_pos = ref_pos + len(v.ref_allele)
    if len(erg) == 0:
        erg = erg_ref
    else:
        erg += ref_genome[nxt_start_pos : ref_end_pos]

    # coordinate on ALT_STRAND
    alt_start_pos = ref_start_pos + offset_start
    alt_end_pos = ref_end_pos + offset_end

    # set this True for testing mode
    TEST_DETAILS = False
    if test_genome:
        full_g = test_genome[alt_start_pos : alt_end_pos]
        
        if erg == full_g:
            print ('pass')
            return
        print ('fail')
        
        if TEST_DETAILS:
            print (alt_start_pos, alt_end_pos)
            print ('erg:\t', erg)
            print ('f_seq:\t', full_g)
            for i, v in enumerate(var_list):
                print (v.line)
            input()
    else:
        # write erg
        print (
            '>%sB-erg-%s-%s' % 
            (chrm, alt_start_pos, alt_end_pos)
        )
        print (erg)

# TODO
def build_index(var_list, main_index, alt_index):
    existed_pos_list = []
    for v in var_list:
        pos = v.alt_pos
        if pos in existed_pos_list:
            continue
        else:
            existed_pos_list.append(pos)
        c_pos = pos - v.offset + v.cor_offset
        if v.strand == MAIN_STRAND:
            # MAIN-MAIN: with replacement
            main_index[pos] = [c_pos, v.vtype, v.ref_allele, v.alt_allele]
            # MAIN-ALT: without replacement
            if alt_index.get(c_pos) == None:
                alt_index[c_pos] = [pos, v.vtype, v.alt_allele, v.ref_allele]
        else:
            # ALT-ALT: with replacement
            alt_index[pos] = [c_pos, v.vtype, v.ref_allele, v.alt_allele]
            # ALT-MAIN: without replacement
            if main_index.get(c_pos) == None:
                main_index[c_pos] = [pos, v.vtype, v.alt_allele, v.ref_allele]
    # TODO
    '''
    There is index imbalance problem due to 
    overlapping or conflicting variants
    '''
    SHOW_IMBALANCE = False
    if SHOW_IMBALANCE and (len(main_index) != len(alt_index)):
        for v in var_list:
            print (v.line)
        print (main_index)
        print (alt_index)
        input()
    return main_index, alt_index

def build_erg(
    main_genome, 
    ref_genome,
    test_genome,
    var_list,
    hap_mode, 
    f_len, 
    mode
):
    '''
    modes:
        'erg' - build erg
        'index' - build indexes for accuracy measurement
    '''
    tmp_var_list = []
    # indexes for 'index' mode
    main_index = {}
    alt_index = {}
    if mode != 'erg' and mode != 'index':
        print ('Error: incorrect mode (%s)' % mode)
        print ('supported modes = [erg, index]')
        exit()

    for vinfo in var_list:
        if len(tmp_var_list) > 0:
            prev_var = tmp_var_list[len(tmp_var_list) - 1]
        else:
            prev_var = 0
        if __debug__:
            print (prev_var)
            print (vinfo.line)
        if len(tmp_var_list) > 0 and \
            vinfo.ref_pos > prev_var.ref_pos + 2 * f_len:
            # write previous vars
            if mode == 'erg':
                write_erg(tmp_var_list, main_genome, f_len, test_genome, ref_genome)
            else:
                build_index(tmp_var_list, main_index, alt_index)
            # reset
            tmp_var_list = [vinfo]
        elif len(tmp_var_list) > 0:
            tmp_var_list.append(vinfo)
        else: # len == 0
            tmp_var_list = [vinfo]
    if mode == 'erg':
        write_erg(tmp_var_list, main_genome, f_len, test_genome, ref_genome)
    else:
        build_index(tmp_var_list, main_index, alt_index)
    if mode == 'index':
        return main_index, alt_index

def read_var(var_fn, remove_redundant):
    '''
    Build a dictionary for the .var file
    '''
    SHOW_REMOVED = False
    removed_count = 0
    var_f = open(var_fn, 'r')
    del_pos = {MAIN_STRAND:[], ALT_STRAND:[]}
    del_allele = {MAIN_STRAND:[], ALT_STRAND:[]}
    var_list = []
    for line in var_f:
        v = VarInfo(line)
        if remove_redundant == False:
            var_list.append(v)
            continue
        if v.is_del():
            vd = del_pos[v.strand]
            if len(vd) > 0 and vd[len(vd) - 1] < v.ref_pos:
                del_pos[v.strand] = []
                del_allele[v.strand] = []
            for i in range(len(v.ref_allele)):
                del_pos[v.strand].append(v.ref_pos + i)
                del_allele[v.strand].append(v.ref_allele[i])
            var_list.append(v)
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
            var_list.append(v)

    # Set this True to show remove information
    SHOW_REMOVE_INFO = False
    if SHOW_REMOVE_INFO:
        print (removed_count, 'redundant vars are removed')
        print ('Num of variants =', len(var_list))

    return var_list

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m', '--main_genome',
        help='main genome file'
    )
    parser.add_argument(
        '-r', '--ref_genome',
        help='ref genome file'
    )
    parser.add_argument(
        '-t', '--test_genome',
        default=None,
        help='test (alt) genome file'
    )
    parser.add_argument(
        '-v', '--var',
        help='var file'
    )
    parser.add_argument(
        '--hap', type=int,
        default=1,
        help='set 1 to enable haplotype-based erg [1]'
    )
    parser.add_argument(
        '--print_mg', type=int,
        default=0,
        help='set 1 to print main genome [0]'
    )
    parser.add_argument(
        '-f', '--flanking_len', type=int,
        default=100,
        help='the length of flanking regions, total length of an erg is 2*f+1 [100]'
    )
    args = parser.parse_args()
    mg_fn = args.main_genome
    rg_fn = args.ref_genome
    tg_fn = args.test_genome
    var_fn = args.var
    hap_mode = args.hap
    print_mg = args.print_mg
    f_len = args.flanking_len

    return mg_fn, rg_fn, tg_fn, var_fn, hap_mode, f_len, print_mg

if __name__ == '__main__':
    mg_fn, rg_fn, tg_fn, var_fn, hap_mode, f_len, print_mg = parse_args()

    ref_genome = read_genome(rg_fn, None)
    main_genome = read_genome(mg_fn, print_mg)
    if tg_fn:
        test_genome = read_genome(tg_fn, None)
    else:
        test_genome = None
    
    var_list = read_var(var_fn, remove_redundant=True)

    build_erg(
        main_genome=main_genome, 
        ref_genome=ref_genome,
        test_genome=test_genome, 
        var_list=var_list,
        hap_mode=hap_mode, 
        f_len=f_len, 
        mode='erg'
    )
