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
            a.ref_allele == self.ref_allele and \
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
    USE_GOLDEN_ERG = True
    if test_genome:
        full_g = test_genome[alt_start_pos : alt_end_pos]

        if USE_GOLDEN_ERG:
            print (
                '>%sB-erg-%s-%s' % 
                (chrm, alt_start_pos, alt_end_pos)
            )
            print (full_g)
        else:
            if erg == full_g:
                print ('pass')
                return
            print ('fail')
        
        if TEST_DETAILS:
            print (alt_start_pos, alt_end_pos)
            print ('erg:\t', erg)
            print ('golden:\t', full_g)
            print ('ref:\t', erg_ref)
            for v in var_list:
                print (v.line)
            input()
        return len(full_g)
    else:
        # write erg
        print (
            '>%sB-erg-%s-%s' % 
            (chrm, alt_start_pos, alt_end_pos)
        )
        print (erg)
        return len(erg)

def build_erg(
    main_genome, 
    ref_genome,
    test_genome,
    var_list,
    hap_mode, 
    f_len
):
    '''
    Builds a pair of erg using genomes and a var file.
    '''
    tmp_var_list = []
    num_erg = 0
    total_len_erg = 0

    for vinfo in var_list:
        if len(tmp_var_list) > 0:
            prev_var = tmp_var_list[len(tmp_var_list) - 1]
        else:
            prev_var = 0
        if len(tmp_var_list) > 0 and \
            vinfo.ref_pos > prev_var.ref_pos + 2 * f_len:
            # write previous vars
            len_erg = write_erg(tmp_var_list, main_genome, f_len, test_genome, ref_genome)
            num_erg += 1
            total_len_erg += len_erg
            # reset
            tmp_var_list = [vinfo]
        elif len(tmp_var_list) > 0:
            tmp_var_list.append(vinfo)
        else: # len == 0
            tmp_var_list = [vinfo]
    len_erg = write_erg(tmp_var_list, main_genome, f_len, test_genome, ref_genome)
    num_erg += 1
    total_len_erg += len_erg
    SHOW_SUMMARY = False
    if SHOW_SUMMARY:
        print ('Num ergs =', num_erg)
        print ('Avg len of an erg =', float(total_len_erg) / num_erg)

def read_var(var_fn, remove_conflict, remove_coexist=False):
    '''
    Build a dictionary for the .var file.

    remove_conflict: 
        If set, removes conflict-represented variants. There might be other type of conflicts, but we only solve the following type here.
        
        Example:
        A   9   INDEL   333711  333728  TA  T
        A   9   SNP     333712  333729  A   T
    
    remove_coexsit:
        If set, removes variant locating on both strands.

        Example:
        A   9   INDEL   10362   10363   C   CT
        B   9   INDEL   10362   10363   C   CT
    '''
    SHOW_REMOVED = False
    count_conflict = 0
    count_coexist = 0
    var_f = open(var_fn, 'r')
    del_pos = {MAIN_STRAND:[], ALT_STRAND:[]}
    del_allele = {MAIN_STRAND:[], ALT_STRAND:[]}
    var_list = []
    for line in var_f:
        v = VarInfo(line)
        if remove_conflict == False:
            var_list.append(v)
            continue
        if remove_coexist and len(var_list) > 0:
            top_v = var_list[len(var_list) - 1]
            if v.samevar(top_v):
                var_list.pop(len(var_list) - 1)
                count_coexist += 1
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
            count_conflict += 1
            if SHOW_REMOVED:
                print ('"%s"' % v.strand, v.ref_pos, d)
        else:
            var_list.append(v)

    # Set this True to show remove information
    SHOW_REMOVE_INFO = False
    if SHOW_REMOVE_INFO:
        print (count_conflict, 'conflict vars are removed')
        print (count_coexist, 'coexisted vars are removed')
        print ('Num of variants =', len(var_list))

    # Set this True to write variants and then exit
    WRITE_VARS = False
    if WRITE_VARS:
        for v in var_list:
            print (v.line)
        exit ()
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
    
    var_list = read_var(var_fn, remove_conflict=True, remove_coexist=True)

    build_erg(
        main_genome=main_genome, 
        ref_genome=ref_genome,
        test_genome=test_genome, 
        var_list=var_list,
        hap_mode=hap_mode, 
        f_len=f_len
    )
