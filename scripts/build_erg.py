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

def write_erg(var_list, main_genome, f_len, test_genome):
    '''
    Write one ERG seq
    '''
    erg = ''
    main_pos_list = []
    erg_start_pos = 0
    offset = 0
    # assume same chromosom
    chrm = var_list[0].chrm

    for i, v in enumerate(var_list):
        if v.strand == MAIN_STRAND:
            # v: ref hapA
            main_pos = v.alt_pos
            main_pos_list.append(main_pos)
            erg_main_allele = v.alt_allele
            erg_alt_allele = v.ref_allele
            offset = -v.offset + v.cor_offset
        else:
            # v: ref hapB
            main_pos = v.alt_pos - v.offset + v.cor_offset
            main_pos_list.append(main_pos)
            erg_main_allele = v.ref_allele
            erg_alt_allele = v.alt_allele
            offset = v.offset - v.cor_offset
        if i == 0:
            erg_start_pos = main_pos - f_len
            erg += main_genome[main_pos - f_len : main_pos]
        else:
            erg += main_genome[
                main_pos_list[i-1] + 1 : main_pos
            ]
        erg += erg_alt_allele
    erg += main_genome[main_pos_list[-1] + 1 : main_pos_list[-1] + 1 + f_len]

    erg_start_pos = erg_start_pos + offset
    erg_end_pos = main_pos_list[-1] + 1 + f_len + offset

    if test_genome:
        full_g = test_genome[erg_start_pos : erg_end_pos]
        if erg == full_g:
            print ('pass')
            return
        print ('not pass')
        print (erg_start_pos, erg_end_pos)
        print ('seq:\t', erg)
        print ('f_seq:\t', full_g)
        for i, v in enumerate(var_list):
            print (v.line)
        input()

    # write erg
    print (
        '>%sB-erg-%s-%s' % 
        (chrm, erg_start_pos, erg_end_pos)
    )
    print (erg)

def build_index(var_list, main_index, alt_index):
    existed_pos_list = []
    for i, v in enumerate(var_list):
        pos = v.alt_pos
        if pos in existed_pos_list:
            continue
        else:
            existed_pos_list.append(pos)
        c_pos = pos - v.offset + v.cor_offset
        if v.strand == MAIN_STRAND:
            main_index[pos] = [c_pos, v.vtype, v.ref_allele, v.alt_allele]
            alt_index[c_pos] = [pos, v.vtype, v.alt_allele, v.ref_allele]
        else:
            alt_index[pos] = [c_pos, v.vtype, v.ref_allele, v.alt_allele]
            main_index[c_pos] = [pos, v.vtype, v.alt_allele, v.ref_allele]
    if len(main_index) != len(alt_index):
        for i, v in enumerate(var_list):
            print (v.line)
        input()
    return main_index, alt_index

def build_erg(
    main_genome, 
    test_genome,
    var_fn, 
    hap_mode, 
    f_len, 
    mode
):
    var_f = open(var_fn, 'r')
    var_list = []
    # indexes for 'index' mode
    main_index = {}
    alt_index = {}
    if mode != 'erg' and mode != 'index':
        print ('Error: incorrect mode (%s)' % mode)
        print ('supported modes = [erg, index]')
        exit()

    for line in var_f:
        vinfo = VarInfo(line)
        if len(var_list) > 0:
            prev_var = var_list[len(var_list) - 1]
        else:
            prev_var = 0
        if __debug__:
            print (prev_var)
            print (vinfo.line)
        if vinfo.strand == MAIN_STRAND:
            if len(var_list) > 0 and \
                vinfo.alt_pos > prev_var.alt_pos + 2 * f_len:
                # write previous vars
                if mode == 'erg':
                    write_erg(var_list, main_genome, f_len, test_genome)
                else:
                    build_index(var_list, main_index, alt_index)
                # reset
                var_list = [vinfo]
            elif len(var_list) > 0:
                var_list.append(vinfo)
            else: # len == 0
                var_list = [vinfo]
            continue
        # var not main strand
        if vinfo.samevar(prev_var):
            # delete a var if it occurs in both strands
            if __debug__:
                print ('samevar')
                print ('A*', prev_var.line)
                print ('B*', vinfo.line)
            del var_list[len(var_list) - 1]
            continue
        else:
            if prev_var == 0:
                var_list.append(vinfo)
            elif vinfo.alt_pos > prev_var.alt_pos + 2 * f_len:
                if mode == 'erg':
                    write_erg(var_list, main_genome, f_len, test_genome)
                else:
                    build_index(var_list, main_index, alt_index)
                var_list = [vinfo]
            else:
                var_list.append(vinfo)
    if mode == 'erg':
        write_erg(var_list, main_genome, f_len, test_genome)
    else:
        build_index(var_list, main_index, alt_index)
    var_f.close()
    if mode == 'index':
        return main_index, alt_index

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-m', '--main_genome',
        help='main genome file'
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
    tg_fn = args.test_genome
    var_fn = args.var
    hap_mode = args.hap
    print_mg = args.print_mg
    f_len = args.flanking_len

    return mg_fn, tg_fn, var_fn, hap_mode, f_len, print_mg

if __name__ == '__main__':
    mg_fn, tg_fn, var_fn, hap_mode, f_len, print_mg = parse_args()
    main_genome = read_genome(mg_fn, print_mg)
    if tg_fn:
        test_genome = read_genome(tg_fn, None)
    else:
        test_genome = None
    build_erg(
        main_genome=main_genome, 
        test_genome=test_genome, 
        var_fn=var_fn, 
        hap_mode=hap_mode, 
        f_len=f_len, 
        mode='erg'
    )
