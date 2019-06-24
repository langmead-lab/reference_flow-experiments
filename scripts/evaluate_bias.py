from build_erg import read_var
from analyze_sam import parse_line
from analyze_diploid_indels import build_all_indexes
import constants
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument(
    '-v', '--fn_var',
    help='target .var file where all the variants are HETs'
)
parser.add_argument(
    '-s', '--fn_sam',
    help='target aligned reads (.sam) file'
)
parser.add_argument(
    '-g', '--fn_golden',
    help='golden log file (.fq.sam) of target reads'
)
parser.add_argument(
    '-vr', '--fn_var_reads',
    help='.var file specifying variants\
         included in the genome where reads are simluated from'
)
parser.add_argument(
    '-vs', '--fn_var_sample',
    help='.var file specifying variants for the reference used for alignment'
)
parser.add_argument(
    '-read_len', '--read_len', type=int, default=100,
    help='read length (100)'
)
parser.add_argument(
    '-step', '--step', type=int, default=1000,
    help='step size (1000)'
)
args = parser.parse_args()
fn_var = args.fn_var
fn_sam = args.fn_sam
fn_golden = args.fn_golden
fn_var_reads = args.fn_var_reads
fn_var_sample = args.fn_var_sample
read_len = args.read_len
step = args.step

def extract_reads_overlapping_hets(
    fn_golden,
    set_var_ref_pos
):
    f_golden = open(fn_golden, 'r')
    # pos_prev = 0
    list_name = []
    for line in f_golden:
        name, info = parse_line(line)
        if name == 'header' or name.count(target_strand) == 0:
            continue
        
        for i in range(info.pos, info.pos+100):
            if i in set_var_ref_pos:
                list_name.append(name)
                # info.print()
                # input()
    return list_name
    # print (len(list_name))
    # print (list_name[:20])

def calc_offsets(
    info,
    sample_main_offset_index,
    sample_alt_offset_index
):
    '''
    Calculates offsets for aligned reads and maps the position 
    to standard reference coordinate system.

    This function is extracted from analyze_diploid_indels -- 
    diploid_compare().
    '''
    i_low = int(info.pos / constants.STEP)
    i_high = math.ceil(info.pos / constants.STEP)
    if info.chrm == constants.MAIN_CHROM or info.chrm == constants.CHROM:
        if i_low >= len(sample_main_offset_index):
            sample_offset_low = sample_main_offset_index[len(sample_main_offset_index) - 1]
        else:
            sample_offset_low = sample_main_offset_index[i_low]
        if i_high >= len(sample_main_offset_index):
            sample_offset_high = sample_main_offset_index[len(sample_main_offset_index) - 1]
        else:
            sample_offset_high = sample_main_offset_index[i_high]
    elif info.chrm == constants.ALT_CHROM:
        if i_low >= len(sample_alt_offset_index):
            sample_offset_low = sample_alt_offset_index[len(sample_alt_offset_index) - 1]
        else:
            sample_offset_low = sample_alt_offset_index[i_low]
        if i_high >= len(sample_alt_offset_index):
            sample_offset_high = sample_alt_offset_index[len(sample_alt_offset_index) - 1]
        else:
            sample_offset_high = sample_alt_offset_index[i_high]
    else:
        print ('Error: invalid chrm', info.chrm, constants.MAIN_CHROM, constants.ALT_CHROM)
        exit()
    sample_offsets = [sample_offset_low, sample_offset_high]
    return sample_offsets

if __name__ == '__main__':
    constants.set_chrom('21')
    constants.set_step(step)
    constants.set_read_len(read_len)
    #: removes INDELS, ALT (homozygous) and TRI-ALLELIC variants
    #: list_var is based on standard reference coordinate system
    list_var = read_var(
        fn_var,
        remove_conflict=True,
        remove_homo_alt=True,
        remove_indel=True,
        remove_tri_allelic=True,
        MAIN_STRAND=constants.MAIN_STRAND,
        ALT_STRAND=constants.ALT_STRAND
    )
    list_var_ref_pos = []
    list_alt_strand = []
    for i in list_var:
        list_var_ref_pos.append(i.ref_pos)
        list_alt_strand.append(i.strand)

    #: dict_var_ref_pos
    #:      key  : POS (wrt to standard reference genome)
    #:      value: (ALT_STRAND, REF_COUNT, ALT_COUNT)
    dict_var_ref_pos = {}
    for i, p in enumerate(list_var_ref_pos):
        dict_var_ref_pos[p] = (list_alt_strand[i], 0, 0)

    all_indexes = build_all_indexes(
        var_reads_fn=fn_var_reads,
        var_sample_fn=fn_var_sample,
        personalized=2,
        step=step,
        MAIN_STRAND=constants.MAIN_STRAND,
        ALT_STRAND=constants.ALT_STRAND
    )
    # main_index, alt_index, reads_main_offset_index, reads_alt_offset_index, sample_main_offset_index, sample_alt_offset_index = all_indexes
    _, _, reads_main_offset_index, reads_alt_offset_index, sample_main_offset_index, sample_alt_offset_index = all_indexes
    
    f_sam = open(fn_sam, 'r')
    for line in f_sam:
        name, info = parse_line(line)
        if name == 'header' or info.is_unaligned():
            continue
        offsets = calc_offsets(
            info,
            sample_main_offset_index,
            sample_alt_offset_index
        )
        # print (name)
        pos_lower = info.pos - offsets[0]
        pos_upper = info.pos - offsets[1]
        set_lower = set(range(pos_lower, pos_lower+constants.READ_LEN))
        set_upper = set(range(pos_upper, pos_upper+constants.READ_LEN))
        range_pos = set_lower.union(set_upper)
        for pos in range_pos:
            if pos in dict_var_ref_pos:
                # print (name)
                # print (dict_var_ref_pos[pos])
                #: read belong to REF
                if (name.count(constants.MAIN_HAP) > 0 and dict_var_ref_pos[pos][0] != constants.MAIN_STRAND) or \
                    (name.count(constants.ALT_HAP) > 0 and dict_var_ref_pos[pos][0] != constants.ALT_STRAND):
                    dict_var_ref_pos[pos] = (dict_var_ref_pos[pos][0], dict_var_ref_pos[pos][1]+1, dict_var_ref_pos[pos][2])
                #: read belong to ALT
                elif (name.count(constants.MAIN_HAP) > 0 and dict_var_ref_pos[pos][0] == constants.MAIN_STRAND) or \
                    (name.count(constants.ALT_HAP) > 0 and dict_var_ref_pos[pos][0] == constants.ALT_STRAND):
                    dict_var_ref_pos[pos] = (dict_var_ref_pos[pos][0], dict_var_ref_pos[pos][1], dict_var_ref_pos[pos][2]+1)
                else:
                    print ('ERROR: haplotype {} not found!'.format(name))
                # print (dict_var_ref_pos[pos])
                # input()
        # input (offsets)
    
    count_ref = [i[1] for i in dict_var_ref_pos.values()]
    count_alt = [i[2] for i in dict_var_ref_pos.values()]

    # for i in list_var_ref_pos:
    #     extract_reads_overlapping_hets(
    #         fn_golden,
    #         i,
    #         constants.MAIN_HAP,
    #         constants.ALT_HAP
    #     )
    # set_var_ref_pos = set(list_var_ref_pos)
    # list_name_A = extract_reads_overlapping_hets(fn_golden, set_var_ref_pos, 'hapA')
    # list_name_B = extract_reads_overlapping_hets(fn_golden, set_var_ref_pos, 'hapB')
    
    
    
    
    # print ('removes CONFLICT, INDEL, ALT and TRI', len(list_var))
    # for i in list_var:
    #     print (i.line)
    # list_var = read_var(fn_var, remove_conflict=True, remove_homo_alt=True, remove_indel=True, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('removes CONFLICT, INDEL and ALT', len(list_var))
    # list_var = read_var(fn_var, remove_conflict=True, remove_homo_alt=False, remove_indel=True, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('removes CONFLICT and INDEL', len(list_var))
    # list_var = read_var(fn_var, remove_conflict=True, remove_homo_alt=False, remove_indel=False, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('removes CONFLICT', len(list_var))
    # list_var = read_var(fn_var, remove_conflict=False, remove_homo_alt=False, remove_indel=False, MAIN_STRAND=MAIN_STRAND, ALT_STRAND=ALT_STRAND)
    # print ('all', len(list_var))