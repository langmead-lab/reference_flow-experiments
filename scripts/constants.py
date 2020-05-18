MAIN_STRAND = 'A'
ALT_STRAND = 'B'
MAIN_HAP = 'hap' + MAIN_STRAND
ALT_HAP = 'hap' + ALT_STRAND

STEP = 1000
READ_LEN = 100

def set_chrom(chrom):
    global CHROM, MAIN_CHROM, ALT_CHROM
    CHROM = chrom
    MAIN_CHROM = chrom + MAIN_STRAND
    ALT_CHROM = chrom + ALT_STRAND

def set_step(step):
    global STEP
    STEP = step

def set_read_len(read_len):
    global READ_LEN
    READ_LEN = read_len
