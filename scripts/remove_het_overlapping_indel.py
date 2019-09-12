import sys

class QuickVcfInfo:
    chrom = ''
    pos = 0
    len_ref = 0
    len_alt = 0
    gt = ''
    line = ''

    def __init__(self, line):
        self.line = line.rstrip()
        field = self.line.split()

        self.chrom = field[0]
        self.pos = int(field[1])
        self.len_ref = len(field[3])
        self.len_alt = len(field[4])
        self.gt = field[9]
    
    def get_var_size(self):
        return abs(self.len_ref - self.len_alt)

    def is_overlapping(self, v):
        '''
        Return values:
            -1 = different chromosome
            0  = v is overlapping or after self
            1  = self is after v
        '''
        if self.chrom != v.chrom:
            return -1
        elif self.pos + self.get_var_size() > v.pos + v.get_var_size():
            return 1
        return 0

tmp_indel = None
tmp_snv = None
for line in sys.stdin:
    if line[0] == '#':
        sys.stdout.write(line)
        continue
    try:
        info = QuickVcfInfo(line)
    except:
        continue
    if info.get_var_size() == 0: #: SNV
        if tmp_snv == None:
            if tmp_indel == None:
                tmp_snv = info
            else:
                is_op = info.is_overlapping(tmp_indel)
                if is_op == 0:
                    tmp_snv = None
                elif is_op == -1 or is_op == 1:
                    tmp_snv = info
                    tmp_indel = None
        elif info.pos > tmp_snv.pos:
            print (tmp_snv.line)
            if tmp_indel == None:
                tmp_snv = info
            else:
                is_op = info.is_overlapping(tmp_indel)
                if is_op == 0:
                    tmp_snv = None
                elif is_op == -1 or is_op == 1:
                    tmp_snv = info
                    tmp_indel = None
                else:
                    sys.stderr.write('Warning: unspecified overlapping condition \n{}\n'.format(info.line))
        #: exceptions (?)
        else:
            sys.stderr.write('Warning: overlapping SNVs \n{}\n{}\n'.format(tmp_snv.line, info.line))
            tmp_snv = info
    elif info.get_var_size() > 0: #: INDEL
        if tmp_snv == None:
            if tmp_indel == None:
                tmp_indel = info
            else:
                is_op = info.is_overlapping(tmp_indel)
                if is_op == 0:
                    pass
                elif is_op == -1 or is_op == 1:
                    tmp_indel = info
                else:
                    sys.stderr.write('Warning: unspecified overlapping condition \n{}\n'.format(info.line))
        else:
            if info.pos > tmp_snv.pos:
                print (tmp_snv.line)
            tmp_snv = None
            if tmp_indel == None:
                tmp_indel = info
            else:
                is_op = info.is_overlapping(tmp_indel)
                if is_op == 0:
                    pass
                elif is_op == -1 or is_op == 1:
                    tmp_indel = info
                else:
                    sys.stderr.write('Warning: unspecified overlapping condition \n{}\n'.format(info.line))
if tmp_snv != None:
    print (tmp_snv.line)
