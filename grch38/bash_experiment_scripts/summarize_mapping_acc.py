import argparse

def organize_accuracy(fn_input, fn_output):
    with open(fn_input, 'r') as f:
        list_tp = []
        list_all = []
        for line in f:
            if line.count('sensitivity_all') > 0:
                line = line.split()
                list_tp.append(int(line[3][1:]))
                list_all.append(int(line[5][:-1]))
    f_out = open(fn_output, 'w')
    f_out.write('{0}\n{1}\n{2}\n'.format(sum(list_tp), sum(list_all), sum(list_tp) / sum(list_all)))
    return


def organize_unaligned(fn_input, fn_output):
    with open(fn_input, 'r') as f:
        list_unaligned = []
        list_all = []
        for line in f:
            if line.count('unaligned            =') > 0:
                line = line.split()
                list_unaligned.append(int(line[3][1:]))
                list_all.append(int(line[5][:-1]))
    f_out = open(fn_output, 'w')
    f_out.write('{0}\n{1}\n{2}\n'.format(sum(list_unaligned), sum(list_all), sum(list_unaligned) / sum(list_all)))
    return


def organize_numincorrect(fn_input, fn_output):
    with open(fn_input, 'r') as f:
        list_tp = []
        list_unaligned = []
        list_all = []
        for line in f:
            if line.count('sensitivity_all') > 0:
                line = line.split()
                list_tp.append(int(line[3][1:]))
            if line.count('unaligned            =') > 0:
                line = line.split()
                list_unaligned.append(int(line[3][1:]))
                list_all.append(int(line[5][:-1]))
        list_numincorrect = [list_all[i] - list_tp[i] - list_unaligned[i] for i, _ in enumerate(list_all)]
    f_out = open(fn_output, 'w')
    f_out.write('{0}\n{1}\n{2}\n'.format(sum(list_numincorrect), sum(list_all), sum(list_numincorrect) / sum(list_all)))
    return

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input',
        help='Input file')
    parser.add_argument(
        '-o', '--output',
        help='Output file')
    parser.add_argument(
        '-m', '--mode', default="sensitivity",
        help='"sensitivity": alignment sensitivity;'
            '"unaligned": number of unaligned reads'
            '"num_incorrect": number of incorrectly aligned reads')
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_args()
    if args.mode == 'sensitivity':
        organize_accuracy(args.input, args.output)
    elif args.mode == 'unaligned':
        organize_unaligned(args.input, args.output)
    elif args.mode == 'num_incorrect':
        organize_numincorrect(args.input, args.output)
    else:
        print ('invalid mode:', args.mode)
        exit()

