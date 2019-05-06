#sen=False
sen=True
state = '0'
with open('results_premerged.log', 'r') as f:
    for line in f:
        if line.startswith('sample'):
            print (line.rstrip())
        elif line.startswith('sensitivity_all') and state == '0' and sen:
            true_first = int(line.split()[3][1:])
            all_first = int(line.split()[5][:-1])
            state = '1'
        elif line.startswith('sensitivity_all') and state == '1' and sen:
            true_second = int(line.split()[3][1:])
            all_second = int(line.split()[5][:-1])
            true_total = true_first + true_second
            all_total = all_first + all_second
            state = '0'
            print ('sensitivity_all      = {0:.4%} ({1} / {2})'.format(true_total/all_total, true_total, all_total))
        elif line.startswith('unaligned            =') and state == '0' and (sen == False):
            #un_first = int(line.split()[2])
            un_first = int(line.split()[3][1:])
            state = '1'
        elif line.startswith('unaligned            =') and state == '1' and (sen == False):
            #un_second = int(line.split()[2])
            un_second = int(line.split()[3][1:])
            un_total = un_first + un_second
            state = '0'
            print ('Unaligned: {0}'.format(un_total))
