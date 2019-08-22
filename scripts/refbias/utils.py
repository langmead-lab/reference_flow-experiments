def get_paths_from_list(fn_list, prefix = '', suffix = ''):
    l = []
    with open(fn_list, 'r') as f:
        for line in f:
            l.append(prefix + line.rstrip() + suffix)
    return l

def get_het_from_list_format(het_in_list):
    hets = [int(i) for i in het_in_list.strip('[]').split(',')]
    #: check if all the positions for a het site are the same
    if len(hets) > 1:
        for e in hets:
            assert e == hets[0]
    return int(hets[0])
