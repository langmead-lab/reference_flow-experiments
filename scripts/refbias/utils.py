import os
def get_paths_from_list(fn_list, prefix = '', suffix = '', auto_prefix = None):
    '''
    @argument fn_list (STR: filename)
        a file that stores a list of ID/path(s)
    @argument prefix (STR)
        prefix added to paths
    @argument suffix (STR)
        suffix added to paths
    @argument auto_prefix (LIST)
        set True to add common dirname from a list as prefix
    @return processed list of path(s)
    '''
    l = []
    with open(fn_list, 'r') as f:
        for line in f:
            l.append(prefix + line.rstrip() + suffix)
    if auto_prefix != None:
        prefix = os.path.dirname(os.path.commonprefix(auto_prefix)) + '/'
        l2 = [prefix + i for i in l]
        return l2
    else:
        return l

def get_het_from_list_format(het_in_list):
    hets = [int(i) for i in het_in_list.strip('[]').split(',')]
    #: check if all the positions for a het site are the same
    if len(hets) > 1:
        for e in hets:
            assert e == hets[0]
    return int(hets[0])
