# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/05/16
content:    Utility functions
'''
# Functions
def mkdirs(newdir):
    """ Create a directory and all parent folders.
        Features:
        - parent directories will be created
        - if directory already exists, then do nothing
        - if there is another filsystem object with the same name, raise an exception
    """
    import os
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("cannot create directory, file already exists: '%s'" % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            mkdirs(head)
        if tail:
            os.mkdir(newdir)


def hamming(seq1, seq2):
    '''Poor man's version of the hamming function'''
    #FIXME: is this how the Distance package does it? Certainly this is the
    # formal definition but one never knows... check with David
    arr1 = np.fromstring(seq1, 'S1')
    arr2 = np.fromstring(seq2, 'S1')
    return (arr1 != arr2).sum()
