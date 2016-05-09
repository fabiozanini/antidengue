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

