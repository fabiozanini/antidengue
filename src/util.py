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
    import numpy as np
    #FIXME: is this how the Distance package does it? Certainly this is the
    # formal definition but one never knows... check with David
    arr1 = np.fromstring(seq1, 'S1')
    arr2 = np.fromstring(seq2, 'S1')
    return (arr1 != arr2).sum()


def write_tree_to_json(tree, output_file,
                       metadata_tree=[],
                       metadata_nodes=[],
                       children_attrname="children"):
    '''Write Biopython tree to JSON'''

    tree_dict = {}
    metadata_tree = [m for m in metadata_tree if m not in ['tree']]
    for field in metadata_tree:
        if hasattr(tree, field):
            tree_dict[field] = getattr(tree, field)


    not_metadata = ['clades', children_attrname]
    metadata_nodes = [m for m in metadata_nodes if m not in not_metadata]
    def convert_to_dict(node):
        d = {}
        for field in metadata_nodes:
            if hasattr(node, field):
                d[field] = getattr(node, field)
        d[children_attrname] = [convert_to_dict(c) for c in node.clades]
        return d

    tree_dict['tree'] = convert_to_dict(tree.root)

    import json
    with open(output_file, 'w') as handle:
        json.dump(tree_dict, handle, indent=1)

