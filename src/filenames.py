# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/05/16
content:    Support module for antidengue
'''
# Modules
import os



# Function
def get_root_data_foldername():
    '''Get the foldername of the data'''
    if 'glass' in os.getenv('HOME'):
        return '/Users/davidglass/antidengue/data/'
    elif 'fabio' in os.getenv('HOME'):
        return os.getenv('HOME')+'university/postdoc/antidengue/data/'


def get_sample_table_filename():
    '''Get the filename of the sample table'''
    root_fn = get_root_data_foldername()
    return root_fn + 'sample_info.csv'


def get_reads_filenames(fmt='fastq'):
    '''Get the filenames of fastq/fasta files'''
    import glob
    fn = get_root_data_foldername() + fmt + '/'
    return glob.glob(fn+'*.'+fmt)
