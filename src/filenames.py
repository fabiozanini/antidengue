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
    else:
        raise ValueError('user not recognized')


def get_sample_table_filename():
    '''Get the filename of the sample table'''
    root_fn = get_root_data_foldername()
    return root_fn + 'sample_info.csv'


def get_reads_foldername(fmt='fastq'):
    '''Get the folder name of the reads'''
    fn = get_root_data_foldername() + fmt + '/'
    return fn


def get_sample_filename(samplename, fmt='fastq'):
    '''Get the filename of a sample'''
    fn = get_reads_foldername() + samplename + '.' + fmt
    return fn


def get_reads_filenames(fmt='fastq'):
    '''Get the filenames of fastq/fasta files'''
    import glob
    fn = get_reads_foldername()
    return glob.glob(fn+'*.'+fmt)


def get_germline_db_filename(genename):
    '''Get the germline database filename'''
    root_fn = get_root_data_foldername() + 'idblastdb/'
    return root_fn + 'IGH' + genename + '_imgt.fasta'


def get_igblast_auxiliary_data_filename():
    '''Get the filename of the human auxiliary file for igblast'''
    root_fn = get_root_data_foldername() + 'idblastdb/optional_file/'
    return root_fn + 'human_gl.aux'


def get_igblast_execfile(variant='n'):
    '''Get the executable of igblastn/igblastp'''
    if 'glass' in os.getenv('HOME'):
        return '/Users/davidglass/Documents/resources/ncbi-igblast-1.4.0/bin/igblast' + variant
    elif 'fabio' in os.getenv('HOME'):
        return '/usr/bin/igblast' + variant
    else:
        raise ValueError('user not recognized')

