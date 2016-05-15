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
        return os.getenv('HOME')+'/university/postdoc/antidengue/data/'
    else:
        raise ValueError('user not recognized')


def get_sample_table_filename():
    '''Get the filename of the sample table'''
    root_fn = get_root_data_foldername()
    return root_fn + 'sample_info.csv'


def get_lineage_table_filename():
    '''Get the filename of the sample table'''
    root_fn = get_root_data_foldername()
    return root_fn + 'lineage_info.csv'


def get_reads_foldername(fmt='fastq'):
    '''Get the folder name of the reads'''
    fn = get_root_data_foldername() + fmt + '/'
    return fn


def get_igblasted_foldername():
    '''Get the folder name of the blasted reads'''
    fn = get_root_data_foldername() + 'igblasted/'
    return fn


def get_by_sample_foldername():
    '''Get the folder name of the blasted reads'''
    fn = get_root_data_foldername() + 'by_sample/'
    return fn


def get_lineages_foldername(options=None):
    '''Get the foldername of the lineages'''
    fn = get_root_data_foldername() + 'lineages'
    if options is not None:
        if 'distance' in options:
            fn += '_'+str(options['distance'])
        if 'alignment' in options:
            fn += '/alignments'
        elif 'tree' in options:
            fn += '/trees'

    fn += '/'
    return fn


def get_lineage_filename(lineagename, options=None):
    '''Get the filename of a lineage'''
    fn = get_lineages_foldername(options=options)
    fn += lineagename
    if (options is not None) and ('alignment' in options):
        fn+= '.fasta'
    elif (options is not None) and ('tree' in options):
        fn+= '.nwk'
    else:
        fn+= '.csv'
    return fn


def get_sample_filename(samplename, fmt='fastq'):
    '''Get the filename of a sample'''
    if fmt in ['fastq', 'fasta']:
        fn = get_reads_foldername(fmt=fmt) + samplename + '.' + fmt
    elif fmt == 'igblasted':
        fn = get_igblasted_foldername() + samplename + '.igblasted'
    elif fmt == 'by_sample':
        fn = get_by_sample_foldername() + samplename + '.csv'
    return fn


def get_reads_filenames(fmt='fastq'):
    '''Get the filenames of fastq/fasta files'''
    import glob
    fn = get_reads_foldername(fmt=fmt)
    return glob.glob(fn+'*.'+fmt)


def get_igblastdb_foldername():
    '''Get the foldername of igblastdb, to be set as IGDATA env var'''
    fn = get_root_data_foldername() + 'igblastdb/'
    return fn


def get_germline_db_filename(genename):
    '''Get the germline database filename'''
    fn = 'IGH' + genename + '_imgt.fasta'
    fn = get_igblastdb_foldername() + fn
    return fn


def get_igblast_auxiliary_data_filename():
    '''Get the filename of the human auxiliary file for igblast'''
    fn = 'optional_file/human_gl.aux'
    fn = get_igblastdb_foldername() + fn
    return fn


def get_igblast_execfile(variant='n'):
    '''Get the executable of igblastn/igblastp'''
    if 'glass' in os.getenv('HOME'):
        return '/Users/davidglass/Documents/resources/ncbi-igblast-1.4.0/bin/igblast' + variant
    elif 'fabio' in os.getenv('HOME'):
        return '/usr/bin/igblast' + variant
    else:
        raise ValueError('user not recognized')


def get_muscle_execfile():
    '''Get the executable of igblastn/igblastp'''
    if 'glass' in os.getenv('HOME'):
        raise NotImplementedError('DAVID PLEASE ADD YOUR MUSCLE PATH TO FILENAMES.PY')
    elif 'fabio' in os.getenv('HOME'):
        return '/usr/bin/muscle'
    else:
        raise ValueError('user not recognized')


def get_fasttree_execfile():
    '''Get the executable of igblastn/igblastp'''
    if 'glass' in os.getenv('HOME'):
        raise NotImplementedError('DAVID PLEASE ADD YOUR FASTTREE PATH TO FILENAMES.PY')
    elif 'fabio' in os.getenv('HOME'):
        return '/usr/bin/FastTree'
    else:
        raise ValueError('user not recognized')


def get_parsed_filename(patient_id, time_point, samplename):
    '''Get the parsed igblast filename from by_sample'''
    filename = "patient_%s_time_point_%s_sample_id_%s" % (patient_id, time_point, samplename)
    filename_with_directory = get_sample_filename(filename, fmt='by_sample')
    return filename_with_directory
            

def format_lineage_filename(patient_id, v_gene, j_gene, CDR3_length):
    '''Get the lineage filename given the lineage info'''
    filename = "patient_%s_Vgene_%s_Jgene_%s_CDR3length_%s" % (patient_id, v_gene, j_gene, CDR3_length)
    fn = get_lineage_filename(filename)
    return fn
