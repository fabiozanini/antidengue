# vim: fdm=indent
'''
author:     David Glass, Fabio Zanini
date:       09/05/16
content:    Process Antibody nucleotide sequences into lineages.
'''
# Modules
from Bio import SeqIO    
import numpy as np
import pandas as pd    
import matplotlib.pyplot as plt
import os
import sys
from subprocess import call

from filenames import (get_sample_table_filename,
                       get_reads_filenames)



# Functions
def run_igblast(filename):
    '''Takes a file from antidengue/data/fasta/ and runs it through igblast
    
    Parameters:
        filename - the name of the fasta
    '''
    if filename.endswith('.fasta'):
        outfile = filename[:-6] + '_out'
    else:
        outfile = filename + '_out'
    
    fasta_path = '/Users/davidglass/antidengue/data/fasta/'
    ncbi_path = '/Users/davidglass/Documents/resources/ncbi-igblast-1.4.0/'
    blast = '/Users/davidglass/Documents/resources/ncbi-igblast-1.4.0/bin/igblastn'
    igblast_dump = '/Users/davidglass/antidengue/data/ig_parse_dump/' + outfile[:-4]

    call(['cp', fasta_path + filename, ncbi_path])
    os.chdir(ncbi_path)
    call([blast, '-out', outfile, '-query', ncbi_path+filename, '-num_alignments_V=1', '-num_alignments_D=1',
          '-num_alignments_J=1', '-evalue=1e-20', '-germline_db_V', ncbi_path+'database/IGHV_imgt.fasta', 
          '-germline_db_D', ncbi_path+'database/IGHD_imgt.fasta', '-germline_db_J',
          ncbi_path+'database/IGHJ_imgt.fasta', '-domain_system', 'imgt', '-auxiliary_data',
          ncbi_path+'optional_file/human_gl.aux'])
    call(['mkdir', igblast_dump])
    call(['mv', filename, igblast_dump])
    call(['mv', outfile, igblast_dump])
    
    return outfile


def get_sample_info_from_csv(csv):
    '''Puts the info from the sample_info csv into a pandas dataframe and returns it
     
    Parameters:
        csv - the csv file
    '''
    
    fn = get_sample_table_filename()
    sample_info = pd.DataFrame().from_csv(fn)
    return sample_info


def run_parse(fasta, blast_out, sample_info):
    '''Runs parse_igblast.py in the igblast_dump folder and sends the parsedfile to the by_patient directory
    
    Parameters:
        fasta - the fasta file
        blast_out - the output of igblast
        dataframe - the pandas dataframe containing patient info
    '''
    accession_num = blast_out[:-4]
    igblast_dump = '/Users/davidglass/antidengue/data/ig_parse_dump/' + accession_num
    destination_directory = '/Users/davidglass/antidengue/data/by_sample'
    patient_id = ""
    time_point = ""
    
    call(['cp', '/Users/davidglass/antidengue/src/parse_igblast.py', igblast_dump])
    os.chdir(igblast_dump)
    call(['python', 'parse_igblast.py', fasta, blast_out])
    
    # find the ID from the sample_info DataFrame
    patient_id = str(sample_info.loc[accession_num][0])
    time_point = str(sample_info.loc[accession_num][4])
    new_filename = patient_id + "_" + time_point + "_" + accession_num
    
    call(['mv', 'parsed_igblast.txt', new_filename])
    call(['mv', new_filename, destination_directory])
    call(['rm', fasta])
    call(['rm', 'parse_igblast.py'])
    


# Classes
class LineageMaker():
    '''Make Antibody lineages from patient sequencing data'''
    def __init__(self):
        self.sample_info = get_sample_info_from_csv()


    def pipeline(self, filename):
        '''Runs a fasta through igblast and then parses the output and properly sorts files
        
        Parameters:
            filename - the name of the fasta
        '''
    
        blast_out = run_igblast(filename)
        run_parse(filename, blast_out, self.sample_info)




# Script
if __name__ == '__main__':

    lm = LineageMaker()

    raw_filenames = get_reads_filenames('fasta')
    for fn in raw_filenames:
        lm.pipeline(fn)
