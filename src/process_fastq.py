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
                       get_germline_db_filename,
                       get_igblast_auxiliary_data_filename,
                       get_igblast_execfile)



# Functions
def run_igblast(filename):
    '''Takes a file from antidengue/data/fasta/ and runs it through igblast
    
    Parameters:
        filename - the name of the fasta
    '''
    import tempfile
    import shutil
    # igblast is terrible and does stuff only in its own folder, so let's make
    # a temp folder for that
    tmp_foldername = tempfile.mkdtemp()

    # Move igblast into the temp
    igblast_tmp = tmp_foldername + 'igblastn'
    shutil.copy(get_igblast_execfile(), igblast_tmp)

    # Move input file into the temp
    filename_tmp = tmp_foldername + 'input.fasta'
    shutil.copy(filename, filename_tmp)

    filename_out_tmp = filename_tmp[:-6] + '_out'

    call([igblast_tmp,
          '-out', filename_out_tmp,
          '-query', filename_tmp,
          '-num_alignments_V=1',
          '-num_alignments_D=1',
          '-num_alignments_J=1',
          '-evalue=1e-20',
          '-germline_db_V', get_germline_db_filename('V'), 
          '-germline_db_D', get_germline_db_filename('D'),
          '-germline_db_J', get_germline_db_filename('J'),
          '-domain_system', 'imgt',
          '-auxiliary_data', get_igblast_auxiliary_data_filename(),
         ])
    
    return filename_out_tmp



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
        self.sample_info = self._get_sample_info_from_csv()


    @staticmethod
    def _get_sample_info_from_csv():
        '''Puts the info from the sample_info csv into a pandas dataframe and returns it
         
        Parameters:
            csv - the csv file
        '''
        
        fn = get_sample_table_filename()
        sample_info = pd.DataFrame().from_csv(fn)
        return sample_info


    def get_reads_filenames(fmt):
        '''Get filenames of reads
        
        Parameters
           - fmt (string): 'fastq' or 'fasta'
        '''
        from filenames import get_reads_filenames as grf
        return grf(fmt=fmt)

    
    def get_sample_filename(samplename, fmt='fastq'):
        '''Get the filename of a sample'''
        from filenames import get_sample_filename as gsf
        return gsf(samplename, fmt=fmt)


    def convert_reads_to_fasta():
        '''Convert reads from FASTQ to FASTA'''
        from filenames import get_reads_foldername as grfn
        from util import mkdirs

        # Make output folder if not existant
        mkdirs(grfn('fasta'))

        # Convert



    def pipeline(self):
        '''Runs a fasta through igblast and then parses the output and properly sorts files
        
        '''
        raw_filenames = self.get_reads_filenames('fasta')
        for fn in raw_filenames:
            blast_out = run_igblast(filename)
            run_parse(filename, blast_out, self.sample_info)




# Script
if __name__ == '__main__':

    lm = LineageMaker()
    #lm.pipeline()

