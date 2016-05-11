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




# Classes
class SampleTable(pd.DataFrame):
    '''Table of antibody sequence samples'''

    @property
    def _constructor(self):
        return SampleTable


    @property
    def _constructor_sliced(self):
        return pd.Series


    @property
    def _constructor_expanddim(self):
        raise NotImplementedError


    @classmethod
    def from_sampletable(cls):
        '''Puts the info from the sample_info csv'''
        from filenames import get_sample_table_filename
        
        fn = get_sample_table_filename()
        return cls.from_csv(fn)


class LineageMaker():
    '''Make Antibody lineages from patient sequencing data'''
    def __init__(self):
        self.sample_table = SampleTable.from_sampletable()


    @staticmethod
    def get_reads_filenames(fmt):
        '''Get filenames of reads
        
        Parameters
           - fmt (string): 'fastq' or 'fasta'
        '''
        from filenames import get_reads_filenames as grf
        return grf(fmt=fmt)

    
    @staticmethod
    def get_sample_filename(samplename, fmt='fastq'):
        '''Get the filename of a sample'''
        from filenames import get_sample_filename as gsf
        return gsf(samplename, fmt=fmt)


    def convert_reads_to_fasta(self, samplenames=None):
        '''Convert reads from FASTQ to FASTA'''
        from filenames import get_reads_foldername as gfn
        from util import mkdirs

        if samplenames is None:
            samplenames = self.sample_table.index

        # Make output folder if not existant
        mkdirs(gfn('fasta'))

        # Convert
        from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI
        for samplename in samplenames:
            fn_fastq = self.get_sample_filename(samplename, fmt='fastq')
            fn_fasta = self.get_sample_filename(samplename, fmt='fasta')
            if not os.path.isfile(fn_fastq):
                print(samplename+' fastq file not found:', fn_fastq)
                continue

            with open(fn_fastq, 'r') as f_fq, open(fn_fasta, 'w') as f_fa:
                for (label, seq, qual) in FGI(f_fq):
                    f_fa.write('>'+label+'\n')
                    f_fa.write(seq+'\n')


    def run_igblast(self, samplenames=None):
        '''Takes a file from antidengue/data/fasta/ and runs it through igblast
        
        Parameters:
            filename - the name of the fasta
        '''

        def get_environment_igblast():
            '''Set environment for igblast call
            
            igblast uses the IGDATA environment variable to find the
            "internal_data" folder with some human V germline nonsense
            '''
            from filenames import get_igblastdb_foldername
            env = dict(os.environ)
            if 'IGDATA' not in env:
                env['IGDATA'] = get_igblastdb_foldername()[:-1]
            return env

        from filenames import get_igblasted_foldername as gfn
        from filenames import (get_igblast_execfile,
                               get_germline_db_filename,
                               get_igblast_auxiliary_data_filename)
        from util import mkdirs

        if samplenames is None:
            samplenames = self.sample_table.index

        # Make output folder if not existant
        mkdirs(gfn())
    
        for samplename in samplenames:
            filename_in = self.get_sample_filename(samplename, fmt='fasta')
            if not os.path.isfile(filename_in):
                print(samplename+' fasta file not found:', filename_in)
                continue

            filename_out = self.get_sample_filename(samplename, fmt='igblasted')
            cll = [get_igblast_execfile(),
                   '-out', filename_out,
                   '-query', filename_in,
                   '-num_alignments_V=1',
                   '-num_alignments_D=1',
                   '-num_alignments_J=1',
                   '-evalue=1e-20',
                   '-organism', 'human',
                   '-ig_seqtype', 'Ig',
                   '-num_threads', '2',
                   '-germline_db_V', get_germline_db_filename('V'), 
                   '-germline_db_D', get_germline_db_filename('D'),
                   '-germline_db_J', get_germline_db_filename('J'),
                   '-domain_system', 'imgt',
                   '-auxiliary_data', get_igblast_auxiliary_data_filename(),
                  ]
            print(' '.join(cll))
            call(cll, env=get_environment_igblast())
    
    
    def parse_igblast(self, samplenames=None):
        '''Runs parse_igblast.py in the igblast_dump folder and sends the parsedfile to the by_patient directory
        
        Parameters:
            fasta - the fasta file
            blast_out - the output of igblast
            dataframe - the pandas dataframe containing patient info
        '''

        from filenames import get_igblastparsed_foldername as gfn
        from util import mkdirs
        from parse_igblast import main as parse_function

        if samplenames is None:
            samplenames = self.sample_table.index

        # Make output folder if not existant
        mkdirs(gfn())

        for samplename in samplenames:
            filename_fa = self.get_sample_filename(samplename, fmt='fasta')
            if not os.path.isfile(filename_fa):
                print(samplename+' fasta file not found:', filename_fa)
                continue
            filename_in = self.get_sample_filename(samplename, fmt='igblasted')
            if not os.path.isfile(filename_in):
                print(samplename+' igblasted file not found:', filename_in)
                continue

            filename_out = self.get_sample_filename(samplename, fmt='igblastparsed')

            # Call the function from the parser module
            parse_function(filename_fa,
                           filename_in,
                           output_filename=filename_out,
                           len_C_cutoff_min=0,
                           len_C_cutoff_max=10000,
                           log_V_Evalue_cutoff=-5.0,
                           log_J_Evalue_cutoff=-5.0)



    def pipeline(self):
        '''Runs a fasta through igblast, parses the output and sorts files'''
        self.convert_reads_to_fasta()
        self.run_igblast()
        self.parse_igblast()





# Script
if __name__ == '__main__':

    lm = LineageMaker()
    lm.pipeline()

