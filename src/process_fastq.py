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
from collections import defaultdict, Counter


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

        # Make output folder if not extant
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

        # Make output folder if not extant
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
        '''

        from filenames import get_by_sample_foldername as gfn
        from util import mkdirs
        from parse_igblast import main as parse_function

        if samplenames is None:
            samplenames = self.sample_table.index

        # Make output folder if not extant
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

            filename_out = self.get_igblast_parsed_filename(samplename)

            # Call the function from the parser module
            parse_function(filename_fa,
                           filename_in,
                           output_filename=filename_out,
                           len_C_cutoff_min=0,
                           len_C_cutoff_max=10000,
                           log_V_Evalue_cutoff=-5.0,
                           log_J_Evalue_cutoff=-5.0)


    def add_sequence_quality_by_sample(self, samplenames=None):
        '''Renames the parsed igblast file by patient, timepoint, and ID. Adds the complete nucleotide
        sequence and quality scores for each sequence. Prepends the file with a header.
        '''

        from filenames import get_reads_foldername

        header = '\t'.join(['#sequence_id', 'abundance',
                            'V-gene', 'D-gene', 'J-gene',
                            'V_E-value', 'D_E-value', 'J_E-value',
                            'FR3_seq', 'CDR3_seq', 'FR4_seq', 'const_seq',
                            'len_FR3', 'len_CDR3', 'len_FR4', 'len_const',
                            'sequence_boundary_indices:FR1_CDR1_FR2_CDR2_FR3_CDR3_FR4',
                            'len_boundary', 'stop_codon_present',
                            'productive_sequence', 'AA_seq_whole_read',
                            'mutation_positions', 'germline_bases', 'derived_bases',
                            'mutation_density', 'V_germline_identity',
                            'leader_seq', 'reads_per_molecule', 'primer_isotype',
                            'sequence_reversed', 'seq', 'quality'])

        if samplenames is None:
            samplenames = self.sample_table.index

        for samplename in samplenames:
            fn_fastq = self.get_sample_filename(samplename, fmt='fastq')
            fn_by_sample = self.get_igblast_parsed_filename(samplename)

            if not os.path.isfile(fn_fastq):
                print(samplename+' fastq file not found:', fn_fastq)
                continue
            
            if not os.path.isfile(fn_by_sample):
                print(samplename+' parsed igblast file not found:', fn_by_sample)
                continue
            
            with open(fn_fastq, 'r') as f_fq:
                fastq_lines = f_fq.readlines()

            with open(fn_by_sample, 'r') as f_bs: 
                parsed_lines = [header] + f_bs.readlines()

            for i in range(1, len(parsed_lines)):
                parsed_lines[i] = (parsed_lines[i].rstrip('\n') + '\t' +
                                   fastq_lines[(i*4)-3].rstrip('\n') + '\t' +
                                   fastq_lines[(i*4)-1].rstrip('\n'))
            
            with open (fn_by_sample, 'w') as fo:
                fo.write('\n'.join(parsed_lines))


    def get_igblast_parsed_filename(self, samplename):
        '''Returns the by_sample filename for a given sample.
        '''
        from filenames import get_parsed_filename as gfn
        patient_id = str(self.sample_table.loc[samplename][0])
        time_point = str(self.sample_table.loc[samplename][4])
        filename = gfn(patient_id, time_point, samplename)
        return filename


    def get_lineage_name(self, patient_id, lineage_info):
        '''Returns the lineage filename for a given lineage.
        '''
        from filenames import format_lineage_filename as gfn
        v_gene = lineage_info[0]
        j_gene = lineage_info[1]
        CDR3_length = lineage_info[2]
        filename = gfn(patient_id, v_gene, j_gene, CDR3_length)
        return filename


    def sort_by_lineage(self, samplenames=None):
        '''Sorts sequences into files based on lineage. This is defined by patient, V gene, 
        J gene, and CDR3 length.
        '''

        def parse_line(line):
            '''Takes a line from the parsed Igblast file and returns an ID string and an array 
            containing V gene, J gene, and CDR3 length.
            '''
            line = line.rstrip('\n')
            if not line:
                return None, None

            split_line = line.split()
            sample_id = ''
            lineage_data = [None] * 3

            sample_id = split_line[0]
            lineage_data[0] = split_line[2]
            lineage_data[1] = split_line[4]
            lineage_data[2] = str(len(split_line[9]))
            
            return sample_id, lineage_data

        from util import mkdirs
        from filenames import get_lineages_foldername as gfn

        if samplenames is None:
            samplenames = self.sample_table.index

        # Make output folder if not extant
        mkdirs(gfn())

        all_samples = defaultdict(dict)
        for samplename in samplenames:
            patient_id = str(self.sample_table.loc[samplename, 'Specimen ID'])
            fn_by_sample = self.get_igblast_parsed_filename(samplename)

            if not os.path.isfile(fn_by_sample):
                print(samplename+' parsed igblast file not found:', fn_by_sample)
                continue

            with open(fn_by_sample, 'r') as f_bs:
                for line in f_bs:
                    if line.startswith('#'):
                        continue
                    sample_id, lineage_data = parse_line(line)
                    if lineage_data is not None:
                        all_samples[patient_id][sample_id] = lineage_data
        

        header = '\t'.join(['#sequence_id', 'abundance',
                    'V-gene', 'D-gene', 'J-gene',
                    'V_E-value', 'D_E-value', 'J_E-value',
                    'FR3_seq', 'CDR3_seq', 'FR4_seq', 'const_seq',
                    'len_FR3', 'len_CDR3', 'len_FR4', 'len_const',
                    'sequence_boundary_indices:FR1_CDR1_FR2_CDR2_FR3_CDR3_FR4',
                    'len_boundary', 'stop_codon_present',
                    'productive_sequence', 'AA_seq_whole_read',
                    'mutation_positions', 'germline_bases', 'derived_bases',
                    'mutation_density', 'V_germline_identity',
                    'leader_seq', 'reads_per_molecule', 'primer_isotype',
                    'sequence_reversed', 'seq', 'quality'])
        lineage_counts = Counter()
        for patient_id, patient_samples in all_samples.items():
            for sample_id, lineage_data in patient_samples.items():
                samplename = sample_id.split('.')[0]
                fn_by_sample = self.get_igblast_parsed_filename(samplename)
                fn_lineage = self.get_lineage_name(patient_id, lineage_data)

                if not os.path.isfile(fn_by_sample):
                    print(samplename+' parsed igblast file not found:', fn_by_sample)
                    continue

                lineage_counts[os.path.basename(fn_lineage)] += 1
                if lineage_counts[os.path.basename(fn_lineage)] == 1:
                    with open(fn_lineage, 'w') as f_l:
                        f_l.write(header+'\n')

                with open(fn_by_sample, 'r') as f_bs, open(fn_lineage, 'a') as f_l:
                    for line in f_bs:
                        if sample_id in line:
                            f_l.write(line)
                            break  


    def pipeline(self, samplenames=None, steps=None):
        '''Runs a fasta through igblast, parses the output and sorts files'''
        if (steps is None) or ('convert_fastq' in steps):
            self.convert_reads_to_fasta(samplenames=samplenames)
        
        if (steps is None) or ('igblast' in steps):
            self.run_igblast(samplenames=samplenames)

        if (steps is None) or ('parse_igblast' in steps):
            self.parse_igblast(samplenames=samplenames)

        if (steps is None) or ('add_sequence_quality' in steps):
            self.add_sequence_quality_by_sample(samplenames=samplenames)

        if (steps is None) or ('sort_by_lineage' in steps):
            self.sort_by_lineage(samplenames=samplenames)


    # FIXME: To be deleted after testing!!
    def temporary_get_samples(self, samplenumber=-1, samplenames=None):
        '''Grabs the first 'samplenumber' samples to run through the pipeline'''
        
        if samplenames is not None:
            return samplenames

        if samplenumber == -1:
            return samplenames

        samplenames = []
        for i in range(samplenumber):
            samplenames.append(self.sample_table.index.iloc[i])
        return samplenames


# Script
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Adds single sample processing functionality')
    parser.add_argument('-s', '--samplenames', nargs='+', default=None,
                        help='a list of sample names to run through the process')
    parser.add_argument('-n', '--samplenumber', type=int, default=-1,
                         help='an integer directing the pipeline to be run on the first X samples')
    parser.add_argument('--steps', nargs='+', default=None,
                        help='pipeline steps to perform')
    args = parser.parse_args()

    lm = LineageMaker()
    samplenames = lm.temporary_get_samples(args.samplenumber, args.samplenames)
    lm.pipeline(samplenames=samplenames, steps=args.steps)
