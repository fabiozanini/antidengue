# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/05/16
content:    Make phylogenetic trees of lineages
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd

from process_fastq import SampleTable



# Classes
class LineageTable(pd.DataFrame):
    '''Table of antibody sequence samples'''

    @property
    def _constructor(self):
        return LineageTable


    @property
    def _constructor_sliced(self):
        return pd.Series


    @property
    def _constructor_expanddim(self):
        raise NotImplementedError


    @classmethod
    def from_lineagetable(cls):
        '''Puts the info from the lineage csv'''
        from filenames import get_lineage_table_filename
        
        fn = get_lineage_table_filename()
        return cls.from_csv(fn)


    @classmethod
    def from_raw_filenames(cls, options=None):
        '''Creates a CSV table from the raw filenames'''
        from filenames import get_lineages_foldername as gfn

        table = []
        foldername = gfn(options=options)
        for fn in os.listdir(foldername):
            if not fn.endswith('.csv'):
                continue
            line = fn[:-4]
            fields = line.split('_')
            names = fields[::2]
            values = fields[1::2]
            datum = dict(zip(names, values))
            datum['lineage_id'] = line
            for key, value in datum.items():
                if 'length' in key:
                    datum[key] = int(value)
            table.append(datum)

        return cls(table).set_index('lineage_id')


class TreeMaker():
    '''Make phylogenetic trees from the antibody lineages'''
    def __init__(self, lineage_options=None):
        if lineage_options is not None:
            self.lineage_options = lineage_options
        else:
            self.lineage_options = {}

        self.sample_table = SampleTable.from_sampletable()
        try:
            self.lineage_table = LineageTable.from_lineagetable()
        except(OSError, IOError):
            self.lineage_table = LineageTable.from_raw_filenames(options=lineage_options)


    def _merge_lineage_options(self, options):
        '''Merge list/dict lineage options with instance-wide options'''
        opts = dict(self.lineage_options)
        if options is None:
            return opts
        elif isinstance(options, dict):
            opts.update(options)
            return opts
        else:
            for opt in options:
                opts[opt] = None
            return opts


    def get_lineage_filename(self, lineagename, options=None):
        '''Get the filenames of the lineages'''
        from filenames import get_lineage_filename as gsf
        return gsf(lineagename, options=self._merge_lineage_options(options))


    def align_lineages(self, lineagenames=None, region='CDR3'):
        '''Align sequences from a lineage in perspective of making a tree'''
        import subprocess as sp
        import tempfile

        from filenames import get_lineages_foldername as gfn
        from filenames import get_muscle_execfile
        from util import mkdirs

        if lineagenames is None:
            lineagenames = self.lineage_table.index

        # Make output folder if not extant
        mkdirs(gfn(options=self._merge_lineage_options(['alignment'])))

        for lineagename in lineagenames:
            fn_lineage = self.get_lineage_filename(lineagename)
            fn_aligned = self.get_lineage_filename(lineagename, options=['alignment'])
            if not os.path.isfile(fn_lineage):
                print(lineagename+' lineage file not found:', fn_lineage)
                continue

            data = pd.read_csv(fn_lineage, usecols=['#sequence_id', 'CDR3_seq'],
                               index_col='#sequence_id',
                               sep='\t')['CDR3_seq']
            data.index.name = 'sequence_id'

            # Copy sequence data into temporary fasta
            tmp_foldername = tempfile.mkdtemp()
            tmp_filename = tmp_foldername + '/'+lineagename+'.fasta'
            with open(tmp_filename, 'w') as f_tmp:
                for sequence_id, seq in data.iteritems():
                    f_tmp.write('>'+sequence_id+'\n')
                    f_tmp.write(seq+'\n')

            call = [get_muscle_execfile(),
                    '-in', tmp_filename,
                    '-out', fn_aligned,
                    '-quiet',
                    '-diags']
            print(' '.join(call))
            sp.call(call)


    def pipeline(self, lineagenames=None):
        '''Run the pipeline'''
        self.align_lineages(lineagenames=lineagenames)



# Script
if __name__ == '__main__':

    tm = TreeMaker()
    tm.pipeline()

