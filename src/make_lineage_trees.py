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


class TreeMaker():
    '''Make phylogenetic trees from the antibody lineages'''
    def __init__(lineage_options=None):
        self.sample_table = SampleTable.from_sampletable()
        self.lineage_table = LineageTable.from_lineagetable()
        self.lineage_options = lineage_options


    def get_lineage_filename(self, lineagename):
        '''Get the filenames of the lineages'''
        from filenames import get_lineage_filename as gsf
        return gsf(samplename, options=self.lineage_options)


    def align_lineage(self, lineagename, region='CDR3'):
        '''Align sequences from a lineage in perspective of making a tree'''
        pass



# Script
if __name__ == '__main__':

    lm = TreeMaker()

