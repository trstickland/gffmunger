import unittest
import os
import gzip
import warnings

from Bio.UniProt.GOA import gafiterator


test_modules_dir     = os.path.dirname(   os.path.realpath( __file__ ) )
data_dir             = os.path.join(      test_modules_dir, 'data' )
test_gaf_file        = os.path.join(      data_dir,         'SAMPLE.gaf.gz' )    # must be valid GAF

sample_gaf_id   = 'PF3D7_1033000.1'  # use with SAMPLE.gaf


class GAF_Tests(unittest.TestCase):

   def test_000_read_gaf_file(self):
      with gzip.open(test_gaf_file, 'rt') as fp:
         found_sample_gene = False
         for annotation in gafiterator(fp):
            self.assertIsNotNone(annotation)
            if annotation['DB_Object_ID'] == sample_gaf_id:
               found_sample_gene = True
         self.assertTrue(found_sample_gene)
