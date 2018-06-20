import unittest
import os
import gffutils
import warnings

test_modules_dir  = os.path.dirname(   os.path.realpath( __file__ ) )
data_dir          = os.path.join(      test_modules_dir, 'data' )
test_gff_file     = os.path.join(      data_dir,         'SAMPLE.gff3' )
test_gff_db_file  = os.path.join(      data_dir,         'gffutils_test.db' )

sample_gff_gene_id   = 'PF3D7_0100100'

expected_db_class    = gffutils.FeatureDB
expected_gene_class  = gffutils.Feature

class GFF_Tests(unittest.TestCase):
   
   test_gff_db    = None
   
   def setUp(self):
      with warnings.catch_warnings():
         warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper",  ResourceWarning,           "gffutils", 133 )
         warnings.filterwarnings("ignore", "generator '_FileIterator\.",         PendingDeprecationWarning, "gffutils", 186 )
         warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper",  ResourceWarning,           "gffutils", 668 )
         self.test_gff_db = gffutils.create_db( test_gff_file,
                                                dbfn                    = test_gff_db_file,
                                                force                   = True,     # overwrite previous testing db
                                                keep_order              = True,
                                                merge_strategy          = 'merge',
                                                sort_attribute_values   = True
                                                )
         
   def tearDown(self):
      #self.test_gff_db.close()
      del self.test_gff_db
      
   def test_000_db_created(self):
      """test that gffutils.create_db was able to create a database from the sample GFF file"""
      self.assertIsNotNone(self.test_gff_db)
      self.assertIsInstance(self.test_gff_db, expected_db_class)

   #@unittest.skipUnless(test_gff_db, "no db available")
   def test_010_gffutils_find_gene(self):
      """test that gffutils can find a sample gene identifier from the GFF"""
      gene = self.test_gff_db[sample_gff_gene_id]
      self.assertIsNotNone(gene)
      self.assertIsInstance(gene, expected_gene_class)
      
      
