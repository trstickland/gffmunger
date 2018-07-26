import unittest
import os
import gffutils
import warnings

from gffmunger.GFFMunger import GFFMunger

test_modules_dir        = os.path.dirname(   os.path.realpath( __file__ ) )
data_dir                = os.path.join(      test_modules_dir, 'data' )
test_gff_file           = os.path.join(      data_dir,         'SMALL_SAMPLE.gff3.gz' )                  # must be valid GFF3 with *no* FASTA
test_gff_and_fasta_file = os.path.join(      data_dir,         'SMALL_SAMPLE_INCL_FASTA.gff3.gz' ) # must be valid GFF3 with FASTA
test_gff_db_file        = os.path.join(      data_dir,         'gffutils_test.db' )

#sample_gff_gene_id      = 'PF3D7_0100100' # use with SAMPLE.gff3
sample_gff_gene_id      = '13J3.01' # use with SMALL_SAMPLE.gff3
sample_gff_featuretypes = ['gene', 'mRNA', 'CDS', 'polypeptide'] 

expected_db_class          = gffutils.FeatureDB
expected_feature_class     = gffutils.Feature
expected_attributes_class  = gffutils.attributes.Attributes

class GFF_Tests(unittest.TestCase):
   
   test_gff_db    = None
   db_available   = False
   
   @classmethod
   def setUpClass(self):
      with warnings.catch_warnings():
         warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper",  ResourceWarning,           "gffutils", 133 )
         warnings.filterwarnings("ignore", "generator '_FileIterator\.",         PendingDeprecationWarning, "gffutils", 186 )
         warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper",  ResourceWarning,           "gffutils", 668 )
         self.test_gff_db = gffutils.create_db( test_gff_file,
                                                dbfn                    = test_gff_db_file,
                                                force                   = True,     # overwrite previous testing db file
                                                merge_strategy          = 'error',
                                                keep_order              = False,    # True doesn't appear to maintain attribute order :-/  (and turning this off may be faster)
                                                sort_attribute_values   = False
                                                )
         if(self.test_gff_db is not None and isinstance(self.test_gff_db, expected_db_class)):
            self.db_available = True
      self.gffmunger = GFFMunger(None)

                        
   def test_000_db_created(self):
      """test that gffutils.create_db was able to create a database from the sample GFF file"""
      #self.assertIsNotNone(self.test_gff_db)
      #self.assertIsInstance(self.test_gff_db, expected_db_class)
      self.assertTrue(self.db_available)
      
   #@unittest.skipUnless(db_available, "no db available")
   def test_010_gffutils_find_gene(self):
      """test that gffutils can find a sample gene identifier from the GFF"""
      if (not self.db_available):
         self.skipTest('no db available')
      this_gene = self.test_gff_db[sample_gff_gene_id]
      self.assertIsNotNone(this_gene)
      self.assertIsInstance(this_gene, expected_feature_class)


   def test_020_gffutils_find_features(self):
      """test that gffutils can find features in GFF"""
      if (not self.db_available):
         self.skipTest('no db available')
      for this_featuretype in sample_gff_featuretypes:
         for this_feature in self.test_gff_db.children(sample_gff_gene_id, featuretype=this_featuretype):
            self.assertIsNotNone(this_feature)
            self.assertIsInstance(this_feature, expected_feature_class)
            # print(this_feature)

   def test_030_gffutils_find_attributes(self):
      """test that gffutils can find attributes in GFF"""
      if (not self.db_available):
         self.skipTest('no db available')
      for this_featuretype in sample_gff_featuretypes:
         for this_feature in self.test_gff_db.children(sample_gff_gene_id, featuretype=this_featuretype):
            #print(this_feature)
            self.assertIsNotNone(this_feature)
            self.assertIsInstance(this_feature, expected_feature_class)
            these_attributes = this_feature.attributes
            self.assertIsNotNone(these_attributes)
            self.assertIsInstance(these_attributes, expected_attributes_class)
            for this_item in sorted(these_attributes.items()):
               self.assertIsNotNone(this_item)
               #print('{this_item[0]}: {this_item[1]}'.format(this_item=this_item))
            
   def test_040_gff_components_no_fasta(self):
      """test separation of GFF3 file without FASTA, into metadata and features """
      self.gffmunger.extract_GFF3_components(test_gff_file)    
      self.assertIsNotNone(self.gffmunger.input_metadata)
      if self.gffmunger.read_features_to_buffer:
         self.assertIsNotNone(self.gffmunger.input_features)
      else:
         self.assertIsNone(self.gffmunger.input_features)
      self.assertIsNone(self.gffmunger.input_fasta)

   def test_045_gff_components_with_fasta(self):
      """test separation of GFF3 file with FASTA, into metadata, features and FASTA data"""
      self.gffmunger.extract_GFF3_components(test_gff_and_fasta_file)    
      self.assertIsNotNone(self.gffmunger.input_fasta)
