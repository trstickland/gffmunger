import argparse
import gffutils
import logging
import os
import pyfaidx
import unittest
import uuid
import warnings

from gffmunger.GFFMunger import GFFMunger

test_modules_dir  = os.path.dirname(   os.path.realpath( __file__ ) )
data_dir          = os.path.join(      test_modules_dir, 'data' )
test_gff_file     = os.path.join(      data_dir,         'SMALL_SAMPLE_INCL_FASTA.gff3.gz' )
test_gff_no_fasta = os.path.join(      data_dir,         'SMALL_SAMPLE.gff3.gz' )
test_fasta_file   = os.path.join(      data_dir,         'SMALL_SAMPLE.fasta' )

# this file has broken mRNA-polypeptide relations, including a simple broken reference (13J3.17:pep),
# which produces two errors (mRNA with no polypeptide; and a polypeptide with no mRNA);
# multiple polypeptides (H25N7.05.01:pep and H25N7.05.02:pep) referring to one mRNA;
# polypeptide missing a 'Dervies_from' attribute (H25N7.09:pep)
broken_gff_file   = os.path.join(      data_dir,         'SMALL_SAMPLE_BROKEN_RELATIONS.gff3.gz' )

expected_num_input_lines      = 2082
expected_num_metadata_lines   = 43

class IO_Tests(unittest.TestCase):
   
   @classmethod
   def setUpClass(self):
      self.output_file  = __file__+'.'+uuid.uuid4().hex+'.gff3'

   @classmethod
   def tearDownClass(self):
      if self.output_file and os.path.exists(self.output_file):
         os.remove(self.output_file)
     
   def test_000_gff3_with_fasta_io(self):
      """check functions for reading and validating a GFF3 file including FASTA, and writing it as output"""
      gffmunger = GFFMunger( None )
      gffmunger.input_file_arg   = test_gff_file
      gffmunger.output_file      = self.output_file
      self.assertEqual(test_gff_file,  gffmunger.get_gff3_source())
      self.assertTrue(gffmunger.validate_GFF3(test_gff_file))
      self.assertEqual(gffmunger.gffutils_db_filename, gffmunger.import_gff3())
      self.assertEqual(gffmunger.gffutils_db_filename, gffmunger.import_gff3(test_gff_file))
      self.assertEqual(expected_num_input_lines,            gffmunger.extract_GFF3_components())
      self.assertEqual(expected_num_input_lines,            gffmunger.extract_GFF3_components(test_gff_file))
      with warnings.catch_warnings():
         warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper", ResourceWarning, "gffutils", 668 )
         try:
            gffmunger.move_annotations()
         except:
            self.fail("Failed to process valid GFF3 file "+test_gff_file)
      warnings.resetwarnings()
      self.assertTrue(gffmunger.export_gff3())
      os.remove(gffmunger.output_file)
      gffmunger.clean_up()

   def test_020_gff3_i0(self):
      """check functions for retreiving input filename, and reading and validating a GFF3 file and separate FASTA file"""
      newmunger = GFFMunger( None )
      newmunger.input_file_arg   = test_gff_no_fasta
      newmunger.fasta_file_arg   = test_fasta_file
      newmunger.output_file      = self.output_file
      self.assertEqual(test_gff_no_fasta,  newmunger.get_gff3_source())
      self.assertTrue(newmunger.validate_GFF3(test_gff_no_fasta))
      self.assertTrue(newmunger.validate_FASTA(test_fasta_file))
      self.assertEqual(newmunger.gffutils_db_filename,   newmunger.import_gff3())
      self.assertEqual(newmunger.gffutils_db_filename,   newmunger.import_gff3(test_gff_no_fasta))
      # when FASTA is in separate file, GFFMunger.extract_GFF3_components() should read only metadata
      self.assertEqual(expected_num_metadata_lines,      newmunger.extract_GFF3_components())
      self.assertEqual(expected_num_metadata_lines,      newmunger.extract_GFF3_components(test_gff_no_fasta))
      with warnings.catch_warnings():
         warnings.filterwarnings("ignore", "unclosed file", ResourceWarning, "gffutils", 0 )
         self.assertIsInstance(newmunger.import_fasta(), pyfaidx.Fasta)
         self.assertIsInstance(newmunger.import_fasta(test_fasta_file), pyfaidx.Fasta)
         try:
            newmunger.move_annotations()
         except:
            self.fail("Failed to process valid GFF3 file "+test_gff_no_fasta)
      warnings.resetwarnings()
      self.assertTrue(newmunger.export_gff3())
      os.remove(newmunger.output_file)
      newmunger.clean_up()


   def test_050_gff_error_handling(self):
      """checks handling of non-fatal errors encountered in GFF"""
      yet_another_munger = GFFMunger( None )
      yet_another_munger.input_file_arg = broken_gff_file
      self.assertTrue(yet_another_munger.validate_GFF3(broken_gff_file)) # should pass validation despire bad relations
      with warnings.catch_warnings():
         warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper", ResourceWarning, "gffutils", 668 )
         yet_another_munger.import_gff3(broken_gff_file)
         try:
            oldloglevel = yet_another_munger.logger.level
            yet_another_munger.logger.setLevel(logging.CRITICAL)
            yet_another_munger.move_annotations()
            yet_another_munger.logger.setLevel(oldloglevel)
         except AssertionError:
            self.fail("AssertionError should not be raised by GFFMunger.move_annotations() when processing annotations in "+broken_gff_file)
      warnings.resetwarnings()
      yet_another_munger.clean_up()
