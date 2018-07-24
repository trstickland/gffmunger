import unittest
import os
import subprocess

from gffmunger.GFFMunger import GFFMunger

test_modules_dir  = os.path.dirname(   os.path.realpath( __file__ ) )
data_dir          = os.path.join(      test_modules_dir, 'data' )
test_gff_file     = os.path.join(      data_dir,         'SAMPLE.gff3.gz' )
bad_gff_file      = os.path.join(      data_dir,         'NOT_GFF.gff3' )
test_fasta_file   = os.path.join(      data_dir,         'SAMPLE.fasta.gz' )
bad_fasta_file    = os.path.join(      data_dir,         'NOT_FASTA.fasta' )

gt_test_arg = '-help' # something guaranteed to be OK with any working install of genometools

class Validator_Tests(unittest.TestCase):
   
   @classmethod
   def setUpClass(self):
      self.gffmunger = GFFMunger(None)
          
   def test_000_gt_installed(self):
      """check genometools is installed and exceuteable"""
      cp = subprocess.run([self.gffmunger.gt_path, gt_test_arg], timeout=3, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
      success = 0 == cp.returncode
      if not success:
         print(cp.stdout)
      self.assertTrue(success)
      
   def test_010_gt_valid_gff_passes(self):
      """check GFFMunger.validate_GFF3 passes a valid GFF3 file"""
      self.assertTrue( self.gffmunger.validate_GFF3(test_gff_file) )

   def test_020_gt_invalid_gff_fails(self):
      """check GFFMunger.validate_GFF3 fails an invalid GFF3 file"""
      self.assertFalse( self.gffmunger.validate_GFF3(bad_gff_file, silent=True) ) # silence validation errors, which are expected

   def test_030_gt_valid_fasta_passes(self):
      """check GFFMunger.validate_FASTA passes a valid FASTA file"""
      self.assertTrue( self.gffmunger.validate_FASTA(test_fasta_file) )

   def test_040_gt_invalid_fasta_fails(self):
      """check GFFMunger.validate_FASTA fails an invalid FASTA file"""
      self.assertFalse( self.gffmunger.validate_FASTA(bad_fasta_file, silent=True) ) # silence validation errors, which are expected
