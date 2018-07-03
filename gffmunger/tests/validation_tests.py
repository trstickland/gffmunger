import unittest
import os
import subprocess

from gffmunger.GFFMunger import GFFMunger

test_modules_dir  = os.path.dirname(   os.path.realpath( __file__ ) )
data_dir          = os.path.join(      test_modules_dir, 'data' )
test_gff_file     = os.path.join(      data_dir,         'SAMPLE.gff3.gz' )
bad_gff_file      = os.path.join(      data_dir,         'NOT_GFF.gff3' )

gt_test_arg = '-help' # something guaranteed to be OK with any working install of genometools

class GFF3Validator_Tests(unittest.TestCase):
   
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
      
   def test_010_gt_valid_passes(self):
      """check genometools passes a valid GFF3 file"""
      self.assertTrue( self.gffmunger.validate_GFF3(test_gff_file) )

   def test_020_gt_invalid_fails(self):
      """check genometools fails an invalid GFF3 file"""
      self.assertFalse( self.gffmunger.validate_GFF3(bad_gff_file, silent=True) ) # silence validation errors, which are expected
