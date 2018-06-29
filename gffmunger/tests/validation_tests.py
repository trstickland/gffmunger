import unittest
import os
import subprocess

test_modules_dir  = os.path.dirname(   os.path.realpath( __file__ ) )
data_dir          = os.path.join(      test_modules_dir, 'data' )
test_gff_file     = os.path.join(      data_dir,         'SAMPLE.gff3.gz' )
bad_gff_file      = os.path.join(      data_dir,         'NOT_GFF.gff3' )

gt_cmd               = '/usr/bin/gt'
test_arg             = '-help' # something guaranteed to be OK with any working install
gff3_validator_tool  = 'gff3validator'

class GFF3Validator_Tests(unittest.TestCase):
   
   #@classmethod
   #def setUpClass(self):
      ## whatever...
                        
   def test_000_gt_installed(self):
      """check genometools is installed and exceuteable"""
      cp = subprocess.run([gt_cmd, test_arg], timeout=3, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
      success = 0 == cp.returncode
      if not success:
         print(cp.stdout)
      self.assertTrue(success)
      
      
   def test_010_gt_valid_passes(self):
      """check genometools passes a valid GFF3 file"""
      cp = subprocess.run([gt_cmd, gff3_validator_tool, test_gff_file], timeout=6, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
      success = 0 == cp.returncode
      if not success:
         print(cp.stdout)
      self.assertTrue(success)

   def test_020_gt_invalid_fails(self):
      """check genometools fails an invalid GFF3 file"""
      cp = subprocess.run([gt_cmd, gff3_validator_tool, bad_gff_file], timeout=6, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
      self.assertFalse(0 == cp.returncode)
