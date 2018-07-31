import argparse
import os
import unittest
import uuid

from gffmunger.GFFMunger import GFFMunger

test_modules_dir  = os.path.dirname(   os.path.realpath( __file__ ) )
data_dir          = os.path.join(      test_modules_dir, 'data' )
test_gff_file     = os.path.join(      data_dir,         'SMALL_SAMPLE_INCL_FASTA.gff3.gz' )

expected_num_input_lines   = 2082

class IO_Tests(unittest.TestCase):
   
   @classmethod
   def setUpClass(self):
      parser = argparse.ArgumentParser()
      parser.add_argument('--verbose',       '-v',    action='store_true',    help = 'Turn on debugging [%(default)s]',                default = False)
      #parser.add_argument('--version',               action='version',       version = str(version))
      parser.add_argument('--no-validate',   '-n',    action='store_true',    help = 'Do not validate the input GFF3 [%(default)s]',   default = False)
      parser.add_argument('--force',         '-f',    action='store_true',    help = 'Force writing of output file, even if it already exists [%(default)s]', default = False)
      parser.add_argument('--input-file',    '-i',                            help = 'Input file [STDIN]')
      parser.add_argument('--output-file',   '-o',                            help = 'Output file [STDOUT]')
      parser.add_argument('--config',        '-c',                            help = 'Config file [%(default)s]',                      default = 'config.yml')
      self.output_file = __file__+'.'+uuid.uuid4().hex+'.gff3'
      self.gffmunger = GFFMunger(parser.parse_args([  '--input-file',   test_gff_file,
                                                      '--output-file',  self.output_file,
                                                      ]
                                                   )
                                 )

   @classmethod
   def tearDownClass(self):
      if self.output_file and os.path.exists(self.output_file):
         os.remove(self.output_file)      

     
   def test_000_gff3_input(self):
      """check functions for retreiving input filename, and reading and validating the file"""
      self.assertEqual(test_gff_file,  self.gffmunger.get_gff3_source())
      self.assertTrue(self.gffmunger.validate_GFF3(test_gff_file))
      self.assertEqual(self.gffmunger.gffutils_db_filename, self.gffmunger.import_gff3())
      self.assertEqual(self.gffmunger.gffutils_db_filename, self.gffmunger.import_gff3(test_gff_file))
      self.assertEqual(expected_num_input_lines,            self.gffmunger.extract_GFF3_components())
      self.assertEqual(expected_num_input_lines,            self.gffmunger.extract_GFF3_components(test_gff_file))

   def test_010_gff3_output(self):
      """checks functions for exporting GFF3 to a file"""
      self.assertIsNotNone(self.gffmunger.output_file)
      self.assertTrue(self.gffmunger.export_gff3())
      os.remove(self.gffmunger.output_file)
