import argparse
import gzip
import logging
import os
import subprocess
import sys
import time
import yaml

from Bio import SeqIO

class GFFMunger:
   def __init__(self,options):
      
      # allow class to be constructed without options being defined
      # intended for example to permit tests, not recommended for normal usage, for which the 'gffmunger' script is provided
      if None == options:
         self.verbose      = False
         self.output_file  = 'gffmunger_output.gff3'
         self.config_file  = 'config.yml'
      else:
         self.verbose      = options.verbose
         self.output_file  = options.output_file
         self.config_file  = options.config
      
      self.logger = logging.getLogger(__name__)
      
      config_fh = open( os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', self.config_file),
                        'r'
                        );
      self.config = yaml.load(config_fh)
      self.gt_path                  = self.config['gt_path']
      self.gff3_validator_tool      = self.config['gff3_validator_tool']
      self.gff3_valiation_timeout   = self.config['gff3_validation_timeout']
      config_fh.close()

      if self.output_file and os.path.exists(self.output_file):
         self.logger.error("The output file already exists, please choose another filename: "+ self.output_file)
         sys.exit(1)
               
      if self.verbose:
         self.logger.setLevel(logging.DEBUG)
      else:
         self.logger.setLevel(logging.ERROR)
         
   def run(self):
      print(__name__+".run() was called")
      
   def open_text_file(self, filename):
      """Opens a possibly gzipped file for reading as text, returns handle"""
      if filename.endswith('.gz'):
         handle = gzip.open(filename, "rt")
      else:
         handle = open(filename, "r")
      return(handle)

   def validate_GFF3(self, gff_filename, silent=False): 
      """Validates GFF3 file.
      Pass path of GFF3 file; if valid, True is returned; if invalid, validator STDERR output is printed and False is returned
      Validator errors are printed to STDOUT; these can be supressed by passing the optional flag 'silent'"""
      cp = subprocess.run( [self.gt_path, self.gff3_validator_tool, gff_filename],
                           timeout=self.gff3_valiation_timeout,
                           stderr=subprocess.PIPE, stdout=subprocess.PIPE)
      if 0 == cp.returncode:
         return True
      if not silent:
         print(cp.stderr)
      return False
   
   def validate_FASTA(self, fasta_filename, silent=False): 
      """Validates FASTA file.
      Pass path of FASTA file; if valid, True is returned; if invalid, validator STDERR output is printed and False is returned
      Validation failure message printed to STDOUT; this can be supressed by passing the optional flag 'silent'"""
      with self.open_text_file(fasta_filename) as handle:
         fasta = SeqIO.parse(handle, "fasta")
         is_fasta = any(fasta)   # False when `fasta` is empty, i.e. wasn't a FASTA file
         if not is_fasta and not silent:
            print(fasta_filename+" is not a valid FASTA file")
      return(is_fasta)
