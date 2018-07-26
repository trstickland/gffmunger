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
   
   def extract_GFF3_components(self, gff_filename):
      """Extracts separate components of a GFF3 file:  metadata, features and FASTA
      Stores these as raw text buffers, as read from the file, unescaped
      - Metadata are lines at the beginning of the GFF3 starting '##'
      - Features start from the first line that does not begin '##', and include following lines that start '#'
        (which are human-readable comments rather than metadata).  
      - FASTA data are read from the first line beginning '>', to the end of the file.  Note that  '>' is
        *specifically* prohibited as the first character of a 
      An artefact of this is that if there is a line such as '## FASTA' that indicate the start of FASTA data,
      this will occur at the end of the feature buffer and not at the start of the FASTA buffer.  This is
      safe as it will be ignored as a comment when parsing features (the GFF3 spec. requires parsers ignore
      comments), whereas it could /(?? check)/ be a problem in the FASTA.  Hope this doesn't confuse anyone."""
      self.input_metadata  = None
      self.input_features  = None
      self.input_fasta     = None
      linenum              = 0
      with self.open_text_file(gff_filename) as f:
         # appends 'new' to 'buf', without error if buf is None
         def append(new, buf=''):
            if new is None:
               return(buf)
            return(buf + new)
         for line in f:
            linenum += 1
            is_metadata = line.startswith('##')
            # line looks like metadata, and haven't read any feature lines yet
            # => add to metadata
            if is_metadata and self.input_features is None:
               #print("METADATA: "+line)
               self.input_metadata = append(self.input_metadata, line)
            # line looks like metadata, but we *have* read feature feature line(s)
            # => treat as a comment, and add to features buffer
            elif is_metadata and self.input_features is not None:
               #print("FEATURE: "+line)
               self.input_features = append(self.input_features, line)
            # not metadata, and haven't read and FASTA yet
            # => add to features
            elif not is_metadata and not line.startswith('>') and self.input_fasta is None:
               #print("FEATURE: "+line)
               self.input_features = append(self.input_features, line)
            # first line of FASTA, or any subsequent line
            # => add to fasta
            elif line.startswith('>') or self.input_fasta is not None:
               #print("FASTA: "+line)
               self.input_fasta = append(self.input_fasta, line)
            else:
               raise ValueError("Unexpected content in line {linenum} in GFF3:\n{line}")
      return(linenum)
      
      
      
