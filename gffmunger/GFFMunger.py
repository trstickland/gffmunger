import argparse
import gffutils
import gzip
import logging
import os
import subprocess
import sys
import time
import uuid
import warnings
import yaml

from Bio import SeqIO

class GFFMunger:
   
   def __init__(self,options):
      # allow class to be constructed without options being defined
      # intended for example to permit tests, not recommended for normal usage, for which the 'gffmunger' script is provided
      if None == options:
         self.verbose         = False
         self.input_file_arg  = '/dev/zero'
         self.output_file     = 'no_such_file'
         self.config_file     = 'config.yml'
      else:
         self.verbose         = options.verbose
         self.input_file_arg  = options.input_file
         self.output_file     = options.output_file
         self.config_file     = options.config
               
      self.logger = logging.getLogger(__name__)
      if self.verbose:
         self.logger.setLevel(logging.DEBUG)
         self.logger.warning("logging seems to be broken at the moment, you'll only see warnings and errors :-/")
      else:
         self.logger.setLevel(logging.ERROR)
      
      
      config_fh = open( os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', self.config_file),
                        'r'
                        );
      self.config = yaml.load(config_fh)
      self.gt_path                  = self.config['gt_path']
      self.gff3_validator_tool      = self.config['gff3_validator_tool']
      self.gff3_valiation_timeout   = self.config['gff3_validation_timeout']
      self.gffutils_db_filename     = str(self.config['gffutils_db_filename']).replace('<uid>',uuid.uuid4().hex)
      self.read_features_to_buffer  = str(self.config['read_features_to_buffer']).lower() == 'true'
            
      config_fh.close()
      self.logger.debug("Using genometools "+self.gt_path+" for validation with the tool "+self.gff3_validator_tool+" (timeout "+str(self.gff3_valiation_timeout)+")")

      if self.input_file_arg:
         self.logger.debug("Reading GFF3 input from "+ self.input_file_arg)
         if not os.path.exists(self.input_file_arg):
            self.logger.error("Input file does not exist: "+ self.input_file_arg)
            sys.exit(1)
      else: #reading from STDIN, but gffutils requires file to parse => create unique temp filename
         self.logger.debug("Reading GFF3 input from STDIN")
         self.temp_input_file = str(self.config['temp_input_file']).replace('<uid>',uuid.uuid4().hex)
         self.logger.debug("Temporary input buffer will be "+ self.temp_input_file)
         if os.path.exists(self.temp_input_file):
            self.logger.error("Something badly wrong :-/   Should have a unique filename for the temporary input buffer, but it already exists: "+ self.temp_input_file)
            sys.exit(1)
            
      if self.output_file:
         self.logger.debug("Writing output to "+ self.output_file)
         if os.path.exists(self.output_file):
            self.logger.error("The output file already exists, please choose another filename: "+ self.output_file)
            sys.exit(1)
      else:
         self.logger.debug("Writing output to STDOUT")
               


            
   def run(self):
      self.get_gff3_source() # sets self.gff3_input_filename
      self.validate_GFF3(self.gff3_input_filename)
      self.import_gff3(self.gff3_input_filename)
      self.extract_GFF3_components(self.gff3_input_filename)
      
      self.export_gff3()
      
      # clean up
      if self.temp_input_file and os.path.exists(self.temp_input_file):
         self.logger.debug("Deleting temporary input buffer "+ self.temp_input_file)
         os.remove(self.temp_input_file)
      if self.gffutils_db_filename and os.path.exists(self.gffutils_db_filename):
         self.logger.debug("Deleting gffutils database file "+ self.gffutils_db_filename)
         os.remove(self.gffutils_db_filename)



   def open_text_file(self, filename):
      """Opens a possibly gzipped file for reading as text, returns handle"""
      if filename.endswith('.gz'):
         handle = gzip.open(filename, "rt")
      else:
         handle = open(filename, "r")
      return(handle)



   def read_in_blocks(self, handle):
      """Generator to read from a file handle by block. Tries to use stat to figure out system blocksize, defaults to 4096"""
      try:
         blocksize = os.stat(__file__).st_blksize
      except:
         blocksize = 4096
      while True:
         data = handle.read(blocksize)
         if not data:
               break
         yield data



   def get_gff3_source(self):
      """Returns name of GFF3 input file for gffutils.
      This may, trivially, return the input file name parameter.
      When input comes from STDIN, this is written to a unique temporary file, and the name of that file is returned
      (but can safely be called more than once; won't attempt to re-read STDIN)"""
      
      # if called previously, self.gff3_input_filename will be defined, so it can be returned
      try:
         self.gff3_input_filename
      # otherwise, assume this is the first call...
      except:
         if self.input_file_arg:
            # given file name as parameter => use that
            self.gff3_input_filename = self.input_file_arg
         elif self.temp_input_file:
            # reading from STDIN, with temporary file as input buffer
            with open(self.temp_input_file, 'w') as f:
               for block in self.read_in_blocks(sys.stdin):
                  f.write( block )
            self.gff3_input_filename = self.temp_input_file
         else:
            # shouldn't be reachable if variables initialized correctly
            raise ValueError("Don't have an input file to read, and can't read from STDIN without name of a temporary file to use as an input buffer")
      
      # return file name
      self.logger.debug("Using GFF3 data from %s", self.gff3_input_filename)
      return(self.gff3_input_filename)



   def validate_GFF3(self, gff_filename, silent=False): 
      """Validates GFF3 file.
      If valid, True is returned; if invalid, validator STDERR output is printed and False is returned
      Validator errors are printed to STDOUT; these can be supressed by passing the optional flag 'silent'"""
      self.logger.debug("Validating GFF3 file "+ gff_filename)
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
      self.logger.debug("Validating FASTA file "+ fasta_filename)
      with self.open_text_file(fasta_filename) as handle:
         fasta = SeqIO.parse(handle, "fasta")
         is_fasta = any(fasta)   # False when `fasta` is empty, i.e. wasn't a FASTA file
         if not is_fasta and not silent:
            print(fasta_filename+" is not a valid FASTA file")
      return(is_fasta)



   def import_gff3(self, gff_filename=None):
      """Optionally path of GFF3 file; otherwise this is retrieved using get_gff3_source()
      Imports GFF3 from the file into gffutils"""
      if not gff_filename:
         gff_filename = self.get_gff3_source()
      self.logger.debug("Importing using gffutils, from GFF3 file "+ gff_filename)
      with warnings.catch_warnings():
         if not self.verbose:
            warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper",  ResourceWarning,           "gffutils", 133 )
            warnings.filterwarnings("ignore", "generator '_FileIterator\.",         PendingDeprecationWarning, "gffutils", 186 )
            warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper",  ResourceWarning,           "gffutils", 668 )
         self.test_gff_db = gffutils.create_db( gff_filename,
                                                dbfn                    = self.gffutils_db_filename,
                                                force                   = True,     # overwrite previous testing db file
                                                merge_strategy          = 'error',
                                                keep_order              = False,    # True doesn't appear to maintain attribute order :-/  (and turning this off may be faster)
                                                sort_attribute_values   = False
                                                )
      return(self.gffutils_db_filename)



   def extract_GFF3_components(self, gff_filename=None):
      """Optionally path of GFF3 file; otherwise this is retrieved using get_gff3_source()
      Extracts separate components from the GFF3 file:  metadata, features and FASTA
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
      if not gff_filename:
         gff_filename = self.get_gff3_source()
      self.logger.debug("extracting metadta, features and (possibly) FASTA from GFF3 file "+ gff_filename)
      
      self.input_metadata  = None
      self.input_features  = None
      self.input_fasta     = None
      linenum              = 0
      with self.open_text_file(gff_filename) as f:
         # wee function that appends 'new' to 'buf', without error if buf is None
         def append(new, buf=''):
            if new is None:
               return(buf)
            if buf is None:
               return(new)
            return(buf + new)
         for line in f:
            linenum += 1
            is_metadata = line.startswith('##')
            if is_metadata and self.input_features is None:
               # line looks like metadata, and haven't read any feature lines yet
               # => add to metadata
               self.input_metadata = append(line, self.input_metadata)
            elif is_metadata and self.input_features is not None:
               # line looks like metadata, but we *have* read feature feature line(s)
               # => treat as a comment, and add to features buffer
               if self.read_features_to_buffer:
                  self.input_features = append(line, self.input_features)
            elif not is_metadata and not line.startswith('>') and self.input_fasta is None:
               # not metadata, and haven't read and FASTA yet
               # => add to features
               if self.read_features_to_buffer:
                  self.input_features = append(line, self.input_features)
            # first line of FASTA, or any subsequent line
            # => add to fasta
            elif line.startswith('>') or self.input_fasta is not None:
               self.input_fasta = append(line, self.input_fasta)
            else:
               raise ValueError("Unexpected content in line {linenum} in GFF3:\n{line}")
            
      return(linenum)



   def export_gff3(self):
      """Writes GFF3 to output file (if previously specified) or STDOUT
      Uses metadata and (if present) FASTA from the GFF3 input; these should be unalatered
      Features are written from the gffutils database, so will refect whatever munging
      was done via the gffutils API
      """
      if self.output_file is not None:
         self.logger.debug("Exporting GFF3 to file "+ self.output_file)
         handle = open(self.output_file, "wt")
      else:
         self.logger.debug("Exporting GFF3 to STDOUT")
         handle = sys.stdout

      handle.write( self.input_metadata )

      ############# TEMPORARY ##############
      #if self.read_features_to_buffer:
      handle.write( self.input_features )
      ######################################   
      
      if self.input_fasta is not None:
         handle.write( self.input_fasta )
         
      handle.close()
      
      return(True)
      
