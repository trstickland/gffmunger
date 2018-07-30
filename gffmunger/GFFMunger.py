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
      # CLI options
      if None == options:
         # passing None allow class to be constructed without options being defined; intended for example to permit tests
         # *not* recommended for normal usage, for which the 'gffmunger' script is provided
         self.verbose         = False
         self.input_file_arg  = '/dev/zero'
         self.output_file     = 'no_such_file'
         self.config_file     = 'config.yml'
      else:
         # this should be the normal case
         self.verbose         = options.verbose
         self.input_file_arg  = options.input_file
         self.output_file     = options.output_file
         self.config_file     = options.config
               
      # set up logger
      self.logger = logging.getLogger(__name__)
      if self.verbose:
         self.logger.setLevel(logging.DEBUG)
         self.logger.warning("logging seems to be broken at the moment, you'll only see warnings and errors :-/")
      else:
         self.logger.setLevel(logging.ERROR)
      
      # options from configuration file
      config_fh = open( os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', self.config_file),
                        'r'
                        );
      self.config = yaml.load(config_fh)
      def config_value_is_true(config_value):
         return( str(config_value).lower() == 'true' )
      try:
         self.keep_attr_value_order       = config_value_is_true(self.config['keep_attr_value_order'])
         self.attr_not_transferred        = self.config['attr_not_transferred']
         self.output_feature_sort         = self.config['output_feature_sort']
         self.gt_path                     = self.config['gt_path']
         self.gff3_validator_tool         = self.config['gff3_validator_tool']
         self.gff3_valiation_timeout      = self.config['gff3_validation_timeout']
         self.gffutils_db_filename        = str(self.config['gffutils_db_filename']).replace('<uid>',uuid.uuid4().hex)
         self.read_features_to_buffer     = config_value_is_true(self.config['read_features_to_buffer'])
      except KeyError as e:
         self.logger.error("required parameter "+str(e)+" missing from configuration in "+self.config_file)
         raise
      config_fh.close()
      
      self.logger.debug("Using genometools "+self.gt_path+" for validation with the tool "+self.gff3_validator_tool+" (timeout "+str(self.gff3_valiation_timeout)+")")

      if self.input_file_arg:
         self.logger.debug("Reading GFF3 input from "+ self.input_file_arg)
         if not os.path.exists(self.input_file_arg):
            self.logger.error("Input file does not exist: "+ self.input_file_arg)
            sys.exit(1)
      #reading from STDIN, but to avoid reading all the data into memory, we need a file to give to gffutils to parse
      #=> create unique temp filename to use an an input buffer
      else:
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
      try:
         self.get_gff3_source() # sets self.gff3_input_filename
         self.validate_GFF3(self.gff3_input_filename)
         self.import_gff3(self.gff3_input_filename)
         self.extract_GFF3_components(self.gff3_input_filename)
         self.move_annotations()
         self.export_gff3()
      except Exception:
         self.clean_up()
         raise
      self.clean_up()

               

   def clean_up(self):
      if self.temp_input_file and os.path.exists(self.temp_input_file):
         self.logger.debug("Deleting temporary input buffer "+ self.temp_input_file)
         os.remove(self.temp_input_file)
      if self.gffutils_db_filename and os.path.exists(self.gffutils_db_filename):
         self.logger.debug("Deleting gffutils database file "+ self.gffutils_db_filename)
         os.remove(self.gffutils_db_filename)
         db_backup_filename = self.gffutils_db_filename+".bak"
         if os.path.exists(db_backup_filename):
            os.remove(db_backup_filename)



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
      except Exception:
         blocksize = 4096
      while True:
         data = handle.read(blocksize)
         if not data:
               break
         yield data



   def get_gff3_source(self):
      """Returns name of GFF3 input file for gffutils.
      This may, trivially, return the input file name parameter.
      When input comes from STDIN, this is written to temporary file, and the name of that file is returned
      (but can safely be called more than once; won't attempt to re-read STDIN)"""
      
      # if called previously, self.gff3_input_filename will be defined, so it can be returned
      try:
         self.gff3_input_filename
      # otherwise, assume this is the first call...
      except Exception:
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
         self.gffutils_db = gffutils.create_db( gff_filename,
                                                dbfn                    = self.gffutils_db_filename,
                                                force                   = True,     # overwrite previous testing db file
                                                merge_strategy          = 'error',
                                                keep_order              = self.keep_attr_value_order,
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
         found_first_feature=False
         found_first_fasta=False
         for line in f:
            linenum += 1
            is_comment = line.startswith('#')
            if is_comment and not found_first_feature:
               # all comments prior to first feature line are assumed to be metadata
               # (or, at least, something to be preserved in the file header)
               # => add to metadata buffer
               self.input_metadata = append(line, self.input_metadata)
            elif is_comment and found_first_feature and not found_first_fasta:
               # this is a comment within the features, i.e. not metadata
               if self.read_features_to_buffer:
                  # if features are being stored as-read in a buffer, the comments should be added too
                  self.input_features = append(line, self.input_features)
            elif not is_comment and not found_first_fasta and not line.startswith('>'):
               # everything that isn't a comment, and occurs vefore the first FASTA line, is a feature
               # => this is a feature
               found_first_feature = True
               if self.read_features_to_buffer:
                  self.input_features = append(line, self.input_features)
            elif found_first_fasta or line.startswith('>'):
               # first line of FASTA, or any subsequent line
               found_first_fasta=True
               self.input_fasta = append(line, self.input_fasta)
            else:
               raise ValueError("Unexpected content in line {linenum} in GFF3:\n{line}")
            
      return(linenum)



   def move_annotations(self):
      """moves annotations from the polypeptide feature to mRNA feature"""
      num_polypeptide=0
      # this list caches all modified Feature objects
      # so they can be used to update the db after the gffutils_db.features_of_type iterator is finished
      modified_feature_cache=[]
      for this_polypeptide in self.gffutils_db.features_of_type('polypeptide'):
         num_polypeptide+=1
         
         # ignore polypeptide, with warning, if 'Derives_from' is missing
         # should be exception??
         if not 'Derives_from' in this_polypeptide.attributes:
            self.logger.warning("Ignoring polypeptide feature without a Derives_from attribute:\n"+str(this_polypeptide)+"\n")
            continue
         # get the Drives_from attribute (asserting presence of single Derives_from attribute)
         for n,derives_from in enumerate(this_polypeptide.attributes.get('Derives_from')):
            if n > 0:
               raise("panicked on encountering a polypeptide with more than one 'Derives_from' attribute")
         # ignore polypeptide in the relation isn't to an mRNA feature
         if not derives_from.endswith('mRNA'):
            self.logger.debug("Ignoring polypeptide feature that doesn't derive from mRNA feature")
            continue
         # get the parent feature (asserting presence of single parent)
         for n,this_parent in enumerate(self.gffutils_db.parents(derives_from)):
            if n > 0:
               raise("panicked on encountering a polypeptide with more than one parent feature")
         # get the related mRNA feature (asserting presence of single related mRNA)
         for n,this_mRNA in enumerate(self.gffutils_db.children(this_parent,featuretype='mRNA')):
            if n > 0:
               raise("panicked on encountering a polypeptide with more than one related mRNA attribute")
         # asssert ID attribute of the retrieved mRNA feature matches the Dervived_from attribute we started with
         # get the ID attribute (asserting presence of single Derives_from attribute)
         for n,mRNA_ID in enumerate(this_mRNA.attributes.get('ID')):
            if n > 0:
               raise("panicked on encountering an mRNA feature with more than one 'ID' attribute")
         if not mRNA_ID == derives_from:
            # self.logger.warning("Identifier mismatch:  mRNA with 'ID' attribute "+mRNA_ID+" was retreived for polypeptide with 'Derives_from' attribute "+derives_from)
            raise("Identifier mismatch:  mRNA with 'ID' attribute "+mRNA_ID+" was retreived for polypeptide with 'Derives_from' attribute "+derives_from)
         
         # create new set of attributes for the mRNA feature
         # these are a copy of those from the polypeptide (hence transferring annotations)...
         self.logger.debug("copying annotations to mRNA feature "+this_mRNA.attributes.get('ID')[0])
         new_mRNA_attributes = dict(this_polypeptide.attributes) # returns copy of this_polypeptide.attributes
         # ...except those attributes that shouldn't be transferred
         for not_copied in self.attr_not_transferred:
            if not_copied in new_mRNA_attributes:
               del new_mRNA_attributes[not_copied]
            if not_copied in this_mRNA.attributes:
               new_mRNA_attributes[not_copied] = this_mRNA.attributes.get(not_copied)
         # assign new attributes to mRNA feature
         this_mRNA.attributes = new_mRNA_attributes
         # cache the ammended Feature object
         modified_feature_cache.append(this_mRNA)
         
         # create new set of attributes for the polypeptide feature
         self.logger.debug("removing annotations from polypeptide feature "+this_polypeptide.attributes.get('ID')[0])
         new_polypeptide_attributes = {}
         # copy attribites that aren't tarnsferred to the mRNA feature
         for preserved_attribute in self.attr_not_transferred:
            if preserved_attribute in this_polypeptide.attributes:
               new_polypeptide_attributes[preserved_attribute] = this_polypeptide.attributes.get(preserved_attribute)
         # assign the enw attribites to the polypeptide feature
         this_polypeptide.attributes = new_polypeptide_attributes
         # cache the ammended mRNA Feature object
         modified_feature_cache.append(this_polypeptide)
                  
      self.logger.info("found "+str(num_polypeptide)+" polypeptide features")
      
      # remove the old versions of all the ammended features from the db
      self.gffutils_db.delete( modified_feature_cache )
      # insert the ammended features into the db
      # (this wee generator exists only to allow the logging durijng iteratingh through the Feature objects;
      # you could comment it out and just have gffutils.update() iterate over modified_feature_cache directly)
      def iterate_over_modified_features(cache):
         for this_feature in cache:
            self.logger.debug("modifying "+this_feature.attributes.get('ID')[0]+" in gffutils db")
            yield this_feature
      self.gffutils_db.update( iterate_over_modified_features(modified_feature_cache) )

      del modified_feature_cache # may be big



   def export_gff3(self):
      """Writes GFF3 to output file (if previously specified) or STDOUT
      Uses metadata and (if present) FASTA from the GFF3 input; these should be unalatered
      Features are written from the gffutils database, so will refect whatever munging
      was done via the gffutils API
      """
      if self.gffutils_db is None:
         raise("Must import some GFF3 data before exporting")
      
      if self.output_file is not None:
         self.logger.debug("Exporting GFF3 to file "+ self.output_file)
         handle = open(self.output_file, "wt")
      else:
         self.logger.debug("Exporting GFF3 to STDOUT")
         handle = sys.stdout

      handle.write( self.input_metadata )
      num_features_written=0
      for this_feature in self.gffutils_db.all_features(order_by=self.output_feature_sort):
         num_features_written+=1
         handle.write( str(this_feature)+"\n" )
      self.logger.info("extracted and wrote "+str(num_features_written)+" features from gffutils db")
      if self.input_fasta is not None:
         handle.write( self.input_fasta )
         
      handle.close()
      
      return(True)
      
