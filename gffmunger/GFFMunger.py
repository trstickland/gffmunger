import gffutils
import gzip
import logging
import os
import re
import subprocess
import sys
import time
import uuid
import warnings
import yaml

from Bio import SeqIO
from pyfaidx import Fasta

class GFFMunger:
   
   def __init__(self,options):
      # CLI options
      if None == options:
         # passing None allow class to be constructed without options being defined; intended for example to permit tests
         # *not* recommended for normal usage, for which the 'gffmunger' script is provided
         self.verbose         = False
         self.quiet           = False
         self.novalidate      = False
         self.force           = False
         self.fasta_file_arg  = None
         self.input_file_arg  = '/dev/zero'
         self.output_file     = 'no_such_file'
         self.config_file     = 'config.yml'
      else:
         # this should be the normal case
         self.verbose         = options.verbose
         self.quiet           = options.quiet
         self.novalidate      = options.no_validate
         self.force           = options.force
         self.fasta_file_arg  = options.fasta_file
         self.input_file_arg  = options.input_file
         self.output_file     = options.output_file
         self.config_file     = options.config
               
      # set up logger
      self.logger = logging.getLogger(__name__)
      if self.verbose:
         self.logger.setLevel(logging.INFO)
         self.logger.warning("logging.debug & logging.info seem to be broken at the moment, you'll only see warnings & above :-/")
      elif self.quiet:
         self.logger.setLevel(logging.CRITICAL)
      else:
         self.logger.setLevel(logging.WARNING)
      
      # options from configuration file
      config_filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', self.config_file)
      try:
         config_fh   = open(config_filename, 'r');
         self.config = yaml.load(config_fh)
      except Exception:
         self.logger.critical("Can't read configuration file"+config_filename)
         raise
      def config_value_is_true(config_value):
         return( str(config_value).lower() == 'true' )
      try:
         self.keep_attr_value_order       = config_value_is_true(self.config['keep_attr_value_order'])
         self.attr_not_transferred        = self.config['attr_not_transferred']
         self.output_feature_sort         = self.config['output_feature_sort']
         self.annotated_feature_types     = self.config['annotated_feature_types']
         self.only_transfer_anot_to_mRNA  = config_value_is_true(self.config['only_transfer_anot_to_mRNA'])
         self.gt_path                     = self.config['gt_path']
         self.gff3_validator_tool         = self.config['gff3_validator_tool']
         self.gff3_valiation_timeout      = self.config['gff3_validation_timeout']
         self.gffutils_db_filename        = str(self.config['gffutils_db_filename']).replace('<uid>',uuid.uuid4().hex)
         self.read_features_to_buffer     = config_value_is_true(self.config['read_features_to_buffer'])
      except KeyError as e:
         self.logger.critical("required parameter "+str(e)+" missing from configuration in "+self.config_file)
         raise
      config_fh.close()

      self.logger.debug("Using genometools "+self.gt_path+" for validation with the tool "+self.gff3_validator_tool+" (timeout "+str(self.gff3_valiation_timeout)+")")

      if self.fasta_file_arg:
         self.logger.info("Reading FASTA from "+ self.fasta_file_arg)
         if not os.path.exists(self.fasta_file_arg):
            self.logger.critical("FASTA file does not exist: "+ self.fasta_file_arg)
            sys.exit(1)

      if self.input_file_arg and "-" != str(self.input_file_arg):
         self.logger.info("Reading GFF3 input from "+ self.input_file_arg)
         if not os.path.exists(self.input_file_arg):
            self.logger.critical("Input file does not exist: "+ self.input_file_arg)
            sys.exit(1)
      #reading from STDIN, but to avoid reading all the data into memory, we need a file to give to gffutils to parse
      #=> create unique temp filename to use an an input buffer
      else:
         if self.input_file_arg:
            self.input_file_arg = None
         self.logger.info("Reading GFF3 input from STDIN")
         self.temp_input_file = str(self.config['temp_input_file']).replace('<uid>',uuid.uuid4().hex)
         self.logger.debug("Temporary input buffer will be "+ self.temp_input_file)
         if os.path.exists(self.temp_input_file):
            self.logger.critical("Something badly wrong :-/   Should have a unique filename for the temporary input buffer, but it already exists: "+ self.temp_input_file)
            sys.exit(1)
            
      if self.output_file and "-" != str(self.output_file):
         self.logger.info("Writing output to "+ self.output_file)
         if not self.force and os.path.exists(self.output_file):
            self.logger.critical("The output file already exists, please choose another filename: "+ self.output_file)
            sys.exit(1)
      else:
         if self.output_file:
            self.output_file = None
         self.logger.debug("Writing output to STDOUT")
      

            
   def run(self):
      try:
         # get GFF3 input, stdin or file; sets self.gff3_input_filename
         self.get_gff3_source() 
         # validate GFF3 if required
         if not self.novalidate:
            self.validate_GFF3(self.gff3_input_filename)
         # import GFF3
         self.import_gff3(self.gff3_input_filename)
         # if FASTA is being read from separate file...
         if self.fasta_file_arg:
            # ...validate if required...
            if not self.novalidate:
               self.validate_FASTA(self.fasta_file_arg)
            # ...and import
            self.import_fasta(self.fasta_file_arg)
         # read GFF3 metadta (and poss. other bits) into text buffer(s)
         self.extract_GFF3_components(self.gff3_input_filename)
         # transfer annotations from polypeptide features to the feature they derived from
         self.move_annotations()
         # write new GFF3 to file or stdout
         self.export_gff3()
         # if GFF3 file was written, validate it if required
         if self.output_file is not None and not self.novalidate:
            self.validate_GFF3(self.output_file)
      # whack temporary files
      except Exception:
         self.clean_up()
         raise
      self.clean_up()



   def clean_up(self):
      if hasattr(self, 'temp_input_file') and self.temp_input_file and os.path.exists(self.temp_input_file):
         self.logger.debug("Deleting temporary input buffer "+ self.temp_input_file)
         os.remove(self.temp_input_file)
      if hasattr(self, 'gffutils_db_filename') and self.gffutils_db_filename and os.path.exists(self.gffutils_db_filename):
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
      self.logger.info("Validating GFF3 file "+ gff_filename)
      cp = subprocess.run( [self.gt_path, self.gff3_validator_tool, gff_filename],
                           timeout=self.gff3_valiation_timeout,
                           stderr=subprocess.PIPE, stdout=subprocess.PIPE)
      if 0 == cp.returncode:
         return True
      if not silent:
         print(gff_filename+" is not valid GFF3:")
         print(cp.stderr)
      return False



   def validate_FASTA(self, fasta_filename, silent=False): 
      """Validates FASTA file.
      Pass path of FASTA file; if valid, True is returned; if invalid, validator STDERR output is printed and False is returned
      Validation failure message printed to STDOUT; this can be supressed by passing the optional flag 'silent'"""
      self.logger.info("Validating FASTA file "+ fasta_filename)
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



   def import_fasta(self, fasta_filename=None):
      """Optionally path of FASTA file; otherwise use self.fasta_file_arg
      Imports FASTA from the file using pyfaidx.Fasta"""
      if not fasta_filename:
         fasta_filename = self.fasta_file_arg
      self.logger.debug("Importing FASTA using pyfaidx.Fasta, from "+ fasta_filename)
      self.faidx = Fasta(fasta_filename)
      return(self.faidx)
      
      
      
   def extract_GFF3_components(self, gff_filename=None):
      """Optionally pass path of GFF3 file; otherwise this is retrieved using get_gff3_source()
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
         
      if self.fasta_file_arg:
         if self.read_features_to_buffer:
            self.logger.debug("extracting metadata and features (but not FASTA) from GFF3 file "+ gff_filename)
         else:
            self.logger.debug("extracting metadata (only) from GFF3 file "+ gff_filename)
      elif self.read_features_to_buffer:
            self.logger.debug("extracting metadata, features and FASTA (if there is any) from GFF3 file "+ gff_filename)
      else:
            self.logger.debug("extracting metadata and FASTA (if there is any) from GFF3 file "+ gff_filename)
      
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
               # everything that isn't a comment, and occurs before the first FASTA line, is a feature
               # => this is a feature
               # if features aren't being read into a buffer, and the FASTA is being read from a separate file,
               # we can bale out right here
               if not self.read_features_to_buffer and self.fasta_file_arg:
                  self.logger.debug("finished reading GFF3 file at the end of the metadata")
                  return(linenum-1)
               found_first_feature = True
               if self.read_features_to_buffer:
                  self.input_features = append(line, self.input_features)
            elif found_first_fasta or line.startswith('>'):
               # first line of FASTA, or any subsequent line
               # if the FASTA is being read from a separate file,
               # we can bale out right here
               if self.fasta_file_arg:
                  self.logger.debug("finished reading GFF3 file at the end of the features")
                  return(linenum-1)
               found_first_fasta=True
               self.input_fasta = append(line, self.input_fasta)
            else:
               raise ValueError("Unexpected content in line {linenum} in GFF3:\n{line}")
            
      self.logger.debug("reached end of GFF3 file")
      return(linenum)



   def move_annotations(self):
      """moves annotations from the polypeptide feature to the feature from which it derives (e.g. mRNA)"""
      num_polypeptide=0
      # this list caches all modified Feature objects
      # so they can be used to update the db after the gffutils_db.features_of_type iterator is finished
      modified_feature_cache=[]
      for this_polypeptide in self.gffutils_db.features_of_type('polypeptide'):
         num_polypeptide+=1
         
         this_derives_from_feature = self.get_derives_from_feature(this_polypeptide)
         # return value of None indicates polypeptide shoukld be ignored, but it's safe to continue
         if this_derives_from_feature is None:
            continue
         
         # log warning if the returned feature type is not one expected to be annotated
         if not this_derives_from_feature.featuretype in self.annotated_feature_types:
            self.logger.warning("Polypeptide %s described as derivant of %s which is an unexpected type: %s",
                                str(this_polypeptide.attributes.get('ID')[0]),
                                str(this_derives_from_feature.attributes.get('ID')[0]),
                                this_derives_from_feature.featuretype,
                                )
         
         # create new set of attributes for the Derives_from feature
         # these are a copy of those from the polypeptide (hence transferring annotations)...
         self.logger.debug("copying annotations to feature from which polypeptide derives "+this_derives_from_feature.attributes.get('ID')[0])
         new_derives_from_feature_attributes = dict(this_polypeptide.attributes) # returns copy of this_polypeptide.attributes
         # ...except those attributes that shouldn't be transferred
         for not_copied in self.attr_not_transferred:
            if not_copied in new_derives_from_feature_attributes:
               del new_derives_from_feature_attributes[not_copied]
            if not_copied in this_derives_from_feature.attributes:
               new_derives_from_feature_attributes[not_copied] = this_derives_from_feature.attributes.get(not_copied)
         # assign new attributes to Derives_from feature
         this_derives_from_feature.attributes = new_derives_from_feature_attributes
         # cache the ammended Feature object
         modified_feature_cache.append(this_derives_from_feature)
         
         # create new set of attributes for the polypeptide feature
         self.logger.debug("removing annotations from polypeptide feature "+this_polypeptide.attributes.get('ID')[0])
         new_polypeptide_attributes = {}
         # make a copy of all attributes that aren't to be transferred to the Derives_from feature
         for preserved_attribute in self.attr_not_transferred:
            if preserved_attribute in this_polypeptide.attributes:
               new_polypeptide_attributes[preserved_attribute] = this_polypeptide.attributes.get(preserved_attribute)
         # assign the new attribites to the polypeptide feature
         this_polypeptide.attributes = new_polypeptide_attributes
         # cache the ammended polypeptide Feature object
         modified_feature_cache.append(this_polypeptide)
                  
      self.logger.info("found "+str(num_polypeptide)+" polypeptide features")
      
      self.check_for_anotations(modified_feature_cache)
      
      # remove the old versions of all the ammended features from the db
      self.gffutils_db.delete( modified_feature_cache )
      # insert the ammended features into the db
      ## (this wee generator exists only to allow the logging whilst iterating through the Feature objects)
      #def iterate_over_modified_features(cache):
         #for this_feature in cache:
            #self.logger.debug("modifying "+this_feature.attributes.get('ID')[0]+" in gffutils db")
            #yield this_feature
      #self.gffutils_db.update( iterate_over_modified_features(modified_feature_cache) )
      self.gffutils_db.update( modified_feature_cache )

      del modified_feature_cache # may be big



   # N.B. this must only return None to indicate when polypeptide should be ignored, but it's safe to continue;
   # raise an exception when there's an error that can't be ignored
   def get_derives_from_feature(self, polypeptide_feature):
      """Pass a gffutils.Feature object representing a polypeptide
      Returns the gffutils.Feature object representing the feature from which the polypeptide
      derives, as specified by the Derives_from attribute.
      Returns None, with a logger warning, if the feature can't be identified for some reason
      Only raises exception on encountering a something so unexpected that we can't safely continue."""
      # ignore polypeptide, with warning, if 'Derives_from' is missing
      if not 'Derives_from' in polypeptide_feature.attributes:
         self.logger.error("Ignoring polypeptide feature without a Derives_from attribute:\n"+str(polypeptide_feature)+"\n")
         return(None)
      # get the polypeptide ID
      num_polypeptide_ID = 0
      for polypeptide_ID in polypeptide_feature.attributes.get('ID'):
         num_polypeptide_ID += 1
      if not 1 == num_polypeptide_ID:
         raise AssertionError("polypeptide "+polypeptide_ID+" must have exactly one 'ID' attribute, found "+str(num_polypeptide_ID)+" in feature line"+str(polypeptide_feature))
      # get the Derives_from attribute (asserting presence of single Derives_from attribute)
      num_derives_from = 0
      for derives_from in polypeptide_feature.attributes.get('Derives_from'):
         num_derives_from += 1
      if not 1 == num_derives_from:
         raise AssertionError("polypeptide must have exactly one 'Derives_from' attribute, found "+str(num_derives_from)+" in "+polypeptide_ID)
         
      # ignore polypeptide in the relation isn't to an mRNA feature
      if self.only_transfer_anot_to_mRNA and not derives_from.endswith('mRNA'):
         self.logger.debug("Ignoring polypeptide feature that doesn't derive from mRNA feature")
         return(None)
      
      # get the parent feature (asserting presence of single parent)
      num_parents = 0
      for this_parent in self.gffutils_db.parents(derives_from):
         num_parents += 1
      if not 1 == num_parents:
         # raise AssertionError("a polypeptide must have exactly one parent feature, found "+str(num_parents)+" parents of "+derives_from)
         self.logger.error("a polypeptide must have exactly one parent feature, found "+str(num_parents)+" parents of "+derives_from+": cannot transfer its annotations")
         return(None)

      # search amongst sibling features to find the one the polypeptide derives from
      num_matches = 0
      derives_from_feature = None
      for this_child in self.gffutils_db.children(this_parent):
         this_child_ID = None
         # get this child's ID
         num_id = 0
         for this_child_ID in this_child.attributes.get('ID'):
            num_id += 1
         # assert one only ID
         if not 1 == num_id:
            raise AssertionError("a feature must have exactly one 'ID' attribute, found "+str(num_id)+" in feature line "+str(this_child))
         # check for match
         if this_child_ID == derives_from:
            num_matches += 1
            derives_from_feature = this_child
      # assert one match
      if not 1 == num_matches:
         # raise AssertionError("a polypeptide Drives_from attribute must match the ID of exactly one sibling feature, found "+str(num_matches)+ " matches for "+derives_from)
         self.logger.error("polypeptide "+num_polypeptide_ID+" apparently derives from "+str(num_matches)+" siblings (should be exactly one): cannot transfer its annotations")
         return(None)
      
      # catch-all assertion ensuring we're returning a Feature
      if derives_from_feature is None:
         raise AssertionError("failed to find sibling feature matching polypeptide Drives_from "+derives_from+" of polypeptide "+polypeptide_ID)
      return(derives_from_feature)



   def check_for_anotations(self, annotated_features):
      """Checks that annotations have been transferred to features that should be annotated
      Pass the list of newly annotated features.
      Logs a warning for each feature without annotation, that should be annotated"""
      annotated_feature_ids = []
      for this_feature in annotated_features:
         annotated_feature_ids.append( this_feature.attributes.get('ID')[0] )
      for annotated_type in self.annotated_feature_types:
         for this_feature in self.gffutils_db.features_of_type(annotated_type):
            if not this_feature.attributes.get('ID')[0] in annotated_feature_ids:
               self.logger.warning("Feature "+str(this_feature.attributes.get('ID')[0])+" ("+this_feature.featuretype+") has no annotation because no derivant polypeptide was found")
   
   
   
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

      # write metadata
      handle.write( self.input_metadata )
      
      # write features
      num_features_written=0
      for this_feature in self.gffutils_db.all_features(order_by=self.output_feature_sort):
         num_features_written+=1
         handle.write( str(this_feature)+"\n" )
      self.logger.info("extracted and wrote "+str(num_features_written)+" features from gffutils db")
      
      # write fasta
      handle.write("##FASTA\n")
      if self.fasta_file_arg is not None:
         # using a separate FASTA file; write sequences from that file
         # if a sequence exists in the FASTA file but is not referenced in the input GFF3, it is *not* written
         # if a sequence is referenced in the input GFF3 but doesn't exist in the FASTA, KeyError is raised
         if self.faidx is None:
            raise AssertionError("When reading dequences from a separate FASTA file, a pxfaidx.Fasta object must be created before export")
         num_seq_written   = 0
         for this_seq_id in self.gffutils_db_sequences():
            try:
               handle.write(  ">"+this_seq_id               +"\n"+
                              str(self.faidx[this_seq_id])  +"\n"
                              )
               num_seq_written+=1
               self.logger.debug("Writing FASTA sequence "+str(num_seq_written)+": "+str(this_seq_id))
            except KeyError:
               self.logger.error("The GFF3 input included sequence "+this_seq_id+" which was not found in the FASTA input:  output will not include this in the FASTA")
      else:
         # using FASTA from the input GFF3 => write whatever was in the input (possibly nowt)
         if self.input_fasta is not None:
            handle.write( self.input_fasta )
         
      handle.close()
      
      return(True)
      


   def gffutils_db_sequences(self):
      """generator that yields each sequence in the gffutils db"""
      sequences_found = []
      for this_gene in self.gffutils_db.features_of_type('gene', order_by='seqid,start'):
         if this_gene.seqid in sequences_found:
            continue
         sequences_found.append(this_gene.seqid)
         yield this_gene.seqid
