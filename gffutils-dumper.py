import os
import gffutils
import warnings
import dumper

# provisioning on ubuntu/xenial64 virtual box:
# apt-get install -y git python3 python-setuptools python3-biopython python3-pip
# pip3 install dumper gffutils


data_dir          = os.path.dirname(   os.path.realpath( __file__ ) )
test_gff_file     = os.path.join(      data_dir,         'SAMPLE.gff3.gz' )
test_gff_db_file  = os.path.join(      data_dir,         'gffutils_test.db' )

with warnings.catch_warnings():
   warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper",  ResourceWarning,           "gffutils", 133 )
   warnings.filterwarnings("ignore", "generator '_FileIterator\.",         PendingDeprecationWarning, "gffutils", 186 )
   warnings.filterwarnings("ignore", "unclosed file <_io\.TextIOWrapper",  ResourceWarning,           "gffutils", 668 )
   test_gff_db = gffutils.create_db(   test_gff_file,
                                       dbfn                    = test_gff_db_file,
                                       force                   = True,     # overwrite previous testing db file
                                       merge_strategy          = 'error',
                                       keep_order              = False,    # True doesn't appear to maintain attribute order :-/  (and turning this off may be faster)
                                       sort_attribute_values   = False
                                       )
      
   for this_feature in test_gff_db.all_features():
      print(dumper.dump(this_feature))
