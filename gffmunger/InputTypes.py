import os
import argparse

class InputTypes:

   def min_length_3(whatever):
      if len(str(whatever)) < 3:
         raise argparse.ArgumentTypeError("must contain at least 3 chars")
         return False
      return whatever
