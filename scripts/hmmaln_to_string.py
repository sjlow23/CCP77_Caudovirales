#!/usr/bin/python


import os
import sys
import argparse

class ParseHmm(object):

  def __init__(self):
    """Initialization."""
    pass
  
  
  def getAlignedMarker(self, hit_name, result_file):
       '''
       Parse the output of Hmmalign

       :param hit_name: gene name
       :param result_file: output file from Hmmalign
       '''
       hit_seq = None
       mask_seq = None
       with open(result_file,'r') as resf:
         for line in resf:
             splitline = line.split(" ", 1)
             if splitline[0] == hit_name.split(" ", 1)[0]:
               rsplitline = line.rsplit(" ", 1)
               hit_seq = rsplitline[-1]
               continue
             if line[0:len("#=GC RF")] == "#=GC RF":
               rsplitline = line.rsplit(" ", 1)
               mask_seq = rsplitline[-1]

         if mask_seq is None:
             raise Exception("Unable to get mask from hmm align result file")

         if hit_seq is None:
             return None

         aligned_marker = ""
         for pos in xrange(0, len(mask_seq)):
             if mask_seq[pos] != 'x':
               continue
             aligned_marker += hit_seq[pos]
         print aligned_marker

if __name__ == '__main__':

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('hitname')
  parser.add_argument('result_file' )

  args = parser.parse_args()

  try:
    hmmparser = ParseHmm()
    hmmparser.getAlignedMarker(args.hitname, args.result_file )
  except SystemExit:
    print "\nControlled exit resulting from an unrecoverable error or warning."
  except:
    print "\nUnexpected error:", sys.exc_info()[0]
    raise

