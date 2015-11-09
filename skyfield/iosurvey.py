"""Routines to explore available datasets
"""

# load low-level file tools
from jplephem.spk import SPK
from jplephem.daf import DAF
# standard tools
from StringIO import StringIO # to create in-memory file
import urllib2 # to retrieve data

def get_summary(url, spk=True):
    ''' simple function to retrieve the header of a BSP file and return SPK object'''
    # connect to file at URL
    bspurl = urllib2.urlopen(url)
    # retrieve the "tip" of a file at URL
    bsptip = bspurl.read(10**5) # first 100kB
    # save data in fake file object (in-memory)
    bspstr = StringIO(bsptip)
    # load into DAF object
    daf = DAF(bspstr)
    # return either SPK or DAF object
    if spk:
      # make a SPK object
      spk = SPK(daf)
      # return representation 
      return spk
    else:
      # return representation 
      return daf