#!/usr/bin/python

"""
Reads a list of intervals (start, stop) and a maf. Produces a new maf
containing the blocks from the original that overlapped the intervals.

NOTE: See maf_extract_ranges_indexed.py which works better / faster for many
      use cases.

NOTE: chromosome/src information in the MAF is ignored by this variant.

usage: %prog interval_file refindex [options] < maf_file
   -m, --mincols=10: Minimum length (columns) required for alignment to be output
"""

import psyco_full

from bx.cookbook import doc_optparse

import bx.align.maf
from bx import intervals
import sys


def __main__():

    # Parse Command Line

    options, args = doc_optparse.parse( __doc__ )

    try:
        range_filename = args[ 0 ]
        refindex = int( args[ 1 ] )
        if options.mincols: mincols = int( options.mincols )
        else: mincols = 10
    except:
        doc_optparse.exit()

    # Load Intervals

    intersecter = intervals.Intersecter()
    for line in file( range_filename ):
        fields = line.split()
        intersecter.add_interval( intervals.Interval( int( fields[0] ), int( fields[1] ) ) )

    # Start MAF on stdout

    out = bx.align.maf.Writer( sys.stdout )

    # Iterate over input MAF

    for maf in bx.align.maf.Reader( sys.stdin ):
        ref = maf.components[ refindex ]
        # Find overlap with reference component
        intersections = intersecter.find( ref.get_forward_strand_start(), ref.get_forward_strand_end() )
        # Keep output maf ordered
        intersections.sort()
        # Write each intersecting block
        for interval in intersections: 
            start = max( interval.start, ref.get_forward_strand_start() )
            end = min( interval.end, ref.get_forward_strand_end() )
            sliced = maf.slice_by_component( refindex, start, end ) 
            good = True
            for c in sliced.components: 
                if c.size < 1: 
                    good = False
            if good and sliced.text_size > mincols: out.write( sliced )
         
    # Close output MAF

    out.close()

if __name__ == "__main__": __main__()
