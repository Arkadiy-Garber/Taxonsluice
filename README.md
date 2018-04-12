# CleanMyOtus
Heuristic-algorithm for decontamination of OTU tables.

More detailed tutorial coming soon!

## Dependencies
-Python3

-BLAST

## Usage

This script takes as input the mothur-formatted tab-delimited OTU table (see sampleOTUtable.txt for specific format of this), as well as a tab-delimited file that maps each sample to a corresponding blank (see sample_blank_map.txt for specific format).

This script will run with only these two inputs, and output two files: 

-an OTU table with only the OTUs whose abundance in the sample-specific and non-specific-blanks is 10% or less of the abundance in the environmental samples.

-an OTU table containing sequences that were flagged due to presence within blanks

