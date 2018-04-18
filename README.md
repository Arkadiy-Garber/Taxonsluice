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

However, if you provide the script with the location of the SILVA database, and a fasta file containing the OTU sequences, an additional output will include an annotated summary of the OTUs that were flagged. This output will include the closest match to the flagged OTUs in SILVA, the study in which those matches originated, and the source of isolation.

### sample command (with the representative 16S sequences in FASTA format and SILVA database provided for classification of potentially-rare OTUs)
    python3 cleanmyotus.py -blank_map mySamplesToBlanks.txt -otu_table myOTUs.txt -seq_file myOTUs.fasta -silva_DB SILVA_128_SSURef_tax_silva.fasta -rare 5 -t 4 -silva_aln 10 -out_folder /path/to/output/directory/

### sample command (simple version without 16S or SILVA database provided)
    python3 cleanmyotus.py -blank_map mySamplesToBlanks.txt -otu_table myOTUs.txt -rare 5 -out_folder /path/to/output/directory/
