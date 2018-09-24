# taxonsluice
Heuristic algorithm for decontamination of OTU tables: The currently accepted input to this program is a mothur-formated OTU table. However, we are currently working on making this tool compatibale with QIIME output and amplicon sequence variants (ASVs). This tool's decontaminanting power is dependent on the presence of multiple control blanks that are specific to different samples or groups of samples within a study. The algorithm compared the relative abundance of an OTU in a blank relative to a sample associated with that blank; next, the relative abundance of that OTU is compared across all other blanks and their associated samples. Based on this two-tier comparison strategy, a decision is made to either keep an OTU in the dataset unchecked, or to flag it; flagged OTUs are then compared against the SILVA database using BLAST; an output is then provided summarizing the closest matches to each OTU that was flagged; this output includes the sequence identity, taxanomy, isolation source (if available), and study where the closest-matching amplicon was sequenced.

The user then makes a decision about which, if any, flagged OTUs to allow back into the full dataset. To inform this decision, the user may look at the relative abundance of the OTU in the blanks vs. samples (two OTU tables are also provided in the output: one with only the unflagged sequences, and one with only the flagged sequences); the information obtainined from the closest match in SILVA may also be used. However, we recognize the limitation of using reference sequences in SILVA to inform any decisions regarding potential contaminants, and strongly caution users to take that reference information with a grain of salt.


More detailed tutorial coming soon! (Please see the 'test' folder for the sample dataset that you can run this program on)

## Dependencies
-Python3

-BLAST

## Usage

This script takes as input the mothur-formatted tab-delimited OTU table (see sampleOTUtable.txt for specific format of this), as well as a tab-delimited file that maps each sample to a corresponding blank (see sample_blank_map.txt for specific format).

This script will run with only these two inputs, and output two files: 

  -an OTU table with only the OTUs whose abundance in the sample-specific and non-specific-blanks is 10% or less of the abundance in the environmental samples.

  -an OTU table containing sequences that were flagged due to presence within blanks

However, if you provide the script with the location of the SILVA database, and a fasta file containing the OTU sequences, an additional output will include an annotated summary of the OTUs that were flagged. This output will include the closest match to the flagged OTUs in SILVA, the study in which those matches originated, and the source of isolation.

The SILVA database can be downloaded from https://www.arb-silva.de/no_cache/download/archive/current/Exports/. The file you should get is "SILVA_132_SSUParc_tax_silva.fasta.gz". As of April 2018, release 132 is the latest release.

### sample command (with the representative 16S sequences in FASTA format and SILVA database provided for classification of potentially-rare OTUs)
    python3 taxonsluice.py -blank_map mySamplesToBlanks.txt -otu_table myOTUs.txt -seq_file myOTUs.fasta -silva_DB SILVA_128_SSURef_tax_silva.fasta -rare 5 -t 4 -silva_aln 10 -out_folder /path/to/output/directory/

### sample command (simple version without 16S or SILVA database provided)
    python3 taxonsluice.py -blank_map mySamplesToBlanks.txt -otu_table myOTUs.txt -rare 5 -out_folder /path/to/output/directory/
