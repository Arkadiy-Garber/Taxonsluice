#!/usr/bin/env python3
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import urllib.request
import ssl

gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)

#TODO: make default setting for rare OTU definition.

# *************************************************************************
# ****************** Function Definitions *********************************
# *************************************************************************


def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x


def digitize(string):
    outStr = ''
    for i in string:
        try:
            int(i)
            outStr += str(i)
        except ValueError:
            pass
    return (int(outStr))


def derep(inList):
    outList = []
    for i in inList:
        if i not in outList:
            outList.append(i)
    return outList


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)


def deAligner(inList):
    outList = []
    for i in inList:
        if i != "-":
            outList.append(i)
    outStr = "".join(outList)
    return outStr


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


parser = argparse.ArgumentParser(
    prog="cleanmyotus.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Gustavo Ramírez^1 and Arkadiy Garber^2; 
    1^University of Rhode Island, Graduate School of Oceanography
    2^University of Southern California, Earth Sciences
    Please send comments and inquiries to arkg@udel.edu

    *******************************************************
    '''))

parser.add_argument('-blank_map', type=str, help='Tab-delimited file denoting which blanks correspond to which samples')

parser.add_argument('-otu_table', type=str,
                    help="Mothur-formatted tab-delimited OTU table, with each row corresponding "
                         "to either a sample or blank. The sample or blank names in this table must match the "
                         "names provided in the blank_map")

parser.add_argument('-seq_file', type=str, help="FASTA file contained the 16S sequences. Header must contain the Otu #"
                                                "(Optional, if you want to classify your potential contaminants)")

parser.add_argument('-silva_DB', type=str, help="SILVA database (provide this if you want to classify your potential "
                                                "contaminants)", default="NA")

parser.add_argument('-rare', type=int,
                    help="remove OTUs that are represented by less than this number of sequences (default = 2)",
                    default=2)

parser.add_argument('-t', type=int, help="number of threads to use for BLASTn (default = 1)", default=1)

parser.add_argument('-silva_aln', type=str,
                    help="how many alignments from the SILVA reference database to keep in the final Flagged OTUs summary file",
                    default="10")

parser.add_argument('-out_folder', type=str,
                    help="Directory to which output files will be written (default = current working directory)",
                    default="out")

args = parser.parse_args()

# *************************************************************************
# ************************** Format Sensing *******************************
# *************************************************************************


blankList = []
BlankDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in BlankDict:
    ls = (i.rstrip().split("\t"))
    if ls[1] not in blankList:
        blankList.append(ls[1])
numBlanks = len(blankList)

# *************************************************************************
# ************ Contaminant and Rare OTU Identification ********************
# *************************************************************************


print("Removing contaminants and \'rare\' OTUs (as defined by the user to be any OTU that is represented by less than "
      + str(args.rare) + " sequences)")
BlankDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
blanks = open(args.blank_map, "r")
for i in blanks:
    ls = (i.rstrip().split("\t"))
    BlankDict[ls[0]]["id"] = 'env'
    BlankDict[ls[0]]["connection"] = ls[1]
    BlankDict[ls[1]]["id"] = 'blank'
    # BlankDict[ls[1]]["connection"] = ls[0]

mothur_output = open(args.otu_table, "r")
OTUdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
counter = 0
for line in mothur_output:
    ls = (line.rstrip().split("\t"))
    if counter == 0:
        for i in ls:
            if re.findall(r'Otu\d', i):
                break
            else:
                counter += 1
        firstOTU = counter
    else:
        sample = ls[1]
        if sample in BlankDict.keys():
            otu = 1
            OTU = stabilityCounter(otu)
            for OTUabund in ls[firstOTU:]:
                BlankDict[sample]["Otu" + str(OTU)] = OTUabund
                otu += 1
                OTU = stabilityCounter(otu)

contaminants = []
rareOTUs = []
count = 0
for i in BlankDict.keys():
    if BlankDict[i]["id"] == "env":
        sample = (i)
        sampleSpecificBlank = (BlankDict[i]["connection"])
        for j in range(1, otu):
            Otu = ("Otu" + stabilityCounter(j))
            sampleOTUabund = BlankDict[sample][Otu]
            blankOTUabund = BlankDict[sampleSpecificBlank][Otu]
            if int(sampleOTUabund) == 0 and int(blankOTUabund) > 0:
                counter = 0
                otherSampleOTUabundTotal = 0
                for otherSamples in BlankDict.keys():
                    if otherSamples != i and BlankDict[otherSamples]["id"] == "env":
                        otherSampleOTUabund = BlankDict[otherSamples][Otu]
                        otherSampleOTUabundTotal += int(otherSampleOTUabund)
                        if int(otherSampleOTUabund) != 0:
                            counter += 1
                if otherSampleOTUabundTotal == 0:
                    count += 1
                    contaminants.append(Otu)
            elif int(sampleOTUabund) < int(args.rare):
                otherSampleOTUabundTotal = int(sampleOTUabund)
                for otherSamples in BlankDict.keys():
                    if otherSamples != i and BlankDict[otherSamples]["id"] == "env":
                        otherSampleOTUabund = BlankDict[otherSamples][Otu]
                        otherSampleOTUabundTotal += int(otherSampleOTUabund)
                if otherSampleOTUabundTotal < int(args.rare):
                    rareOTUs.append(Otu)
contaminants = (derep(sorted(contaminants)))
rareOTUs = (derep(sorted(rareOTUs)))

# *************************************************************************
# ****************** Sample-Specific Blank Identification *****************
# *************************************************************************


BlankDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
blanks = open(args.blank_map, "r")
for i in blanks:
    ls = (i.rstrip().split("\t"))
    BlankDict[ls[0]]["id"] = 'env'
    BlankDict[ls[0]]["connection"] = ls[1]
    BlankDict[ls[1]]["id"] = 'blank'

# ****************************************************************************
# ********************* OTU Table Transformation *****************************
# ****************************************************************************


count = 0
mothur_output = open(args.otu_table, "r")
OTUdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
counter = 0
for line in mothur_output:
    ls = (line.rstrip().split("\t"))
    if counter == 0:
        for i in ls:
            if re.findall(r'Otu\d', i):
                break
            else:
                counter += 1
        firstOTU = counter
    else:
        sample = ls[1]
        if sample in BlankDict.keys():
            otu = 1
            OTU = stabilityCounter(otu)
            for OTUabund in ls[firstOTU:]:
                Otu = ("Otu" + str(OTU))
                if Otu in contaminants or Otu in rareOTUs:
                    otu += 1
                    OTU = stabilityCounter(otu)
                else:
                    BlankDict[sample][Otu] = OTUabund
                    otu += 1
                    OTU = stabilityCounter(otu)

# ***********************************************************************
# ******************** Algorithm Implementation *************************
# ***********************************************************************


print("Identifying potential contaminants...")
cleanOTUs = []
discardedOTUs = []
flaggedOTUs = []
for sample in BlankDict.keys():
    blank = (BlankDict[sample]["connection"])
    if blank != "EMPTY":
        for i in range(1, int(otu)):
            i = stabilityCounter(i)
            Otu = ("Otu" + str(i))
            if Otu not in contaminants and Otu not in rareOTUs:
                env_abund = BlankDict[sample]["Otu" + str(i)]  # abundance of OTU in sample.
                blank_abund = BlankDict[blank]["Otu" + str(i)]  # abundance of OTU in sample-specific blank
                if int(blank_abund) == 0 and int(env_abund) > 0:  # If the OTU is not present in sample-specific blank
                    if numBlanks > 1:
                        counter = 0  # setting memory variable
                        for sample2 in BlankDict.keys():
                            if BlankDict[sample2]["id"] == "blank" and sample2 != sample:  # Looking for OTU abundances in
                                # non-sample specific blanks
                                nonSpecificBlank_abund = BlankDict[sample2][
                                    "Otu" + str(i)]  # abundance of OTU in non-sample
                                # specific blank
                                if int(nonSpecificBlank_abund) == 0:  # If OTU is not present in this blank
                                    pass  # move on to next blank
                                else:
                                    if int(nonSpecificBlank_abund) * 10 >= int(env_abund):  # If the OTU is more than
                                        # 10% abundant in blank, compared with the sample
                                        if i not in flaggedOTUs:
                                            flaggedOTUs.append("Otu" + i)
                                            # Discard OTU if it has not been discarded already.
                                        counter += 1
                                        break
                        if counter == 0:  # if all the non-sample specific blanks were looked at and OTU not discarded
                            if i not in cleanOTUs:
                                cleanOTUs.append("Otu" + i)  # OTU is clean
                            counter = 0  # re-setting memory variable
                    else:
                        cleanOTUs.append("Otu" + i)
                else:  # IF OTU is present in sample-specific blank
                    counter = 0  # setting memory variable
                    if int(blank_abund) * 10 <= int(env_abund):  # If OTU is 10x more abundant in sample than in blank
                        if numBlanks > 1:
                            for sample2 in BlankDict.keys():  # Looking for OTU in non-sample specific blanks
                                if BlankDict[sample2]["id"] == "blank" and sample != sample:
                                    nonSpecificBlank_abund = BlankDict[sample2]["Otu" + str(i)]  # abundance of OTU in
                                    # non-sample specific blank
                                    if int(nonSpecificBlank_abund) == 0:  # if OTU not present in this blank
                                        pass  # do nothing
                                    else:  # if it is...
                                        # Is the OTU more than 10% abundant in the blank, compared to the sample
                                        if int(nonSpecificBlank_abund) * 10 >= int(env_abund):
                                            if i not in flaggedOTUs:  # if it is, discard OTU
                                                flaggedOTUs.append("Otu" + i)
                                            break

                            if counter == 0:  # if all the non-sample specific blanks were looked at and OTU not discarded
                                if i not in cleanOTUs:
                                    cleanOTUs.append("Otu" + i)  # OTU is clean
                                counter = 0  # re-setting memory variable
                        else:
                            cleanOTUs.append("Otu" + i)
                    else:
                        if i not in flaggedOTUs:
                            flaggedOTUs.append("Otu" + i)
            else:
                discardedOTUs.append("Otu" + i)

CleanOTUs_final = []
RareOTUs_final = []
for i in cleanOTUs:
    if i not in flaggedOTUs:
        CleanOTUs_final.append(i)

for i in rareOTUs:
    if i not in contaminants:
        RareOTUs_final.append(i)

finalFlagged = sorted(derep(flaggedOTUs))
finalClean = sorted(derep(CleanOTUs_final))
finalRare = sorted(derep(RareOTUs_final))
finalCont = sorted(derep(contaminants))

# ***************************************************************************
# **************************** Writing outfiles *****************************
# ***************************************************************************


print("Preparing file with flagged OTUs")
outFlagged = open(args.out_folder + "/flaggedOTUs.csv", "w")
outFlagged.write("label" + "," + "Group" + "," + "numOtus" + ",")
for i in finalFlagged:
    outFlagged.write(i + ",")
outFlagged.write("\n")

for i in BlankDict.keys():
    if BlankDict[i]["id"] == "env":
        sampleOTUs = BlankDict[i]
        outFlagged.write("0.03" + "," + i + "," + str(len(finalFlagged)) + ",")
        lastOTU = (lastItem(finalFlagged))
        length = (digitize(lastOTU))
        for j in range(1, length + 1):
            Otu = "Otu" + str(stabilityCounter(j))
            if Otu in finalFlagged:
                outFlagged.write(BlankDict[i][Otu] + ",")
        outFlagged.write("\n")
outFlagged.close()

print("preparing file with non-flagged OTUs")
outclean = open(args.out_folder + "/cleanOTUs.csv", "w")
outclean.write("label" + "," + "Group" + "," + "numOtus" + ",")
for i in finalClean:
    outclean.write(i + ",")
outclean.write("\n")

for i in BlankDict.keys():
    if BlankDict[i]["id"] == "env":
        sampleOTUs = BlankDict[i]
        outclean.write("0.03" + "," + i + "," + str(len(finalClean)) + ",")
        lastOTU = (lastItem(finalClean))
        length = (digitize(lastOTU))
        for j in range(1, length + 1):
            Otu = "Otu" + str(stabilityCounter(j))
            if Otu in finalClean:
                outclean.write(BlankDict[i][Otu] + ",")
        outclean.write("\n")
outclean.close()

# ***************************************************************************
# ************************ Sequence identification **************************
# ***************************************************************************

if args.silva_DB != "NA":
    seqs = open(args.seq_file, "r")
    seqs = fasta(seqs)

    print("Writing flagged OTUs to FASTA file")
    outfile = open(args.out_folder + "/flaggedOTUs.fasta", "w")
    for i in finalFlagged:
        for j in seqs.keys():
            if re.findall(i, j):
                outfile.write(">" + j + "\n")
                outfile.write(seqs[j] + "\n")
    outfile.close()

    print("Building BLAST database file from the SILVA database...")
    os.system("makeblastdb -dbtype nucl -in " + args.silva_DB +
              " -out " + args.silva_DB)

    print("comparing flagged OTUs against the SILVA database...")
    os.system("blastn -query " + args.out_folder + "/flaggedOTUs.fasta" + " -out " + args.out_folder + "/OTUblast.txt" +
              " -outfmt 6 -evalue 1E-50 -db " + args.silva_DB)
    print("Done with BLAST. Now parsing BLAST output file")
    os.system("rm " + args.silva_DB + "*.nhr")
    os.system("rm " + args.silva_DB + "*.nin")
    os.system("rm " + args.silva_DB + "*.nsq")

    nameDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    silva = open(args.silva_DB, "r")
    for i in silva:
        if re.match(r'>', i):
            id = i.split(" ")[0][1:]
            id = id.split(".")[0]
            nameDict[id] = i[1:]

    blastDict = defaultdict(lambda: defaultdict(list))
    blast = open(args.out_folder + "/OTUblast.txt", "r")
    count = 0
    for i in blast:
        ls = i.rstrip().split("\t")
        otu = ls[0]
        match = ls[1].split(".")[0]
        identity = ls[2]
        cov = ls[3]
        e = ls[10]
        if count <= int(args.silva_aln):
            blastDict[otu]["match"].append(match)
            blastDict[otu]["identity"].append(identity)
            blastDict[otu]["cov"].append(cov)
            blastDict[otu]["e"].append(e)
            count += 1
        count = 0

    print("Identifying potential contaminants and writing summary file")
    outfile = open(args.out_folder + "/flaggedOtus_annotated.csv", "w")
    outfile.write(
        "OTU" + "," + "silva_match" + "," + "perc_id" + "," + "perc_cov" + "," + "e_value" + "," + "source" + "\n")
    for i in blastDict.keys():
        for j in range(0, 10):
            match = blastDict[i]["match"][j]
            identity = blastDict[i]["identity"][j]
            cov = blastDict[i]["cov"][j]
            e = blastDict[i]["e"][j]
            head = (nameDict[match].rstrip())
            fp = urllib.request.urlopen("https://www.ebi.ac.uk/ena/data/view/%s&display=text" % match, context=gcontext)
            mybytes = fp.read()
            mystr = mybytes.decode("utf8")
            fp.close()
            lines = mystr.split('\n')
            count = 0
            string = ''
            for j in lines:
                if re.findall(r'isolation_source', j) or re.findall(r'tissue_type', j):
                    count += 1
                    string = j.split("\"")[1]
                    string = replace(string, [","], ";")
                elif count > 0:
                    if not re.findall(r'/', j):
                        string += (" " + j.split("                   ")[1])
                        string = string[0:len(string) - 1]
                        string = replace(string, [","], ";")
                    else:
                        count = 0
            outfile.write(i + "," + head + "," + identity + "," + cov + "," + e + "," + string + "\n")
    os.system("rm " + args.out_folder + "/OTUblast.txt")
    os.system("rm " + args.out_folder + "/flaggedOTUs.fasta")

print("All done! Thank you for using taxonsluice.")
