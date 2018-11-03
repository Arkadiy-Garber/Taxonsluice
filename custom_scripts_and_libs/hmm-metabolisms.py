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


def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


def deExt(string):
    ls = string.split(".")
    ls2 = ls[0:len(ls)-1]
    outstring = "".join(ls2)
    return outstring


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:].split(" ")[0]
                seq = ''
            else:
                header = i[1:].split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


parser = argparse.ArgumentParser(
    prog="hmm-metabolisms.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Arkadiy Garber;
    University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu

    *******************************************************
    '''))

parser.add_argument('-hmm_dir', type=str, help='directory of HMMs')
parser.add_argument('-bin_dir', type=str, help="directory of bins")
parser.add_argument('-bin_ext', type=str, help="extension for bins (do not include the period)")
parser.add_argument('-out', type=str, help="basename of output file", default="out")
parser.add_argument('-m', type=str, help="location HMMs meta-data file (hmm-meta.csv). This is optional, only if you "
                                         "want to filter results by the suggested \'trusted\' bitscore thresholds for "
                                         "each HMM", default="NA")
parser.add_argument('-mode', type=str, help="would you like all hits to an HMM, or just one (all/one)", default="all")
args = parser.parse_args()


# binDir = "/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/magnetite_reducer/bin"
# hmmDir = "/Users/arkadiygarber/Desktop/metabolic-hmms-banfield-ag/"
bins = os.listdir(args.bin_dir)
HMMs = os.listdir(args.hmm_dir)


bits = open(args.m, "r")
bitDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
print("\nreading in HMM bitscore cut-offs...")
for i in bits:
    ls = i.rstrip().split("\t")
    bitDict[ls[0]]["bit"] = ls[1]
    bitDict[ls[0]]["process"] = ls[2]
    bitDict[ls[0]]["element"] = ls[3]

redundancyDict = defaultdict(list)
out = open(args.out + ".csv", "w")
out.write("bin" + "," + "gene" + "," + "process" + "," + "substrate" + "," + "ORF" + "," + "evalue" + "," + "bitscore" + "," + "sequence" + "\n")
for i in bins:
    HMMdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    if not re.match(r'^\.', i) and i != (args.out + ".csv") and lastItem(i.split(".")) == args.bin_ext:
        print("analyzing " + i)
        fastaFile = open(args.bin_dir + "/" + i, "r")
        fastaFile = fasta(fastaFile)
        os.system("mkdir " + args.bin_dir + "/" + i + "-HMM")
        for hmm in HMMs:
            if not re.match(r'^\.', hmm):
                os.system("hmmsearch "
                          "--tblout " + args.bin_dir + "/" + i + "-HMM/" + i + "__" + hmm +
                          " -o " + args.bin_dir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt " +
                          args.hmm_dir + "/" + hmm + " " +
                          args.bin_dir + "/" + i)
                os.system("rm " + args.bin_dir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt")

        HMMresults = os.listdir(args.bin_dir + "/" + i + "-HMM")
        for result in HMMresults:
            if not re.match(r'^\.', result):
                file = open(args.bin_dir + "/" + i + "-HMM/" + result)
                for line in file:
                    line = line.rstrip()
                    if not re.match(r'^#', line):
                        ls = (delim(line))
                        orf = ls[0]
                        gene = ls[2]
                        evalue = ls[4]
                        bitscore = ls[5]
                        hmmFileName = (result.split("__")[1])
                        hmmFileName = hmmFileName.split(".hm")[0]
                        print(hmmFileName)
                        hmmName = deExt(result.split("__")[1])
                        threshold = bitDict[hmmFileName]["bit"]
                        process = bitDict[hmmFileName]["process"]
                        element = bitDict[hmmFileName]["element"]
                        if float(threshold) > 0:
                            if float(bitscore) > float(threshold):
                                if orf not in HMMdict.keys():
                                    HMMdict[orf]["hmm"] = hmmName
                                    HMMdict[orf]["evalue"] = evalue
                                    HMMdict[orf]["bitscore"] = bitscore
                                    HMMdict[orf]["process"] = process
                                    HMMdict[orf]["element"] = element
                                else:
                                    if float(bitscore) > float(HMMdict[orf]["bitscore"]):
                                        HMMdict[orf]["hmm"] = hmmName
                                        HMMdict[orf]["evalue"] = evalue
                                        HMMdict[orf]["bitscore"] = bitscore
                                        HMMdict[orf]["process"] = process
                                        HMMdict[orf]["element"] = element
                                    else:
                                        pass
                        else:
                            if float(evalue) <= float(1E-30):
                                if orf not in HMMdict.keys():
                                    HMMdict[orf]["hmm"] = hmmName
                                    HMMdict[orf]["evalue"] = evalue
                                    HMMdict[orf]["bitscore"] = bitscore
                                    HMMdict[orf]["process"] = process
                                    HMMdict[orf]["element"] = element
                                else:
                                    if float(bitscore) > float(HMMdict[orf]["bitscore"]):
                                        HMMdict[orf]["hmm"] = hmmName
                                        HMMdict[orf]["evalue"] = evalue
                                        HMMdict[orf]["bitscore"] = bitscore
                                        HMMdict[orf]["process"] = process
                                        HMMdict[orf]["element"] = element
                                    else:
                                        pass

        if args.mode == "one":
            for key in HMMdict.keys():
                if HMMdict[key]["hmm"] not in redundancyDict[i]:

                    redundancyDict[i].append(HMMdict[key]["hmm"])

                    out.write(i + "," + HMMdict[key]["hmm"] + "," + HMMdict[key]["process"] + "," + HMMdict[key]["element"] +
                              "," + key + "," + HMMdict[key]["evalue"] + "," + HMMdict[key]["bitscore"] + "," + fastaFile[key]
                              + "\n")
            print(redundancyDict[i])
            os.system("rm -r " + args.bin_dir + "/" + i + "-HMM")

        if args.mode == "all":
            for key in HMMdict.keys():
                out.write(
                    i + "," + HMMdict[key]["hmm"] + "," + HMMdict[key]["process"] + "," + HMMdict[key]["element"] +
                    "," + key + "," + HMMdict[key]["evalue"] + "," + HMMdict[key]["bitscore"] + "," + fastaFile[key]
                    + "\n")
            os.system("rm -r " + args.bin_dir + "/" + i + "-HMM")

print("Finished")