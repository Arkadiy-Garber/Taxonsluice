# Author: Arkadiy Garber
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import urllib.request
import ssl
from urllib.error import HTTPError


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


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
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


# *******************************************************************************************************************
# *******************************************************************************************************************
# *******************************************************************************************************************


# blast = open("/Users/arkadiygarber/Desktop/flocs.blast", "r")
# scaffs = open("/Users/arkadiygarber/Desktop/scaffolds1362AB-renamed.fasta", "r")
# scaffs = fasta(scaffs)
# flocs = open("/Users/arkadiygarber/Desktop/JDF16S.fna", "r")
# flocs = fasta(flocs)
#
# blastDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
# for i in blast:
#     ls = i.rstrip().split("\t")
#     if ls[1] not in blastDict.keys():
#         blastDict[ls[1]]["otu"] = ls[0]
#         blastDict[ls[1]]["id"] = ls[2]
#         blastDict[ls[1]]["cov"] = float(ls[3])
#
# out = open("/Users/arkadiygarber/Desktop/FLOCSscaffs.fasta", "w")
# for i in scaffs.keys():
#     if i in blastDict.keys():
#         out.write(">" + i + "\n")
#         out.write(scaffs[i] + "\n")


# *******************************************************************************************************************
# **************************** ANALYZING RESULTS FROM GHOSTKOALA ****************************************************
# *******************************************************************************************************************

flocs = open("/Users/arkadiygarber/Desktop/JdF2018_last_efforts/KO/FLOCS_ko_definition.txt", "r")
flocsKO = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
flocsKOcounts = defaultdict(list)
for i in flocs:
    ls = (i.rstrip().split("\t"))
    if len(ls) > 1:
        ko = ls[1]
        if ko != "":
            flocsKO[ko] = ls[2]
            flocsKOcounts[ko].append(ls[2])

flocscopyDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in flocsKO.keys():
    flocscopyDict[flocsKO[i]] = len(flocsKOcounts[i])


plankton = open("/Users/arkadiygarber/Desktop/JdF2018_last_efforts/KO/PLANKTON_ko_definition.txt", "r")
planktonKO = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
planktonKOcounts = defaultdict(list)
for i in plankton:
    ls = (i.rstrip().split("\t"))
    if len(ls) > 1:
        ko = ls[1]
        if ko != "":
            planktonKO[ko] = ls[2]
            planktonKOcounts[ko].append(ls[2])

PlanktoncopyDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in planktonKO.keys():
    PlanktoncopyDict[planktonKO[i]] = len(planktonKOcounts[i])

# *******************************************************************************************************************
# *******************************************************************************************************************
# *******************************************************************************************************************


db = open("/Users/arkadiygarber/GhostKoalaParser/KO_Orthology_ko00001.txt", "r")
dbDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in db:
    ls = (i.rstrip().split("\t"))
    ko = (ls[3].split("  ")[0])
    dbDict[ko]["family"] = ls[1]
    dbDict[ko]["path"] = ls[2]

Dictionary = defaultdict(list)
for i in flocsKO.keys():
    if i not in planktonKO.keys():
        Dictionary[dbDict[i]["path"]].append(flocsKO[i])


# *******************************************************************************************************************
# ********************** CREATING A BED FILE FROM ORF CALLS ON JDF ASSEMBLY *****************************************
# *******************************************************************************************************************

out = open("/Users/arkadiygarber/Desktop/genes_exclusive_to_FLOCS.csv", "w")
out.write("KEGG category" + "," + "gene" + "," + "gene copies in FLOCS" + "," + "gene copies in PLANKTON" + "\n")
for i in Dictionary.keys():
    for j in Dictionary[i]:
        out.write(replace(i, [","], ";") + "," + replace(j, [","], ";") + "," + str(flocscopyDict[j]) + "," + str(PlanktoncopyDict[j]) + "\n")

prots = open("/Users/arkadiygarber/Desktop/JdF2018_last_efforts/scaffolds1362AB-renamed-proteins.faa", "r")
out = open("/Users/arkadiygarber/Desktop/JdF2018_last_efforts/scaffolds1362AB-renamed-proteins.bed", "w")
for i in prots:
    if re.match(r'^>', i):
        i = i.rstrip()[1:]
        ls = i.split(" # ")
        print(i)
        contig = (ls[0].split("_")[0])
        start = ls[1]
        stop = ls[2]
        out.write(str(contig) + "\t" + str(start) + "\t" + str(stop) + "\n")

# *******************************************************************************************************************
# ****************** PULLING OUT CLEAN OTUS FROM  OTU TABLE OUTPUT OF TAXONSLUICE ***********************************
# *******************************************************************************************************************

otu = open("/Volumes/Seagate_Backup_Plus_Drive/Ongoing Research Projects/JdF2018/Algorithm/cleanOTUs.csv", "r")
otuList = []
count = 0
for i in otu:
    if count == 0:
        ls = i.rstrip().split(",")
        for j in ls:
            if re.findall(r'Otu0', j):
                otuList.append(j)
        count += 1
    else:
        break

amplicons = open("/Volumes/Seagate_Backup_Plus_Drive/Ongoing Research Projects/JdF2018/Algorithm/JDF16S.fna", "r")
amplicons = fasta(amplicons)
out = open("/Users/arkadiygarber/Desktop/CleanOTUs.fna", 'w')
for i in amplicons.keys():
    otu = i.split("_")[7].split("|")[0]
    if otu in otuList:
        out.write(">" + i + "\n")
        out.write(amplicons[i] + "\n")

# *******************************************************************************************************************
# ******************* ANALYZING BLAST RESULTS OF JDF rRNA READS AGAINST AMPLICONS **********************************
# *******************************************************************************************************************

blast = open("/Users/arkadiygarber/Desktop/JDF16Scombined.blast", "r")
Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in blast:
    ls = i.rstrip().split("\t")
    read = ls[0]
    otu = ls[1]
    otu = otu.split("_")[7].split("|")[0]
    id = ls[2]
    cov = ls[3]
    pos1 = ls[8]
    pos2 = ls[9]
    start = sorted([pos1, pos2])[0]
    end = sorted([pos1, pos2])[1]
    if int(cov) > 100 and int(start) < 242:
        if read not in Dict.keys():
            Dict[read]["otu"] = otu
            Dict[read]["id"] = id
            Dict[read]["cov"] = cov
            Dict[read]["start"] = start
            Dict[read]["end"] = end
        else:
            if float(id) > float(Dict[read]["id"]):
                Dict[read]["otu"] = otu
                Dict[read]["id"] = id
                Dict[read]["cov"] = cov
                Dict[read]["start"] = start
                Dict[read]["end"] = end
            else:
                if int(cov) > int(Dict[read]["cov"]):
                    Dict[read]["otu"] = otu
                    Dict[read]["id"] = id
                    Dict[read]["cov"] = cov
                    Dict[read]["start"] = start
                    Dict[read]["end"] = end
                else:
                    pass

count = 0
OTUdict = defaultdict(list)
for i in Dict.keys():
    OTUdict[Dict[i]["otu"]].append(i)
    count += 1
print(count)

out = open("/Users/arkadiygarber/Desktop/OTUabund_inFluids.csv", "w")
out.write("OTU" + "," + "rel_abund" + "\n")
for i in sorted(OTUdict.keys()):
    print(i)
    prop = (len(OTUdict[i])/count)
    print(prop)
    out .write(i + "," + str(prop) + "\n")
    perc = prop*100
    perc = ("%.3f" % perc + "%")
    print("")

for i in range(1, 10):
    print(i)


# **************************************************************************************
# **************************************************************************************
# **************************************************************************************


# meta = open("/Users/arkadiygarber/Desktop/JdF2018_last_efforts/Supplemental/JdF_Bin_FLOCS_summary-ArchaeaRemoved.csv")
# metaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
# for i in meta:
#     ls = i.rstrip().split(",")
#     if ls[5] != "plankton":
#         metaDict[ls[0]] = ls[5]
#
#
# silva = open("/Users/arkadiygarber/Desktop/OTU_diversity_in_Fluids/SILVA_132_SSURef_Nr99_tax_silva_trunc.ge1200bp.le2000bp.0.97.fixed.fasta")
# silvaDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
# for i in silva:
#     if re.match(r'^>', i):
#         ls = i.rstrip().split(" ")
#         silvaDict[ls[0][1:]] = ls[1]
#
# blast = open("/Users/arkadiygarber/Desktop/OTU_diversity_in_Fluids/ssu_silva.blast", "r")
# binDict = defaultdict(lambda: defaultdict(lambda: 'no_match'))
# for i in blast:
#     ls = i.rstrip().split("\t")
#     if ls[0] not in binDict.keys():
#         binDict[ls[0]]["name"] = silvaDict[ls[1]]
#         binDict[ls[0]]["id"] = ls[2]
#         binDict[ls[0]]["cov"] = ls[3]


# **************************************************************************************
# **************************************************************************************
# **************************************************************************************

# cov = open("/Users/arkadiygarber/Desktop/OTU_diversity_in_Fluids/ssu.coverage", "r")
# count = 0
# for i in cov:
#     ls = i.rstrip().split("\t")
#     if ls[2] != "totalAvgDepth":
#         count += float(ls[2])
#
# out = open("/Users/arkadiygarber/Desktop/JDF_TAXA_reconstructed16_COV.csv", "w")
# out.write("16S" + "," + "relative_proportion_in_fluids" + "," + "SILVA_match" + "," + "similarity" + "\n")
# cov = open("/Users/arkadiygarber/Desktop/OTU_diversity_in_Fluids/ssu.coverage", "r")
# for i in cov:
#     ls = i.rstrip().split("\t")
#     name = ls[0]
#     if ls[2] != "totalAvgDepth" and int(ls[1]) > 200:
#         if re.findall(r'FLOCS', ls[0]):
#             bin = (ls[0].split("&&")[0])
#             out.write(metaDict[bin] + "," + str(float(ls[2])/float(count)) + "," + str(binDict[name]["name"]) + "," + str(binDict[name]["id"]) +  "\n")
#         else:
#             contig = ls[0].split("&&")
#             out.write("BG" + "&&" + str(contig[1]) + "," + str(float(ls[2])/float(count)) + "," + str(binDict[name]["name"]) + "," + str(binDict[name]["id"]) + "\n")


# **************************************************************************************
# **************************************************************************************
# **************************************************************************************


# count = 0
# idx = open("/Users/arkadiygarber/Desktop/OTU_diversity_in_Fluids/ssu.idxstats")
# for i in idx:
#     ls = i.rstrip().split("\t")
#     if ls[0] != "*" and int(ls[1]) > 200:
#         correctedReads = float(float(ls[2]) / float(ls[1])) * 1000
#         count += correctedReads
#
# out = open("/Users/arkadiygarber/Desktop/OTU_diversity_in_Fluids/JDF_TAXA_reconstructed16_READS.csv", "w")
# out.write("16S" + "," + "relative_proportion_in_fluids" + "," + "SILVA_match" + "," + "similarity" + "\n")
# idx = open("/Users/arkadiygarber/Desktop/OTU_diversity_in_Fluids/ssu.idxstats")
# for i in idx:
#     ls = i.rstrip().split("\t")
#     if ls[0] != "*" and int(ls[1]) > 200:
#         correctedReads = float(float(ls[2])/float(ls[1]))*1000
#         if re.findall(r'FLOCS', ls[0]):
#             bin = (ls[0].split("&&")[0])
#             out.write(metaDict[bin] + "," + str(correctedReads/count)+ "," + str(binDict[ls[0]]["name"]) + "," + str(binDict[ls[0]]["id"]) + "\n")
#         else:
#             contig = ls[0].split("&&")
#             out.write("BG" + "&&" + str(contig[1]) + "," + str(correctedReads/count) + "," + str(binDict[ls[0]]["name"]) + "," + str(binDict[ls[0]]["id"]) + "\n")


# **************************************************************************************
# *** ANALYZING RESULTS OF AMPLICON COMPARISONS AGAINST JDF AND NP METAGENOME READS ****
# **************************************************************************************


NPreadsBlast = open('/Users/arkadiygarber/Desktop/NPss_FLOCS.blast', "r")
BlastDictNP = defaultdict(list)
for i in NPreadsBlast:
    ls = i.rstrip().split("\t")
    BlastDictNP[ls[1]].append(ls[0])

JDFreadsBlast = open('/Users/arkadiygarber/Desktop/JdF_FLOCS.blast', "r")
BlastDictJDF = defaultdict(list)
for i in JDFreadsBlast:
    ls = i.rstrip().split("\t")
    BlastDictJDF[ls[1]].append(ls[0])


# otu = open("/Users/arkadiygarber/Desktop/OTUtable.csv", "r")
# OTUdict = defaultdict(lambda: defaultdict(lambda: 'no_match'))
# counter = 0
# for i in otu:
#     ls = i.rstrip().split("\t")
#     print(ls)
#     count = 0
#     if counter == 0:
#         for j in ls:
#             otu = (remove(j, ["\ufeff"]))
#             OTUdict[count]["otu"] = otu
#             counter += 1
#             count += 1
#     else:
#         pass
#     count = 0
#     for j in ls:
#         otu = (remove(j, ["\ufeff"]))
#         if re.findall(r'Otu', otu):
#             OTUdict[count]["otu"] = otu
#             count += 1
#         else:
#             pass


count = 0
mothur_output = open("/Volumes/Seagate_Backup_Plus_Drive/Ongoing Research Projects/JdF2018/Algorithm/cleanOTUs.csv", "r")
otus = open("/Users/arkadiygarber/Desktop/OTU_diversity_in_Fluids/OTUs.csv")
OTUlist = []
for i in otus:
    ls = i.rstrip().split("\t")
    for j in ls:
        OTUlist.append(remove(j, ["\ufeff"]))
print(OTUlist)

BlankDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
OTUdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
counter = 0
for line in mothur_output:
    ls = (line.rstrip().split(","))
    if counter == 0:
        for i in ls:
            if re.findall(r'Otu\d', i):
                break
            else:
                counter += 1
        firstOTU = counter
    else:
        sample = ls[1]
        otu = 0
        # OTU = stabilityCounter(otu)
        for OTUabund in ls[firstOTU:]:
            # Otu = ("Otu" + str(OTU))
            BlankDict[sample][OTUlist[otu]] = OTUabund
            otu += 1
            # OTU = stabilityCounter(otu)


for i in BlankDict.keys():
    print(i)
    print(BlankDict[i])


outfile = open("/Users/arkadiygarber/Desktop/OTUs_in_NorthPond_and_JdF-summary.csv", "w")
outfile.write("OTU" + "," + "Fluids Metagenome" + "," + "Total Metagenome Reads Recruited" + "," + "1026B_D_Bas" + "," + "1301A_W_Bas" +
              "," "1362A_D_Bas2" + "," "1362A_D_Pyr" + "," "1362A_W_EnrBas" + "," "1362A_W_EnrPyr" + "," +
              "1362B_D_Bas" + "," + "1362B_D_Pyr" + "," "\n")

for i in BlastDictNP.keys():
    otu = i.split("_")[7].split("|")[0]
    if len(BlastDictNP[i]) > 9:
        outfile.write(otu + "," + "North Pond" + "," + str(len(BlastDictNP[i])) + "," + str(BlankDict["1026B_D_Bas"][otu])
                      + "," + str(BlankDict["1301A_W_Bas"][otu]) + "," + str(BlankDict["1362A_D_Bas2"][otu]) + ","
                      + str(BlankDict["1362A_D_Pyr"][otu]) + "," + str(BlankDict["1362A_W_EnrBas"][otu]) + "," +
                      str(BlankDict["1362A_W_EnrPyr"][otu]) + "," + str(BlankDict["1362B_D_Bas"][otu]) + "," +
                      str(BlankDict["1362B_D_Pyr"][otu]) + "\n")
        print(i)
        print(len(BlastDictNP[i]))
        print("")

for i in BlastDictJDF.keys():
    otu = i.split("_")[7].split("|")[0]
    if len(BlastDictJDF[i]) > 9:
        outfile.write(otu + "," + "Juan de Fuca" + "," + str(len(BlastDictJDF[i])) + "," + str(BlankDict["1026B_D_Bas"][otu])
                      + "," + str(BlankDict["1301A_W_Bas"][otu]) + "," + str(BlankDict["1362A_D_Bas2"][otu]) + ","
                      + str(BlankDict["1362A_D_Pyr"][otu]) + "," + str(BlankDict["1362A_W_EnrBas"][otu]) + "," +
                      str(BlankDict["1362A_W_EnrPyr"][otu]) + "," + str(BlankDict["1362B_D_Bas"][otu]) + "," +
                      str(BlankDict["1362B_D_Pyr"][otu]) + "\n")
        print(i)
        print(len(BlastDictJDF[i]))
        print("")
