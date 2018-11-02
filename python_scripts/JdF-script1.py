from collections import defaultdict
import re, os


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


# ****************************************************************************************************
# ************************** unbinned contigs identification *****************************************
# ****************************************************************************************************

allScaffs = open("/Users/arkadiygarber/Desktop/JdF1362B-anvi/scaffolds-fixed.fasta", "r")
allScaffs = fasta(allScaffs)

scaffDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
DASscaffs = open("/Users/arkadiygarber/Desktop/JdF1362B-anvi/DASbins.das", "r")
for i in DASscaffs:
    ls = i.rstrip().split("\t")
    scaffDict[ls[0]] = ls[1]


outfile = open("/Users/arkadiygarber/Desktop/JdF1362B-anvi/unbinned_scaffolds.fasta", "w")
count = 0
for i in allScaffs.keys():
    if i not in scaffDict.keys():
        print(i)
        print(len(allScaffs[i]))
        print("")
        outfile.write(">" + i + "\n")
        outfile.write(allScaffs[i] + "\n")
        count += 1
print(count)


# ****************************************************************************************************
# *********************** checkm-pulled SSU genes blast against amplicons ****************************
# ****************************************************************************************************


outcsv = open("/Users/arkadiygarber/Desktop/JdF1362B-anvi/SSU_seqs/Reconstructed_SSU/OTUs_to_bins_to_Full16S.csv", "w")
outcsv.write("bin" + "," + "otu" + "," + "identity" + "," + "length of alignment" + "," + "reconstructed 16S length" + "\n")
outfasta = open("/Users/arkadiygarber/Desktop/JdF1362B-anvi/SSU_seqs/Reconstructed_SSU/OTUs_to_bins_to_Full16S.fasta", "w")


seqs = open("/Users/arkadiygarber/Desktop/JdF1362B-anvi/SSU_seqs/ssu.fna", "r")
seqs = fasta(seqs)

binList = []
dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
OTUs = open("/Users/arkadiygarber/Desktop/JdF1362B-anvi/SSU_seqs/JdF1362B-CORKs-ssu.blast", "r")
for i in OTUs:
    ls = i.rstrip().split("\t")
    e = ls[10]
    align = ls[3]
    id = ls[2]
    otu = ls[0]
    bin = ls[1]
    if float(id) >= 99.0 and int(align) >= 300 and bin not in binList:
        outfasta.write(">" + bin + "\n")
        outfasta.write(seqs[bin] + "\n")
        outcsv.write(bin + "," + otu + "," + str(id) + "," + str(align) + "," + str(len(seqs[bin])) + "\n")
        binList.append(bin)
        print(otu)
        print(bin)
        print(e)
        print(align)
        print(id)
        print(seqs[bin])
        print(len(seqs[bin]))
        print("")

binsA = os.listdir("/Volumes/Seagate Backup Plus Drive/Juan de Fuca/JdF1362A-anvi/manual_bins")
binsB = os.listdir("/Volumes/Seagate Backup Plus Drive/Juan de Fuca/JdF1362B-anvi/manual_bins_contigs")


# ****************************************************************************************************
# ************************************ creating files files ******************************************
# ****************************************************************************************************


blast = open("/Volumes/Seagate Backup Plus Drive/Juan de Fuca/JdF1362AB-anvi/JDF16S.blast", "r")
ssu = open("/Volumes/Seagate Backup Plus Drive/Juan de Fuca/JdF1362AB-anvi/ssu.fna", "r")
ssu = fasta(ssu)
dict = defaultdict(lambda: defaultdict(list))
for i in blast:
    i = i.rstrip()
    ls = i.split("\t")
    dict[ls[1]]["otu"].append(ls[0])
    dict[ls[1]]["id"].append(ls[2])
    dict[ls[1]]["aln"].append(ls[3])

out = open("/Users/arkadiygarber/Desktop/final_matches.csv", "w")
out.write("Bin&&contig" + "," + "otu" + "," + "id" + "," + "aln" + "," + "16S_len" + "," + "Seq" + "\n")
for i in dict.keys():
    print(i)
    for j in range(0, len(dict[i]["otu"])):
        out.write(i + "," + dict[i]["otu"][j] + "," + dict[i]["id"][j] + "," + dict[i]["aln"][j] + "," +
                  str(len(ssu[i])) + "," + str(ssu[i]) + "\n")
    print(dict[i]["otu"])
    print(dict[i]["id"])
    print(dict[i]["aln"])
    print("")