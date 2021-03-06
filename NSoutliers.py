#!/usr/bin/env python

import sys

def usage():
    if len(sys.argv) < 5:
        print("Usage: <program> <inputFile> <annFile> <syncFile> <outFile>")
        sys.exit(0)

def make_NSdict():
    NSd = {}
    with open("/home/mrood/WH-BH/data/enrichmentGuide2.txt", 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split('\t')
            NSd[line[1]] = int(line[5])
    print("There are {0} samples.".format(len(NSd)))
    return NSd

def gen_outliers(inputFile, NSd):
    outliers = []
    patA = set()
    patB = set()
    patC = set()
    patD = set()
    patE = set()
    for samp in NSd:
        column = int(NSd[samp])
        with open(inputFile, 'r') as infile:
            next(infile)
            for line in infile:
                line = line.strip().split('\t')
                value = line[column]
                if value != 'DQ' and value != 'Inf' and value != 'NaN' and value != '-':
                    NS = float(line[column])
                    if NS > 1.0:
                        outliers.append((samp, line[0], value))
                        if samp[0] == 'A':
                            patA.add(line[0])
                        elif samp[0] == 'B':
                            patB.add(line[0])
                        elif samp[0] == 'C':
                            patC.add(line[0])
                        elif samp[0] == 'D':
                            patD.add(line[0])
                        elif samp[0] == 'E':
                            patE.add(line[0])
    print("There are {0} outliers.".format(len(outliers)))
    return outliers, patA, patB, patC, patD, patE





def summary_dict(inputFile):
    d = {}
    with open(inputFile, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split('\t')
            if int(line[1]) in d:
                d[int(line[1])].append((line[0],float(line[2]),float(line[3]),float(line[4])))
            else:
                d[int(line[1])] = [(line[0],float(line[2]),float(line[3]),float(line[4]))]
    print('There are {0} snps in the dictionary'.format(len(d)))
    return d

def simp_dict(d):
    d2 = {}
    for pos in d:
        fstVals = []
        for tup in d[pos]:
            fst = tup[1]
            fstVals.append(fst)
        #print('for {pos} there are {l} entries'.format(pos=pos, l=len(fstVals)))
        rInd = fstVals.index(max(fstVals))
        d2[pos] = [d[pos][rInd]]
    return d2

def make_geneDict():
    geneDict = {}
    with open("/home/mrood/scripts/importantFiles/PerGeneSummary.txt", 'r') as dat:
        for line in dat:
            line = line.strip().split('\t')
            start = int(line[4])
            stop = int(line[5])
            gene = line[0]
            geneDict[gene] = [start,stop,line[1],line[7]]
    return geneDict

def annotate_snp(d2, geneDict):
    for snp in d2:
        geneList = []
        commonList = []
        desList = []
        for gene in geneDict:
            if snp >= geneDict[gene][0] and snp <= geneDict[gene][1]:
                geneList.append(gene)
                commonList.append(geneDict[gene][2])
                desList.append(geneDict[gene][3])
        d2[snp].append(geneList)
        d2[snp].append(commonList)
        d2[snp].append(desList)
    return d2

def find_freq(syncFile, d2):
    freqDict = {}
    with open(syncFile, 'r') as sync:
        for line in sync:
            line = line.strip().split('\t')
            pos = int(line[1])
            ref = line[2]
            freqs = []
            samps = line[3:]
            if pos in d2:
                cS1 = samps[0].split(":")
                cS1 = [int(i) for i in cS1]
                maa = cS1.index(max(cS1))
                for samp in samps:
                    counts = samp.split(":")
                    counts = [int(x) for x in counts]
                    cov = sum(counts)
                    if cov < 10:
                        af = "NA"
                    else:
                        #if ref == "A":
                        if maa == 0:
                            ma = "A"
                            af = 1 - float(counts[0])/cov
                        #elif ref == "T":
                        elif maa == 1:
                            ma = "T"
                            af = 1 - float(counts[1])/cov
                        #elif ref == "C":
                        elif maa == 2:
                            ma = "C"
                            af = 1 - float(counts[2])/cov
                        #elif ref == "G":
                        elif maa == 3:
                            ma = "G"
                            af = 1 - float(counts[3])/cov
                    af2 = "{0:.2f}".format(af)
                    if float(af2) == 1.00:
                        af2 = "1.0"
                    freqs.append(af2)
                freqDict[pos] = freqs
                d2[pos].append(ma)
                d2[pos].append(freqs)
    return d2

def make_snpDict(annFile, d2):
    snpDict = {}
    with open(annFile, 'r') as ann:
        for line in ann:
            line = line.strip().split()
            snpDict[int(line[0])] = line[1:]
    for snp in d2:
        try:
            d2[snp].append(snpDict[snp])
        except KeyError:
            d2[snp].append(["NA","NA","NA","NA","NA","NA"])
    return d2

def make_catDict(d2):
    catDict = {}
    with open("/home/mrood/scripts/importantFiles/NumberedGeneSet.txt", 'r') as catFile:
        for line in catFile:
            line = line.strip().split('\t')
            cat = line[0]
            catDict[cat] = [i.replace("c","") for i in line[2].split()]
        for snp in d2:
            if len(d2[snp][1]) > 0:
                gene = d2[snp][1][0].replace("c","")
                catList = []
                for cat in catDict:
                    if gene in catDict[cat]:
                        catList.append(cat)
                catList = [int(i) for i in catList]
            d2[snp].append(sorted(catList))
    return d2

def make_catDict2():
    catDict = {}
    with open("/home/mrood/scripts/importantFiles/NumberedGeneSet.txt", 'r') as catFile:
        for line in catFile:
            line = line.strip().split('\t')
            cat = line[0]
            catDict[cat] = [i.replace("c","") for i in line[2].split()]
    return catDict
  
def write_outfile2(inFileName,outFileName,catDict,geneDict):
    with open(inFileName, 'r') as inFile, open(outFileName, 'w') as outFile:
        for line in inFile:
            line = line.strip()
            gene = line.replace("c", "")
            geneCatList = []
            for cat in catDict:
                if gene in catDict[cat]:
                    geneCatList.append(cat)
            outFile.write('%s\t%s\t%s\n' %
            (line,
            "\t".join(geneDict[line][2:]),
            ",".join(geneCatList)))



def write_outfile(d2,outFile):
    with open(outFile, 'w') as outfile:
        for snp in d2:
            outfile.write('%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
            ("\t".join([str(i) for i in d2[snp][0]]),
            snp,
            ", ".join(d2[snp][1]),
            ", ".join(d2[snp][2]),
            ", ".join(d2[snp][3]),
            d2[snp][4],
            ", ".join([str(i) for i in d2[snp][5]]),
            "\t".join(d2[snp][6]),
            ", ".join([str(i) for i in d2[snp][7]]))
            )

#usage()
#inputFile, annFile, syncFile, outFile = sys.argv[1:]
#d = summary_dict(inputFile)
#d2 = simp_dict(d)
#geneDict = make_geneDict()
#d2 = annotate_snp(d2,geneDict)
#d2 = find_freq(syncFile,d2)
#d2 = make_snpDict(annFile,d2)
#d2 = make_catDict(d2)
#write_outfile(d2,outFile)
