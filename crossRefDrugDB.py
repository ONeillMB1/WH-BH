#!/usr/bin/env python

import sys

def usage():
    if len(sys.argv) < 5:
        print("Usage: <program> <inputFile> <annFile> <syncFile> <outFile>")
        sys.exit(0)

def make_geneDict():
    geneDict = {}
    with open("/home/mrood/scripts/importantFiles/PerGeneSummary.txt", 'r') as dat:
        for line in dat:
            line = line.strip().split('\t')
            start = int(line[4])
            stop = int(line[5])
            gene = line[0]
            geneDict[gene] = [start,stop,line[6],line[1],line[7]]
    return geneDict

def cleanup_mutList():
    geneMutDict = {}
    with open("/home/mrood/WH-BH/global_data/DownloadDB.txt", 'r') as infile:
        for i, line in enumerate(infile):
            line = line.strip().split("\t")
            if i > 0:
                if line[10] != "" and line[10] != "?" and line[10] != "See Note":
                    try:
                        mut = int(line[10])
                    except ValueError:
                        mut = line[10] 
                    geneMutDict[int(line[0])] = (line[1], mut)
                else:
                    print line[0]
                    continue
    return geneMutDict

def make_mutList(geneMutDict, geneDict):
    geneDict["MTB000019"] = [1471846, 1473382, "+", "rrs", "Ribosomal RNA 16S"]
    geneDict["Rv2427A"] = [2725571, 2726087, "-", "oxyR", "hypothetical protein"]
    mutations = []
    for key in geneMutDict:
        gene = geneMutDict[key][0]
        start,stop,strand = geneDict[gene][0:3]
        if isinstance(geneMutDict[key][1], int):
            nt = geneMutDict[key][1]
            try:
                if strand == "-":
                    if nt > 0:
                        pos = stop - nt + 1
                    else:
                        print("mut on neg strand and upstream")
                        pos = stop - nt
                else:
                    if nt > 0:
                        pos = start + nt - 1
                    else:
                        pos = start + nt
                        print("mut on pos strand and upstream")
                mutations.append(pos) 
            except KeyError:
                print("{0} not in gene dictionary".format(gene))
        else:
            print geneMutDict[key][1]
            muts = geneMutDict[key][1].split("-")
            print muts
            if len(muts) == 2:
                nts = range(int(muts[0]),int(muts[1]) + 1)
                print nts
            else: 
                muts2 = geneMutDict[key][1].strip('"').split(",")
                print muts2
                try:
                    nts = [ int(x) for x in muts2 ]
                    print nts
                except ValueError:
                    nts = []
            for i in nts:
                if strand == "-":
                    if i > 0:
                        pos = stop - i + 1
                    else:
                        print("mut on neg strand and upstream")
                        pos = stop - i
                else:
                    if i > 0:
                        pos = start + i - 1
                    else:
                        pos = start + i
                        print("mut on pos strand and upstream")
                mutations.append(i)
    return mutations            

def write_overlap(mutations):
    mutSet = set(mutations)
    print len(mutSet)
    print len(mutations)
    with open("/home/mrood/WH-BH/global_data/150520_global_backup/masked_global_pass_syn-nonsyn.snps", 'r') as infile, open("/home/mrood/WH-BH/global_data/150520_global_backup/crossRefDrugResMutGlobal.txt", 'w') as outfile:
        for line in infile:
            line = line.strip().split('\t')
            pos = int(line[2])
            if pos in mutSet:
                outfile.write("\t".join(line) + '\n')

