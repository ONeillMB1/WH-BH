#!/usr/bin/env python

import sys

def usage():
    if len(sys.argv) < 2:
        print("Usage: <program> <inputFile> ")
        sys.exit(0)


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

def make_sampList(inputFile):
    """From an input file extract all the samples to be processed"""
    samples = []
    with open(inputFile, 'r') as inFile:
        for i, line in enumerate(inFile):
            if i > 0:
                line = line.strip().split('\t')
                samples.append(line[0])
    print("{0} samples to be processed.".format(len(samples)))
    return samples

def make_snpList(samples):
    """Make a set of all SNPs called in the subsampled pileups for a given sample and write them to file"""
    for samp in samples:
        outFileName = samp + '.sub.snps'
        snps = set()
        for i in range(1,10):
            snpFileName = samp + '_mq20_rand{i}.bygene.snps'.format(i=i)
            #print("Adding snps from {0}".format(snpFileName))
            with open(snpFileName, 'r') as snpFile:
                for i, line in enumerate(snpFile):
                    line = line.strip().split('\t')
                    if len(line) > 1:
                        try:
                            pos = int(line[1])
                            #print pos
                            snps.add(pos)
                        except ValueError:
                            print("ValueError encountered in {0}".format(snpFileName))
                            print i
        with open(outFileName, 'w') as outFile:
            for var in snps:
                outFile.write(str(var) + '\n')

def extract_pileup(samp):
    snps = set()
    snpFileName = samp + '.sub.snps'
    print snpFileName
    with open(snpFileName, 'r') as snpFile:
        for line in snpFile:
            line = line.strip()
            snps.add(int(line))
    print("There are {0} snps for {1}".format(len(snps),snpFileName))
    pileupFileName = '/home/mrood/WH-BH/rawData/' + samp + '_Q20_filtered.pileup'
    outFileName = samp + '_sub_snps.pileup'
    print('Reading {0}'.format(pileupFileName))
    with open(pileupFileName, 'r') as pileupFile, open(outFileName, 'w') as outFile:
        for line in pileupFile:
            line = line.strip().split('\t')
            pos = int(line[1])
            if pos in snps:
                outFile.write("\t".join(line) + '\n')


def extract_baseq(samp):
    snps = set()
    snpFileName = samp + '.sub.snps'
    print snpFileName
    with open(snpFileName, 'r') as snpFile:
        for line in snpFile:
            line = line.strip()
            snps.add(int(line))
    print("There are {0} snps for {1}".format(len(snps),snpFileName))
    baseqFileName = samp + '.baseq.ext'
    outFileName = samp + 'baseq.snps'
    print('Reading {0}'.format(baseqFileName))
    with open(baseqFileName, 'r') as baseqFile, open(outFileName, 'w') as outFile:
        next(baseqFile)
        for line in baseqFile:
            line = line.strip().split('\t')
            pos = int(line[1])
            if pos in snps:
                outFile.write("\t".join(line) + '\n')

def make_catDict():
    catDict = {}
    with open("/home/mrood/scripts/importantFiles/TBGeneSet.txt", 'r') as catFile:
        for line in catFile:
            line = line.strip().split('\t')
            cat = line[0]
            catDict[cat] = [i.replace("c","") for i in line[2].split()]
    return catDict

def annotate_snp(samp, geneDict, catDict):
    inFile = samp + 'baseq.snps'
    snpDict = {}
    with open(inFile, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            snp = int(line[1])
            maf = float(line[7])/float(line[3])
            snpDict[snp] = [line, maf, [], []]
            for gene in geneDict:
                if snp >= geneDict[gene][0] and snp <= geneDict[gene][1]:
                    snpDict[snp][2].append(gene)
    outFile = samp + 'baseq.snps.gene'
    for snp in snpDict:
        if len(snpDict[snp][2]) < 1:
            snpDict[snp][3].append(["yes"] + ["no"]*32) #yes for 'all'
        else:
            for i, geneX in enumerate(snpDict[snp][2]):
                gene = geneX.replace("c", "")
                if i == 0:
                    catList = ["yes"] #yes for 'all'
                if i > 0:
                    catList = ["no"] #only count each snp once for 'all'
                for cat in sorted(catDict, key=catDict.get):
                    if gene in catDict[cat]:
                        catList.append("yes")
                    else:
                        catList.append("no")
                snpDict[snp][3].append(catList)
    outFile = samp + 'baseq.snps.ann'
    with open(outFile, 'w') as outfile:
        outfile.write('chrom\tpos\tref\treads_all\treads_pp\t\
matches\tmatches_pp\tmismatches\tmismatches_pp\t\
rms_baseq\trms_baseq_pp\trms_baseq_matches\trms_baseq_matches_pp\t\
rms_baseq_mismatches\trms_baseq_mismatches_pp\talt_allele_freq\t\
gene\tall\tTL.info\tTraSH.invitro\tTraSH.noness\tTL.CHP\t\
TL.cell\tTL.inter\tTL.reg\tTL.vir\tTL.LIP\tTraSH.invivo\t\
COG.Q\tTIM\tCOG.S\tCOG.O\tCOG.K\tCOG.F\tCOG.L\tCOG.V\tCOG.I\t\
COG.R\tCOG.J\tCOG.G\tCOG.P\tCOG.C\tCOG.M\tCOG.T\tCOG.E\tCOG.A\t\
COG.H\tCOG.D\tCOG.U\tCOG.N\n') #header
        for snp in snpDict:
            for i, gene in enumerate(snpDict[snp][2]):
                outfile.write('%s\t%f\t%s\t%s\n' %
                ("\t".join(snpDict[snp][0]),
                snpDict[snp][1],
                snpDict[snp][2][i],
                "\t".join(snpDict[snp][3][i]))
                )
    return snpDict
