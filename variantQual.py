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
