#!/usr/bin/env python

import sys
import numpy as np

def make_geneDict():
    geneDict = {}
    with open("/home/mrood/scripts/importantFiles/PerGeneSummary.txt", 'r') as dat:
        for line in dat:
            line = line.strip().split('\t')
            start = int(line[4])
            stop = int(line[5])
            gene = line[0]
            geneDict[gene] = [start,stop]
    return geneDict

def make_dicts():
    mqDict = {}
    dpDict = {}
    qualDict = {}
    with open("MQmean.sort", 'r') as infile:
        for i, line in enumerate(infile):
            if i > 1:
                line = line.strip().split()
                pos, metric = int(line[0]), float(line[1])
                if pos in mqDict:
                    mqDict[pos].append(metric)
                else:
                    mqDict[pos] = [metric]
    with open("DPmean.sort", 'r') as infile:
        for i, line in enumerate(infile):
            if i > 1:
                line = line.strip().split()
                pos, metric = int(line[0]), float(line[1])
                if pos in dpDict:
                    dpDict[pos].append(metric)
                else:
                    dpDict[pos] = [metric]
    with open("QUALmean.sort", 'r') as infile:
        for i, line in enumerate(infile):
            if i > 1:
                line = line.strip().split()
                pos, metric = int(line[0]), float(line[1])
                if pos in qualDict:
                    qualDict[pos].append(metric)
                else:
                    qualDict[pos] = [metric]
    return mqDict, dpDict, qualDict

def make_masterDict():
    masterDict = {}
    for metric in ["mq", "dp", "qual"]:
        masterDict[metric] = {}
    return masterDict

def add_mq(masterDict, geneDict, mqDict):
    for gene in geneDict:
        start = int(geneDict[gene][0])
        stop = int(geneDict[gene][1])
        mqList = []
        mqList2 = []
        for i in mqDict:
            if i >= start and i <= stop:
                mqList.append(mqDict[i][0])
		try:
                    mqList2.append(mqDict[i][1])
                except IndexError:
                    continue
            else:
                continue
        if len(mqList) > 0:
            mqA = np.array(mqList)
            mq = mqA.mean()
        else:
            mq = "NA"
        if len(mqList2) > 0:
            mq2A = np.array(mqList2)
            mq2 = mq2A.mean()
        else:
            mq2 = "NA"
        masterDict["mq"][gene] = [mq,mq2]
    return masterDict

def add_dp(masterDict, geneDict, dpDict):
    for gene in geneDict:
        start = int(geneDict[gene][0])
        stop = int(geneDict[gene][1])
        dpList = []
        dpList2 = []
        for i in dpDict:
            if i >= start and i <= stop:
                dpList.append(dpDict[i][0])
		try:
                    dpList2.append(dpDict[i][1])
                except IndexError:
                    continue
            else:
                continue
        if len(dpList) > 0:
            dpA = np.array(dpList)
            dp = dpA.mean()
        else:
            dp = "NA"
        if len(dpList2) > 0:
            dp2A = np.array(dpList2)
            dp2 = dp2A.mean()
        else:
            dp2 = "NA"
        masterDict["dp"][gene] = [dp,dp2]
    return masterDict

def add_qual(masterDict, geneDict, qualDict):
    for gene in geneDict:
        start = int(geneDict[gene][0])
        stop = int(geneDict[gene][1])
        qualList = []
        qualList2 = []
        for i in qualDict:
            if i >= start and i <= stop:
                qualList.append(qualDict[i][0])
		try:
                    qualList2.append(qualDict[i][1])
                except IndexError:
                    continue
            else:
                continue
        if len(qualList) > 0:
            qualA = np.array(qualList)
            qual = qualA.mean()
        else:
            qual = "NA"
        if len(qualList2) > 0:
            qual2A = np.array(qualList2)
            qual2 = qual2A.mean()
        else:
            qual2 = "NA"
        masterDict["qual"][gene] = [qual,qual2]
    return masterDict

def write_file(masterDict):
    with open("GeneSummaryStat.txt", 'w') as outfile:
        outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
        ('gene',
        'mq', 'mq2',
        'qual', 'qual2',
        'dp', 'dp2')
        )
        for gene in sorted(masterDict['mq'], key=masterDict['mq'].get):
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' %
            (gene,
            str(masterDict['mq'][gene][0]),
            str(masterDict['mq'][gene][1]),
            str(masterDict['qual'][gene][0]),
            str(masterDict['qual'][gene][1]),
            str(masterDict['dp'][gene][0]),
            str(masterDict['dp'][gene][1])
            ))

def make_TD():
    tdDict = {}
    with open('/home/mrood/WH-BH/data/global/global.genes.tajD', 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            gene = line[0]
            td = line[3]
            tdDict[gene] = td
    return tdDict

def add_TD(tdDict):
    with open("GeneSummaryStat.txt", 'r') as infile, open("GeneTDSummaryStat.txt", 'w') as outfile:
        for i, line in enumerate(infile):
            line = line.strip().split('\t')
            if i == 0:
                outfile.write("\t".join(line) + '\t' + "td" + '\n')
            else:
                outfile.write("\t".join(line) + '\t' + tdDict[line[0]] + '\n')


#geneDict = make_geneDict()
#mqDict, dpDict, qualDict = make_dicts()
#masterDict = calc_gene(geneDict, mqDict, dpDict, qualDict)

tdDict = make_TD()
add_TD(tdDict)
