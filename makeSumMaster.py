#!/usr/bin/env python

sampleFile = raw_input("sample file path: ")
parameters = raw_input("parameters, i.e. mq20, mq20_mc3, mq30: ")
outfileName = raw_input("outfile name: ")

def make_geneList():
    geneList = []
    with open("/home/mrood/scripts/importantFiles/genes_ordered.txt", 'r') as infile:
        for line in infile:
            gene = line.strip()
            geneList.append(gene)
    return geneList

def make_geneDict():
    geneDict = {}
    with open("/home/mrood/scripts/importantFiles/PerGeneSummary.txt", 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            geneDict[line[0]] = [line[1], line[3], line[7]]
    #print geneDict["Rv0001"] #test to make sure generated 
    return geneDict

def make_sampleList(sampleFile):
    samples = []
    with open(sampleFile, 'r') as infile:
        for line in infile:
            line = line.strip()
            samples.append(line)
    print('There are {0} samples to be processed.'.format(len(samples)))
    return samples

def make_dicts(geneList, samples):
    d = {}
    stats = ["theta", "pi", "TD", "piNpiS"]
    for sample in samples:
        d[sample] = {}
        for stat in stats:
            d[sample][stat] = {}
        for gene in geneList:
            d[sample]["theta"][gene] = ["-", "-"]
            d[sample]["pi"][gene] = "-"
            d[sample]["TD"][gene] = "-"
            d[sample]["piNpiS"][gene] = ["-","-","-"]
    d["global"] = {}
    for stat in stats:
        d["global"][stat] = {}
    for gene in geneList:
        d["global"]["theta"][gene] = ["-", "-", "-"]
        d["global"]["pi"][gene] = "-"
        d["global"]["TD"][gene] = "-"
        d["global"]["piNpiS"][gene] = ["-","-","-","-","-","-","-"]
    return d

def pop_dg(d):
    with open("global/global.genes.tajD", 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            gene = line[0]
            try:
                cov = float(line[2])
            except ValueError:
                cov = line[2]
            if type(cov) == float and cov >= 0.5:
                d["global"]["TD"][gene] = line[3]
    with open("global/global.genes.pi", 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            gene = line[0]
            try:
                cov = float(line[2])
            except ValueError:
                cov = line[2]
            if type(cov) == float and cov >= 0.5:
                d["global"]["pi"][gene] = line[3]
    passGenes = []
    with open("global/global.genes.theta", 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            gene = line[0]
            try:
                cov = float(line[2])
            except ValueError:
                cov = line[2]
            if type(cov) == float and cov >= 0.5:
                d["global"]["theta"][gene] = line[1:]
                passGenes.append(gene)
    print("global sample has " + str(len(passGenes)) + " genes with >= 0.5 coverage")
    with open("global/global_syn-nonsyn.pi", 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            gene = line[0]
            if gene in passGenes:
                if float(line[5]) == 0:
                    if float(line[6]) == 0:
                        NS = "NaN"
                    else:
                        NS = "DQ"
                elif float(line[6]) == 0:
                    NS = "Inf"
                else:
                    NS = float(line[5])/float(line[6])
                d["global"]["piNpiS"][gene] = [line[1],line[2],line[3],line[4],line[5],line[6],str(NS)]
    return d

def pop_di(d, samples):
    for sample in samples:
        with open('tajD/' + sample + '_' + parameters + '_tajD.genes.mean', 'r') as file:
            next(file)
            for line in file:
                line = line.strip().split('\t')
                gene = line[0]
                try:
                    cov = float(line[1])
                except ValueError:
                    cov = line[1]
                if type(cov) == float and cov >= 0.5:
                    d[sample]["TD"][gene] = line[12]
        with open('pi/' + sample + '_' + parameters + '_pi.genes.mean', 'r') as infile:
            next(infile)
            for line in infile:
                line = line.strip().split('\t')
                gene = line[0]
                try:
                    cov = float(line[1])
                except ValueError:
                    cov = line[1]
                if type(cov) == float and cov >= 0.5:
                    d[sample]["pi"][gene] = line[12]
        passGenes = []
        with open("theta/" + sample + '_' + parameters + '_theta.genes.mean', 'r') as infile3:
            for line in infile3:
                line = line.strip().split('\t')
                gene = line[0]
                try:
                    cov = float(line[1])
                except ValueError:
                    cov = line[1]
                if type(cov) == float and cov >= 0.5:
                    d[sample]["theta"][gene] = [line[1],line[12]]
                    passGenes.append(gene)
        print(sample + " has " + str(len(passGenes)) + " genes with >= 0.5 coverage")
        #with open('piNpiS/' + sample + '_med.calc', 'r') as infile4:
        #    next(infile4)
        #    for line in infile4:
        #        line = line.strip().split('\t')
        #        gene = line[0]
        #        if gene in passGenes:
        #            d[sample]["piNpiS"][gene] = [line[1],line[3],line[5]]
        with open('piNpiS/' + sample + '_' + parameters + '_med.calc', 'r') as infile4:
            next(infile4)
            for line in infile4:
                line = line.strip().split('\t')
                gene = line[0]
                if gene in passGenes and line[1] != "NA":
                    if float(line[1]) == 0:
                        if float(line[2]) == 0:
                            NS = "NaN"
                        else:
                            NS = "DQ"
                    elif float(line[2]) == 0:
                        NS = "Inf"
                    else:
                        NS = float(line[1])/float(line[2])
                    d[sample]["piNpiS"][gene] = [line[1],line[2],str(NS)]
    return d


def write_files(geneList, geneDict, d, samples):
    with open(outfileName, 'w') as outfile:
        outfile.write(
        'gene.symbol\tgene.name\tgene.length\tgene.des\tg.snps\tg.coverage\tg.theta\tg.pi\tg.TD\tg.NSlen\tg.Slen\tg.NSsnps\tg.Ssnps\tg.piN\tg.piS\tg.piN.piS\t')
        for sample in samples:
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %
            (sample+'.cov',
            sample+'.theta',
            sample+'.pi',
            sample+'.tajD',
            sample+'.piN',
            sample+'.piS',
            sample+'.piNpiS'))
        
        #a1.coverage\ta1.theta\ta1.pi\ta1.TD\ta1.piN\ta1.piS\ta1.piN.piS\ta2.coverage\ta2.theta\ta2.pi\ta2.TD\ta2.piN\ta2.piS\ta2.piN.piS\ta3.coverage\ta3.theta\ta3.pi\ta3.TD\ta3.piN\ta3.piS\ta3.piN.piS\tb2.coverage\tb2.theta\tb2.pi\tb2.tajD\tb2.piN\tb2.piS\tb2.piN.piS\tc1.coverage\tc1.theta\tc1.pi\tc1.tajD\tc1.piN\tc1.piS\tc1.piN.piS\tc2.coverage\tc2.theta\tc2.pi\tc2.tajD\tc2.piN\tc2.piS\tc2.piN.piS\n')
        for gene in geneList:
            try:
                outfile.write(gene + '\t' + "\t".join(geneDict[gene]) + '\t')
            except KeyError:
                print('KeyError: ' + gene)
                outfile.write(gene + '\t\t\t')
            for sample in samples:
                outfile.write("%s\t%s\t%s\t%s\t" %
                ("\t".join(d[sample]["theta"][gene]),
                d[sample]["pi"][gene],
                d[sample]["TD"][gene],
                "\t".join(d[sample]["piNpiS"][gene])
                ))
            outfile.write('\n')

geneList = make_geneList()
geneDict = make_geneDict()
samples = make_sampleList(sampleFile)
d = make_dicts(geneList, samples)
dg = pop_dg(d)
dgi = pop_di(dg, samples)
write_files(geneList, geneDict, dgi, samples)
