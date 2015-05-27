#!/usr/bin/env python

import sys

def usage():
    if len(sys.argv) < 2:
        print("Usage: <program> <sampleFile>")
        sys.exit(0)

def make_sampleInfo(sampleFile):
    prefixes = []
    sampDict = {}
    with open(sampleFile, 'r') as infile:
        for line in infile:
            line=line.strip().split('\t')
            #pat,samp,prefix,study,studySamp = line.split()
            prefix = line[0]
            samp = line[7]
            pat = line[9]
            print(pat + '\t' + samp + '\t' + prefix)
            sampDict[prefix] = [pat, samp]
            prefixes.append(prefix)
    return prefixes, sampDict

def make_dataDicts(prefixes):
    d = {}
    stats = ["pi", "theta", "TD"]
    for stat in stats:
        d[stat] = {}
        for prefix in prefixes:
            d[stat][prefix] = {}
    print("There are " + str(len(d)) + " stats in dictionary")
    if len(d["pi"]) == len(d["theta"] and len(d["pi"]) == len(d["TD"]):
        print("There are " + str(len(d["pi"])) + " samples for each stat.")
    return d, stats

def make_pi(prefixes, d):
    for prefix in prefixes:
        for i in range(1,10):
            piFileName = 'pi/' + prefix + '_mq20_rand{i}_filtered_c50_w100K_s10K_n10K.pi'.format(i=i)
            print piFileName
            with open(piFileName, 'r') as piFile:
                for line in piFile:
                    line = line.strip()
                    chrom,window,snps,cov,value = line.split()
                    if window in d["pi"][prefix]:
                        d["pi"][prefix][window].append(value)
                    else:
                        d["pi"][prefix][window] = [value]

def make_theta(prefixes, d):
    for prefix in prefixes:
        for i in range(1,10):
            thetaFileName = 'theta/' + prefix + '_mq20_rand{i}_filtered_c50_w100K_s10K_n10K.theta'.format(i=i)
            with open(thetaFileName, 'r') as thetaFile:
                for line in thetaFile:
                    line = line.strip()
                    chrom,window,snps,cov,value = line.split()
                    if window in d["theta"][prefix]:
                        d["theta"][prefix][window].append(value)
                    else:
                        d["theta"][prefix][window] = [value]
                
def make_TD(prefixes, d):
    for prefix in prefixes:
        for i in range(1,10):
            TDFileName = 'TD/' + prefix + '_mq20_rand{i}_filtered_c50_w100K_s10K_n10K.tajD'.format(i=i)
            with open(TDFileName, 'r') as TDFile:
                for line in TDFile:
                    line = line.strip()
                    chrom,window,snps,cov,value = line.split()
                    if window in d["TD"][prefix]:
                        d["TD"][prefix][window].append(value)
                    else:
                        d["TD"][prefix][window] = [value]

def write_files(prefixes, sampDict, d, stats):
    for stat in stats:
        with open(stat + '_mq20.out2', 'w') as outfile:
            #outfile.write(
            #'prefix\twindow\tvalue\n'
            #)
            for prefix in prefixes:
                for window in d[stat][prefix]:
                    #print type(sampDict[prefix][0])
                    #print type(sampDict[prefix][1])
                    #print type(window)
                    #print type(d[stat][prefix][window])
                    outfile.write("%s\t%s\t%s\t%s\t%s\n" %
                    (sampDict[prefix][0],
                    sampDict[prefix][1],
                    prefix,
                    window,
                    "\t".join(d[stat][prefix][window])))
                    
########################
usage()
script, sampleFile = sys.argv

prefixes, sampDict = make_sampleInfo(sampleFile)
d, stats = make_dataDicts(prefixes)
make_pi(prefixes, d)
make_theta(prefixes, d)
make_TD(prefixes, d)
write_files(prefixes, sampDict, d, stats)
