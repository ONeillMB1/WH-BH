#!/usr/bin/env python

def make_patDict():
    pD = {}
    with open("/home/mrood/WH-BH/data/enrichmentGuide.txt", 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split()
            pat = line[1][0]
            if pat in pD:
                pD[pat].append(line[0])
            else:
                pD[pat] = [line[0]]
    pD.pop("B", None)
    print('There are {0} patients to be processed.'.format(len(pD)))
    return pD

def make_snpDict(pD):
    for pat in pD:
        outfileName = pat + '_SNPannotation.txt'
        sD = {}
        for samp in pD[pat]:
            #for i in range(1,10):
            filename = samp+'_Q20_filtered_syn-nonsyn.snps'
            with open(filename) as data:
                for line in data:
                    line = line.strip().split('\t')
                    gene = line[0]
                    pos = int(line[2])
                    ref = line[3]
                    snpType = line[9]
                    strand = line[10]
                    codon = line[11]
                    aa = line[12]
                    if pos in sD:
                        continue
                    else:
                        sD[pos] = [gene,ref,snpType,strand,codon,aa]
        with open(outfileName, 'w') as outfile:
            for key in sorted(sD, key=sD.get):
                outfile.write(str(key) + '\t' + '\t'.join(sD[key]) + '\n')


pD = make_patDict()
print(pD)
make_snpDict(pD)

