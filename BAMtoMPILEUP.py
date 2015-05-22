#!/usr/bin/python

from subprocess import call
import sys

def usage():
    if len(sys.argv) < 2:
        print("Usage: <program> <RGfile>")
        sys.exit(0)

def realign(inputList, pat):
    """Run realign for all samples of a patient"""
#    print("Processing sample " + pat)
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -R /home/mrood/data/Mtb/ref.v3/MtbNCBIH37Rv.fa -T RealignerTargetCreator -I ' + " -I ".join(inputList) + ' -o {pat}.intervals'.format(pat=pat), shell=True)
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -R /home/mrood/data/Mtb/ref.v3/MtbNCBIH37Rv.fa -T IndelRealigner -I ' + " -I ".join(inputList) + ' -targetIntervals {pat}.intervals -nWayOut .realn.bam'.format(pat=pat), shell=True)

def callable_loci(samp):
    """Run callable loci on samples"""
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T CallableLoci -I {samp}.realn.bam -summary {samp}_defaults.summary -o {samp}_defaults.bed -R /home/mrood/data/Mtb/ref.v3/MtbNCBIH37Rv.fa'.format(samp=samp), shell=True)
    call('java -Xmx4g -jar /opt/PepPrograms/RGAPipeline/GenomeAnalysisTK.jar -T CallableLoci -I {samp}.realn.bam -summary {samp}_strict.summary -o {samp}_strict.bed -R /home/mrood/data/Mtb/ref.v3/MtbNCBIH37Rv.fa -frlmq 0.04 -mmq 20'.format(samp=samp), shell=True)

def ind_mpileup(samp):
    """Make mpileup for individual samples"""
    call('samtools mpileup -B -q 20 -s -O -l /home/mrood/data/Mtb/diversity_scales/rev2.0/process/remReg/150424_includeRegions.bed -f /home/mrood/data/Mtb/ref.v3/MtbNCBIH37Rv.fa {samp}.ready.realn.bam > {samp}_Q20.pileup'.format(samp=samp), shell=True)
#    call('samtools mpileup -B -u -v -S -q 20 -Q 20 -f /home/mrood/data/Mtb/ref.v3/MtbNCBIH37Rv.fa {samp}.realn.bam > {samp}'.format(samp=samp), shell=True)

def rem_repReg(samp):
    """Remove PPE, phage, rep, etc."""
    call('perl /home/peplab/src/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input {samp}_30.pileup --gtf /home/mrood/scripts/PPEremoval/combinedIntervals.gtf --output {samp}_30.noRep.pileup'.format(samp=samp), shell=True)

def rem_indels(samp):
    """Remove indels"""
    call('perl /home/peplab/src/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input {samp}_Q20.pileup --output {samp}_Q20.indelreg.gtf'.format(samp=samp), shell=True)
    call('perl /home/peplab/src/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input {samp}_Q20.pileup --gtf {samp}_Q20.indelreg.gtf --output {samp}_Q20.noIndel.pileup'.format(samp=samp), shell=True)

def rem_indels_mpileup(pat):
    """Remove indels"""
    #call('perl /home/peplab/src/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input {pat}.mpileup --output {pat}.indelreg.gtf'.format(pat=pat), shell=True)
    call('perl /home/peplab/src/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input {pat}.mpileup --gtf {pat}.indelreg.gtf --output {pat}.noIndel.mpileup'.format(pat=pat), shell=True)

def rem_SB(pat):
    call('perl /home/peplab/src/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input {pat}.noIndel.mpileup --gtf /home/mrood/data/Mtb/diversity_scales/rev2.0/process/remReg/SBandTBFail.gtf --output {pat}_filtered.mpileup'.format(pat=pat), shell=True)

def get_RG(inFileName):
    patDict = {}
    with open(inFileName, 'r') as inFile:
        for line in inFile:
            RGID,RGSM,RGLB,RGPL,pair,subm,model,samp,strat,patient=line.strip().split('\t')
            if patient in patDict:
                patDict[patient].append(RGID)
            else:
                patDict[patient] = [RGID]
    print("There is/are " + str(len(patDict)) + " patient in this dataset.")
    for pat in patDict:
        print("There are " + str(len(patDict[pat])) + " samples for patient " + pat + ".") 
    return patDict

def pat_mpileup(inputList, pat):
    """Make mpileup for patient"""
    print("Processing patient " + pat)
    print(" ".join(inputList))
#    call('samtools mpileup -B -f /home/mrood/data/Mtb/ref.v3/MtbNCBIH37Rv.fa -q 20 -Q 20 ' + " ".join(inputList) + ' > {pat}.mpileup'.format(pat=pat), shell=True)
    call('samtools mpileup -B -q 20 -l /home/mrood/data/Mtb/diversity_scales/rev2.0/process/remReg/150424_includeRegions.bed -f /home/mrood/data/Mtb/ref.v3/MtbNCBIH37Rv.fa ' + " ".join(inputList) + ' > {pat}.mpileup'.format(pat=pat), shell=True)
    
def pat_sync(pat):
    """Make sync for patient"""
    print("Processing patient " + pat)
    call('java -ea -Xmx7g -jar /home/peplab/src/popoolation2_1201/mpileup2sync.jar --input {pat}_filtered.mpileup --output {pat}_Q20.sync --fastq-type sanger --min-qual 20 --threads 4'.format(pat=pat), shell=True)
    call('java -ea -Xmx7g -jar /home/peplab/src/popoolation2_1201/mpileup2sync.jar --input {pat}_filtered.mpileup --output {pat}_Q30.sync --fastq-type sanger --min-qual 30 --threads 4'.format(pat=pat), shell=True)

#Convert the relevant bam files into an mpileup (by population, 
#i.e. all samples from #31 together):
#samtools mpileup -B -f ./<path to ref>/MtbNCBIH37Rv.fa -q 20 \
#-Q 20 ./<path to bams - sep by space>/ > {Prefix}.pileup 

#Convert the mpileups to sync files:
#java -ea -Xmx7g -jar \
#/home/peplab/src/popoolation2_12-1/mpileup2syn.jar \
#--input {Prefix}.mpileup --output {Prefix}.sync --fastq-type sanger \
#--min-qual 20 --threads 4

Prefix = raw_input("Prefix: ")
def snp_frequency_diff(Prefix):
    #Calculate allele frequency differences
    call('perl /opt/PepPrograms/popoolation2_1201/snp-frequency-diff.pl --input {Prefix}_Q20.sync --output {Prefix}_Q20_mc6 --min-count 6 --min-coverage 10 --max-coverage 900'.format(Prefix=Prefix), shell=True)
    call('perl /opt/PepPrograms/popoolation2_1201/snp-frequency-diff.pl --input {Prefix}_Q30.sync --output {Prefix}_Q30_mc6 --min-count 6 --min-coverage 10 --max-coverage 900'.format(Prefix=Prefix), shell=True)

def fisher_test(Prefix):
    #Estimate the significance of allele frequency differences
    call('perl /opt/PepPrograms/popoolation2_1201/fisher-test.pl --input {Prefix}_Q20.sync --output {Prefix}_Q20_mc6.fet --min-count 6 --min-coverage 10 --max-coverage 900 --min-covered-fraction 1 --window-size 1 --step-size 1 --suppress-noninformative'.format(Prefix=Prefix), shell=True)
    call('perl /opt/PepPrograms/popoolation2_1201/fisher-test.pl --input {Prefix}_Q30.sync --output {Prefix}_Q30_mc6.fet --min-count 6 --min-coverage 10 --max-coverage 900 --min-covered-fraction 1 --window-size 1 --step-size 1 --suppress-noninformative'.format(Prefix=Prefix), shell=True)

def fst_sliding(Prefix):
    #Calculate Fst values using a sliding-window approach
    call('perl /opt/PepPrograms/popoolation2_1201/fst-sliding.pl --input {Prefix}_Q20.sync --output {Prefix}_Q20_mc6_p10K.fst --min-count 4 --min-coverage 10 --max-coverage 900 --min-covered-fraction 1 --window-size 1 --step-size 1 --suppress-noninformative --pool-size 10000'.format(Prefix=Prefix), shell=True)
    call('perl /opt/PepPrograms/popoolation2_1201/fst-sliding.pl --input {Prefix}_Q30.sync --output {Prefix}_Q30_mc6_p10K.fst --min-count 4 --min-coverage 10 --max-coverage 900 --min-covered-fraction 1 --window-size 1 --step-size 1 --suppress-noninformative --pool-size 10000'.format(Prefix=Prefix), shell=True)

#for samp in samples:
#    print("Processing sample: " + samp)
#    remove_ambMap(samp)
#    callable_loci(samp)

#snp_frequency_diff(Prefix)
#fisher_test(Prefix)
#fst_sliding(Prefix)

def make_fetDict(Prefix):
    fetDict = {}
    with open('{Prefix}_mc6.fet'.format(Prefix=Prefix), 'r') as fetFile:
        for line in fetFile:
            line = line.strip().split('\t')
            pos = int(line[1])
            fetDict[pos] = {}
            paircomps = line[5:]
            for i in paircomps:
                key,fet = i.split("=")
                fetDict[pos][str(key)] = float(fet) 
    print('fetDict made')
    return fetDict

def make_alleleFreqDict(Prefix):
    afDict = {}
    with open('{Prefix}_mc6_pwc'.format(Prefix=Prefix), 'r') as afFile:
        for i, line in enumerate(afFile):
            line = line.strip().split('\t')
            if i == 0:
                keys = []
                preKeys = line[8:]
                for spot in preKeys:
                    prekey = spot.split(":")[1]
                    key = prekey.replace('-',':')
                    print key
                    keys.append(key)
            elif i > 0:
                pos = int(line[1])
                afDict[pos] = {}
                paircomps = line[8:]
                for position, item in enumerate(paircomps):
                    afDict[pos][str(keys[position])] = item
    print('afDict made')
    return afDict

def make_figureFile(Prefix, fetDict, afDict):
    with open("{Prefix}_mc6_p10K.fst".format(Prefix=Prefix), 'r') as fstFile, open("{Prefix}_mc6_summary.txt".format(Prefix=Prefix), 'w') as outFile:
        outFile.write("%s\t%s\t%s\t%s\t%s\n" %
        ("key",
        "position",
        "fst",
        "afd",
        "fet")
        )
        for line in fstFile:
            line = line.strip().split('\t')
            pos = int(line[1])
            paircomps = line[5:]
            if pos in fetDict:
                for i in paircomps:
                    key = str(i.split("=")[0])
                    fst = float(i.split("=")[1])
                    outFile.write("%s\t%i\t%f\t%s\t%f\n" %
                    (key,
                    pos,
                    fst,
                    afDict[pos][key],
                    fetDict[pos][key])
                    )
            else:
                continue
    print('file written')                


#############################
usage()
script, inFileName = sys.argv

patDict = get_RG(inFileName)
for pat in patDict:
    idList = patDict[pat]
    inputRealnList = []
    inputMpileupList = []
    for i in idList:
        ready = i + '.ready.bam'
        realn = i + '.ready.realn.bam'
        inputRealnList.append(ready)
        inputMpileupList.append(realn)
#    realign(inputRealnList, pat)
#    pat_mpileup(inputMpileupList, pat)   
#    rem_indels_mpileup(pat)
#    rem_SB(pat)
#    pat_sync(pat)
#    map(ind_mpileup, inputMpileupList)
#    snp_frequency_diff(pat)
#    fisher_test(pat)
#    fst_sliding(pat)
        
with open(inFileName, 'r') as inFile:
    for line in inFile:
        RGID,RGSM,RGLB,RGPL,pair,subm,model,samp,strat,patient=line.strip().split('\t')
#        RGID,RGSM,RGLB,RGPL,pair=line.strip().split('\t')
#        print("Processing sample " + RGID)
#        callable_loci(RGID)         
#        ind_mpileup(RGID)
#        rem_repReg(RGID)
#        rem_indels(RGID)

fetDict = make_fetDict(Prefix)
afDict = make_alleleFreqDict(Prefix)
make_figureFile(Prefix, fetDict, afDict)
