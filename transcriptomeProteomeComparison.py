###
### this script computes the distribution of correlation coefficients between  transcript and protein expression
###

import sys,numpy
import scipy,scipy.stats
import matplotlib,matplotlib.pyplot

def annotationReader():

    '''
    this function reads the annotation file and creates a dictionary between systematic gene names and protein names
    '''

    annotation={}
    with open(annotationFile,'r') as f:
        header=f.readline()
        for line in f:
            vector=line.split('\t')
            proteinName=vector[1].replace('"','')
            systematicName=vector[-1].replace('"','').replace('\n','')
            annotation[proteinName]=systematicName

    return annotation

def getProteinNames(fileName):

    '''
    this function reads the protein expression file and returns the protein names
    '''

    proteinNames=[]
    with open(fileName,'r') as f:
        next(f)
        for line in f:
            organism=line.split('\t')[0]
            proteinName=line.split('\t')[1].split(';')[0]
            if organism == 'DESVH':
                proteinNames.append(proteinName)

    return proteinNames

def log2FoldTransformer(series):

    '''
    this function converts an absolute time series data into a log2 fold-change time series
    '''

    transformedSeries=[numpy.log2(element/series[0]) for element in series]
    print(series)
    print(transformedSeries)

    return transformedSeries

def proteomeFileReader(inputFileName,uniqueProteins,proteome):

    '''
    this function reads proteome expression file
    '''

    with open(inputFileName,'r') as f:
        header=f.readline()
        headerVector=header.split('\t')
        for line in f:
            vector=line.split('\t')
            for i in range(len(headerVector)):
                if genotypes[0] in headerVector[i] or genotypes[1] in headerVector[i]:
                    if '-AVG' not in headerVector[i] and '/' not in headerVector[i]:

                        # defining data
                        identifiers=headerVector[i].split('-')
                        workingGenotype=identifiers[0]
                        workingGenotype=workingGenotype.replace('D','')
                        workingTime=identifiers[1]
                        workingReplicate=identifiers[2]
                        proteinName=vector[1].split(';')[0]
                        value=float(vector[i])

                        # associating data to a variable
                        if proteinName in uniqueProteins:
                            geneName=annotation[proteinName]
                            if workingGenotype not in proteome.keys():
                                proteome[workingGenotype]={}
                            if workingReplicate not in proteome[workingGenotype].keys():
                                proteome[workingGenotype][workingReplicate]={}
                            if workingTime not in proteome[workingGenotype][workingReplicate].keys():
                                proteome[workingGenotype][workingReplicate][workingTime]={}
                            if geneName not in proteome[workingGenotype][workingReplicate][workingTime].keys():
                                proteome[workingGenotype][workingReplicate][workingTime][geneName]=value

    return proteome

def proteomeReader():

    '''
    this function reads protein expression data and returns a dictionary
    '''

    # f.1. checking that proteins were detected in both ST and SR
    STproteins=getProteinNames(stDataFile)
    print('%s proteins detected in ST conditions...'%len(STproteins))
    
    SRproteins=getProteinNames(srDataFile)
    print('%s proteins detected in SR conditions...'%len(STproteins))

    commonProteins=[]
    for element in STproteins:
        if element in SRproteins:
            commonProteins.append(element)
    uniqueProteins=list(set(commonProteins))
    print('%s proteins found in both conditions...'%len(uniqueProteins))

    # f.2. building the expression variable
    proteome={}
    proteome=proteomeFileReader(stDataFile,uniqueProteins,proteome)
    proteome=proteomeFileReader(srDataFile,uniqueProteins,proteome)

    # f.3. converting protein names to gene names
    listOfProteomicsNames=[annotation[element] for element in uniqueProteins]
                                
    return proteome,listOfProteomicsNames

def transcriptomeReader():

    '''
    this function reads transcriptome expression data and returns it in a dictionary format
    '''

    transcriptome={}
    listOfTranscriptomicsNames=[]
            
    with open(transcriptomeDataFile,'r') as f:
        header=f.readline()
        identifiers=header.split()[1:]

        for line in f:
            vector=line.split('\t')
            for i in range(len(identifiers)):
        
                # defining the data
                workingGenotype=identifiers[i].split('_')[-1].split('-')[0]
                workingLabel=identifiers[i].split('_')[-2]
                if workingLabel == 'L1':
                    workingTime='ST1'
                elif workingLabel == 'LS1':
                    workingTime='SR1'
                elif workingLabel == 'L3':
                    workingTime='ST3'
                elif workingLabel == 'LS3':
                    workingTime='SR3'
                else:
                    workingTime=None
                workingReplicate=identifiers[i].split('-')[-1]
                geneName=vector[0].replace('_','')
                value=float(vector[i+1])

                # associating data to a variable
                if 'DVU' in geneName:
                    if geneName not in listOfTranscriptomicsNames:
                        listOfTranscriptomicsNames.append(geneName)
                    if workingTime is not None:
                        if workingGenotype not in transcriptome.keys():
                            transcriptome[workingGenotype]={}
                        if workingReplicate not in transcriptome[workingGenotype].keys():
                            transcriptome[workingGenotype][workingReplicate]={}
                        if workingTime not in transcriptome[workingGenotype][workingReplicate].keys():
                            transcriptome[workingGenotype][workingReplicate][workingTime]={}
                        if geneName not in transcriptome[workingGenotype][workingReplicate][workingTime].keys():
                            transcriptome[workingGenotype][workingReplicate][workingTime][geneName]=value                
    
    return transcriptome,listOfTranscriptomicsNames

##
## MAIN
##

# 0. user defined variables
stDataFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dv-mmp_20160830-L1vsL3.txt'
srDataFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dv-mmp_20160830-LS1vsLS3.txt'
figureDir='/Users/alomana/gDrive2/tmp/'
annotationFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dvh-protein-annotations.txt'
transcriptomeDataFile='/Users/alomana/gDrive2/projects/enigma/data/cufflinks_fpkms_matrix.v2.txt'
histogramFigureFile='/Users/alomana/gDrive2/tmp/histogramFigure.pdf'
correlatedGenesFile='/Users/alomana/gDrive2/tmp/correlated.txt'
anticorrelatedGenesFile='/Users/alomana/gDrive2/tmp/anticorrelated.txt'

genotypes=['WT','744']
conditionNames=['ST1','SR1','ST3','SR3']
transcriptReplicateNames=['2','4']
proteinReplicateNames=['1','2','3']


# 1. read data
print()
print('reading data...')

# 1.1. read identifiers
annotation=annotationReader()
print('annotation data read completed')

# 1.2. read proteomics data
proteome,listOfProteomicsNames=proteomeReader()
print('%s proteins quantified.'%len(listOfProteomicsNames))

# 1.3. read transcriptome data
transcriptome,listOfTranscriptomicsNames=transcriptomeReader()
print('%s transcripts quantified.'%len(listOfTranscriptomicsNames))

# 2. compute correlation
print()
print('computing correlations...')

# 2.1. find the intersection of both protein and transcript quantification
intersection=list(set(listOfProteomicsNames) & set(listOfTranscriptomicsNames))
print('%s elements with both transcript and protein quantification.'%(len(intersection)))

# 2.2. format for correlation computation
rhos=[]
threshold=10 # threshold for maximal expression
correlatedGenes={}
anticorrelatedGenes={}

for geneName in intersection:
    
    coefficients=[]
    lowAbundance=False
    print(geneName)

    for genotype in genotypes:
        
        # take the 3 transcript replicates
        transcriptValues=[]
        for replicate in transcriptReplicateNames:
            series=[]
            for conditionName in conditionNames:
                value=transcriptome[genotype][replicate][conditionName][geneName]
                series.append(value)
            if max(series) < threshold:
                lowAbundance=True
            transformedSeries=[numpy.log2((element+1)/(series[0]+1)) for element in series]
            transcriptValues.append(transformedSeries)
            print(transformedSeries)
                
        # take the 2 protein replicates
        proteinValues=[]
        for replicate in proteinReplicateNames:
            series=[]
            for conditionName in conditionNames:
                value=proteome[genotype][replicate][conditionName][geneName]
                series.append(value)
            transformedSeries=[element-series[0] for element in series]
            proteinValues.append(transformedSeries)
            print(series)
            print(transformedSeries)
        print()
            
        # compute the 6 correlation coefficients and take the median
        for seriesT in transcriptValues:
            for seriesP in proteinValues:
                pearsonC=scipy.stats.pearsonr(seriesT,seriesP)[0]
                coefficients.append(pearsonC)
                
    # adding the median value
    if lowAbundance is not True:
        average=numpy.median(coefficients)
        rhos.append(average)

        # saving gene names depending on their performance
        if average < -0.5:
            anticorrelatedGenes[geneName]=average
        if average > 0.5:
            correlatedGenes[geneName]=average

# 2.3. saving anticorrelated and correlated genes into files
print('saving into files anticorrelated and correlated genes...')
f=open(anticorrelatedGenesFile,'w')
f.write('Gene ID\tPearson Correlation Coefficient\n')
sortedGenes=sorted(anticorrelatedGenes,key=anticorrelatedGenes.__getitem__)
for element in sortedGenes:
    line=element+'\t'+str(anticorrelatedGenes[element])+'\n'
    f.write(line)
    print(element)
f.close()

f=open(correlatedGenesFile,'w')
f.write('Gene ID\tPearson Correlation Coefficient\n')
sortedGenes=sorted(correlatedGenes,key=correlatedGenes.__getitem__,reverse=True)
for element in sortedGenes:
    line=element+'\t'+str(correlatedGenes[element])+'\n'
    f.write(line)
f.close()  

# 3. make a figure with distribution of coefficients
print()
print('making a distribution figure...')
n,bins,tempo=matplotlib.pyplot.hist(rhos,bins=33,range=(-1,1),normed=True,color='black')
matplotlib.pyplot.xlim([-1.1,1.1])
#matplotlib.pyplot.ylim([-0.05,1.5])
matplotlib.pyplot.xlabel('rho')
matplotlib.pyplot.ylabel('p(rho)')

# 3.1. saving the figure to a file
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(histogramFigureFile)
matplotlib.pyplot.clf()


