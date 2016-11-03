###
### this script reads protein quantification data and select the genes that have specific trends.
### also it evaluates through permutations if the number of found proteins for a specific trend is expected by chance. 
###

import sys,numpy,random,copy
import scipy,scipy.stats
import matplotlib,matplotlib.pyplot
import multiprocessing, multiprocessing.pool

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

def bootstrapWTdown(tempo):

    '''
    this function is a wrapper to quantify the number of patterns (WT down) found in a random proteome
    '''

    passed=0
    for randomProtein in range(len(listOfProteins)):
        trajectoryWT=random.choice(originalTrajectories)
        trajectory744=random.choice(originalTrajectories)
        flag,gap=checkWTdown(sampleTimePoints,trajectoryWT,trajectory744)
        if flag == True:
            passed=passed+1

    return passed

def bootstrapWTup(tempo):

    '''
    this function is a wrapper to quantify the number of patterns (WT up) found in a random proteome
    '''

    passed=0
    for randomProtein in range(len(listOfProteins)):
        trajectoryWT=random.choice(originalTrajectories)
        trajectory744=random.choice(originalTrajectories)
        flag=checkWTup(sampleTimePoints,trajectoryWT,trajectory744)
        if flag == True:
            passed=passed+1

    return passed

def bootstrap744down(tempo):

    '''
    this function is a wrapper to quantify the number of patterns (744 down) found in a random proteome
    '''

    passed=0
    for randomProtein in range(len(listOfProteins)):
        trajectoryWT=random.choice(originalTrajectories)
        trajectory744=random.choice(originalTrajectories)
        flag=check744down(sampleTimePoints,trajectoryWT,trajectory744)
        if flag == True:
            passed=passed+1

    return passed

def checkWTdown(sampleTimePoints,trajectoryWT,trajectory744):

    '''
    this function returns a binary for downregulation in the WT
    '''

    success=0

    # f.1. condition 1: checking that there is significant regulation in WT ST1 vs SR3
    x=trajectoryWT[:3]
    y=trajectoryWT[-3:]

    tempo,pvalueA=scipy.stats.shapiro(x)
    tempo,pvalueB=scipy.stats.shapiro(y)
    if min(pvalueA,pvalueB) < 0.05:
        statistic,pvalue=scipy.stats.mannwhitneyu(x,y)
    else:
        statistic,pvalue=scipy.stats.ttest_ind(x,y)
    if pvalue < 0.05:
        success=success+1

    # f.2. condition 2: slope > 1 fold change in range (1/3) for WT
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(sampleTimePoints,trajectoryWT)
    if slope < -(1/3):
        success=success+1

    # f.3. condition 3: slope slope < 1 fold change in range (1/3) for 744
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(sampleTimePoints,trajectory744)
    if slope > -(1/3):
        success=success+1

    # f.4. condition 4: median WT < median 744 at SR3
    x=trajectoryWT[-3:]
    y=trajectory744[-3:]
    if numpy.median(x) < numpy.median(y):
        success=success+1

    if success == 4:
        flag=True
    else:
        flag=False

    # f.5. computing the gap for ranking
    gap=min(y)-max(x)

    return flag,gap
    
def checkWTup(sampleTimePoints,trajectoryWT,trajectory744):

    '''
    this function returns a binary for upregulation in the WT
    '''

    success=0

    # f.1. condition 1: checking that there is significant upregulation in WT ST1 vs SR3
    x=trajectoryWT[:3]
    y=trajectoryWT[-3:]

    tempo,pvalueA=scipy.stats.shapiro(x)
    tempo,pvalueB=scipy.stats.shapiro(y)
    if min(pvalueA,pvalueB) < 0.05:
        statistic,pvalue=scipy.stats.mannwhitneyu(x,y)
    else:
        statistic,pvalue=scipy.stats.ttest_ind(x,y)
    if pvalue < 0.05:
        success=success+1

    # f.2. condition 2: slope > 1 fold change in range (1/3) for WT
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(sampleTimePoints,trajectoryWT)
    if slope > 1/3:
        success=success+1

    # f.3. condition 3: slope slope < 1 fold change in range (1/3) for 744
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(sampleTimePoints,trajectory744)
    if slope < 1/3:
        success=success+1

    # f.4. condition 4: median WT > median 744 at SR3
    x=trajectoryWT[-3:]
    y=trajectory744[-3:]
    if numpy.median(x) > numpy.median(y):
        success=success+1

    if success == 4:
        flag=True
    else:
        flag=False

    return flag

def check744down(sampleTimePoints,trajectoryWT,trajectory744):

    '''
    this function returns a binary for downregulation in the 744
    '''

    success=0

    # f.1. condition 1: checking that there is significant downregulation in 744 ST1 vs SR3
    x=trajectory744[:3]
    y=trajectory744[-3:]

    tempo,pvalueA=scipy.stats.shapiro(x)
    tempo,pvalueB=scipy.stats.shapiro(y)
    if min(pvalueA,pvalueB) < 0.05:
        statistic,pvalue=scipy.stats.mannwhitneyu(x,y)
    else:
        statistic,pvalue=scipy.stats.ttest_ind(x,y)
    if pvalue < 0.05:
        success=success+1

    # f.2. condition 2: slope < -1 fold change in range (-1/3) for 744
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(sampleTimePoints,trajectory744)
    if slope < -(1/3):
        success=success+1

    # f.3. condition 3: slope > -1 fold change in range (-1/3) for WT
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(sampleTimePoints,trajectoryWT)
    if slope > -(1/3):
        success=success+1

    # f.4. condition 4: median WT > median 744 at SR3
    x=trajectoryWT[-3:]
    y=trajectory744[-3:]
    if numpy.median(x) > numpy.median(y):
        success=success+1

    if success == 4:
        flag=True
    else:
        flag=False

    return flag

def essentialityReader():

    '''
    This file returns a list of essential genes
    '''

    essentialGenes=[]
    with open(essentialityFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split()
            dvuName=vector[0]
            
            proteinName=None
            try:
                proteinName=list(annotation.keys())[list(annotation.values()).index(dvuName)]
            except:
                pass
            if proteinName != None and proteinName not in essentialGenes:
                essentialGenes.append(proteinName)

    return essentialGenes

def expressionFileReader(inputFileName,uniqueProteins,proteome):

    '''
    this function reads expression input files
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
                        geneName=vector[1]
                        value=float(vector[i])

                        # associating data to a variable
                        if geneName in uniqueProteins:
                            if workingGenotype not in proteome.keys():
                                proteome[workingGenotype]={}
                            if workingReplicate not in proteome[workingGenotype].keys():
                                proteome[workingGenotype][workingReplicate]={}
                            if workingTime not in proteome[workingGenotype][workingReplicate].keys():
                                proteome[workingGenotype][workingReplicate][workingTime]={}
                            if geneName not in proteome[workingGenotype][workingReplicate][workingTime].keys():
                                proteome[workingGenotype][workingReplicate][workingTime][geneName]=value

    return proteome

def expressionReader():

    '''
    this function reads the protein expression data and returns a dictionary
    '''

    # f.1. checking that proteins were detected in both ST and SR
    STproteins=getNames(stDataFile)
    print('%s proteins detected in ST conditions...'%len(STproteins))
    
    SRproteins=getNames(srDataFile)
    print('%s proteins detected in SR conditions...'%len(STproteins))

    commonProteins=[]
    for element in STproteins:
        if element in SRproteins:
            commonProteins.append(element)
    uniqueProteins=list(set(commonProteins))
    print('%s proteins found in both conditions...'%len(uniqueProteins))

    # f.2. building the expression variable
    proteome={}
    proteome=expressionFileReader(stDataFile,uniqueProteins,proteome)
    proteome=expressionFileReader(srDataFile,uniqueProteins,proteome)
                            
    return proteome,uniqueProteins

def getNames(fileName):

    '''
    this function reads the expression file and return the protein names
    '''

    proteinNames=[]
    with open(fileName,'r') as f:
        next(f)
        for line in f:
            organism=line.split('\t')[0]
            proteinName=line.split('\t')[1]
            if organism == 'DESVH':
                proteinNames.append(proteinName)

    return proteinNames

def histogramBuilder(observedValue,distribution,label,histogramFigureFile):

    '''
    this function creates a histogram from the permutation analysis
    '''

    # preparing bins
    top=max([max(distribution),observedValue])
    bottom=min([min(distribution),observedValue])
    interval=top-bottom
    epsilon=int(0.1*interval)
    theBins=list(range(bottom-epsilon,top+epsilon))

    # making the figure
    matplotlib.pyplot.hist(distribution,bins=theBins,normed=True,color='black')
    matplotlib.pyplot.axvline(x=observedValue,linewidth=2,color='r',ls='--')
    matplotlib.pyplot.xlabel('count')
    matplotlib.pyplot.ylabel('p(count)')
    matplotlib.pyplot.title(label)

    # saving the figure to a file
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(histogramFigureFile)
    matplotlib.pyplot.clf()

    return None

def panelGrapher(targets,label):

    '''
    this function makes a panel figure with target genes
    '''
    
    # f.1. dealing with representative genes panel
    figureFileName=figureDir+label+'.pdf'
    if label == 'representative':
        fig,axes=matplotlib.pyplot.subplots(nrows=5,ncols=2,figsize=(8.5,11))
        for i in range(len(targets)):
            indexRow=int(i/2)
            indexCol=i-2*indexRow
            singleFigureGrapher(targets[i],axes,indexRow,indexCol,label)
        # completing the figure
        fig.tight_layout()
        matplotlib.pyplot.savefig(figureFileName)
        matplotlib.pyplot.clf()

    # f.2. dealing with systematic target genes
    if label == 'systematic':
        for i in range(len(targets)):
            if i%10 == 0:
                figureFileName=figureDir+label+'.'+str(int(i/10)).zfill(2)+'.pdf'
                fig,axes=matplotlib.pyplot.subplots(nrows=5,ncols=2,figsize=(8.5,11))

            localIndex=i-int(i/10)*10
            indexRow=int(localIndex/2)
            indexCol=localIndex-2*indexRow
            singleFigureGrapher(targets[i],axes,indexRow,indexCol,label)

            if (i+1)%10 == 0:
                fig.tight_layout()
                matplotlib.pyplot.savefig(figureFileName)
                matplotlib.pyplot.clf()
        
        # completing last page of figures
        fig.tight_layout()
        for index1 in range(1,5):
            for index2 in range(1,3):
                axes[index1,index2-1].axis('off')
        matplotlib.pyplot.savefig(figureFileName)
        matplotlib.pyplot.clf()
    
    return None

def setBoxColors(bp,theColor):

    '''
    this function access the elements of a boxplot and colors them appropriately
    '''

    matplotlib.pyplot.setp(bp['boxes'],color=theColor,alpha=0.33)
    matplotlib.pyplot.setp(bp['caps'],color='None')
    matplotlib.pyplot.setp(bp['whiskers'],color=theColor,ls='-',alpha=0.33)
    matplotlib.pyplot.setp(bp['fliers'],markeredgecolor=theColor,marker='+')
    matplotlib.pyplot.setp(bp['medians'],color=theColor,alpha=0.33)

    return None

def singleFigureGrapher(protein,axes,indexRow,indexCol,label):

    '''
    this function builds single panels of a figure
    '''
    
    # f.0. variable to get trios
    csize=3
    x=numpy.array(sampleTimePoints)
    
    # f.1. 744
    genotype='744'
    theColor='#F67770'
    trajectory744=[]
    for condition in conditionNames:
        for replicate in replicates:
            trajectory744.append(proteome[genotype][replicate][condition][protein])
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,trajectory744)
    model=slope*x+intercept
    cx=[x[i] for i in range(0,len(x),csize)]
    cy=[list(trajectory744[i:i+csize]) for i in range(0,len(trajectory744),csize)]
    bp=axes[indexRow,indexCol].boxplot(cy,positions=cx,patch_artist=True)
    setBoxColors(bp,theColor)
    axes[indexRow,indexCol].plot(x,model,'-',color=theColor,lw=3,label='744')

    # f.2. WT
    genotype='WT'
    theColor='#1FBFC3'
    trajectoryWT=[]
    for condition in conditionNames:
        for replicate in replicates:
            trajectoryWT.append(proteome[genotype][replicate][condition][protein])
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,trajectoryWT)
    model=slope*x+intercept
    cx=[sampleTimePoints[i] for i in range(0,len(x),csize)]
    cy=[list(trajectoryWT[i:i+csize]) for i in range(0,len(trajectoryWT),csize)]
    bp=axes[indexRow,indexCol].boxplot(cy,positions=cx,patch_artist=True)
    setBoxColors(bp,theColor)
    axes[indexRow,indexCol].plot(x,model,'-',color=theColor,lw=3,label='WT')

    begin=cy[0]
    end=cy[-1]
    foldChange=numpy.mean(end)-numpy.mean(begin)
    tempo,pvalueA=scipy.stats.shapiro(begin)
    tempo,pvalueB=scipy.stats.shapiro(end)
    if min(pvalueA,pvalueB) < 0.05:
        statistic,pvalue=scipy.stats.mannwhitneyu(begin,end)
    else:
        statistic,pvalue=scipy.stats.ttest_ind(begin,end)

    if pvalue > 0.05:
        formattedSignificance=''
    elif pvalue <= 0.05 and pvalue > 0.01:
        formattedSignificance='*'
    elif pvalue < 0.01:
        formattedSignificance='**'
    else:
        print('error from plotter at significances')

    # f.3. making the figures nicer
    formattedFC="{0:.2f}".format(foldChange)
    title=annotation[protein]+', log$_2$FC='+formattedFC+', '+formattedSignificance
    if protein in essentialGenes:
        titleColor='red'
    else:
        titleColor='black'
    axes[indexRow,indexCol].set_title('%s'%title,color=titleColor)
    
    matplotlib.pyplot.xlim([0.5,4.5])
    if indexRow == 4:
        axes[indexRow,indexCol].set_xticklabels(conditionNames)
    elif indexRow == 0 and label == 'systematic' and protein in ['Q3V892','Q72FN6']:
        axes[indexRow,indexCol].set_xticklabels(conditionNames)
    else:
        axes[indexRow,indexCol].set_xticklabels([])

    if indexRow == 0:
        axes[indexRow,indexCol].legend(loc=3)

    return None

###
### MAIN
###

print()
print('WELCOME TO PROTEOME TREND ANALYSIS.')

# 0. user defined functions
stDataFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dv-mmp_20160830-L1vsL3.txt'
srDataFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dv-mmp_20160830-LS1vsLS3.txt'
figureDir='/Users/alomana/gDrive2/tmp/'
annotationFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dvh-protein-annotations.txt'
visualTargetsFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/proteins-with-trend.txt'
essentialityFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dvu_essential-genes.txt'

genotypes=['WT','744']
replicates=['1','2','3']
conditionNames=['ST1','SR1','ST3','SR3']

iterations=100 # for the permutation analysis. 10,000 is a good size
numberOfThreads=4
sampleTimePoints=[1,1,1,2,2,2,3,3,3,4,4,4]

# 1. read data
print()
print('reading  data...')

# 1.1. read identifiers
annotation=annotationReader()

# 1.2. read expression data
proteome,listOfProteins=expressionReader()
listOfProteins.sort()
print('%s proteins found.'%len(listOfProteins))

# 1.3. reading essentiality
essentialGenes=essentialityReader()

# 2. classifying proteins
print()
print('classifying proteins...')

allFlags=[]
lastTimeGap={}
for protein in listOfProteins:
    # selecting the trajectories for WT
    trajectoryWT=[]
    genotype='WT'
    for condition in conditionNames:
        for replicate in replicates:
            trajectoryWT.append(proteome[genotype][replicate][condition][protein])

    # selecting the trajectories for 744
    trajectory744=[]
    genotype='744'
    for condition in conditionNames:
        for replicate in replicates:
            trajectory744.append(proteome[genotype][replicate][condition][protein])

    # checking for pattern
    flag,gap=checkWTdown(sampleTimePoints,trajectoryWT,trajectory744)
    if flag == True:
        lastTimeGap[protein]=gap
systematicTargets=sorted(lastTimeGap,key=lastTimeGap.__getitem__,reverse=True)
print('%s systematic targets identified.'%len(systematicTargets))

# 3. plotting
print()
print('plotting...')

# 3.1. plotting representative genes ordered by gap
print('plotting representative panel...')
selectedTargets=['Q72WJ8','Q72CT6','Q72WF7','Q72A54','Q726R4','Q72AR0','Q72FA9','Q725Z7','Q727Q1','Q72AF9']
panelGrapher(selectedTargets,'representative')

# 3.2. plotting all genes ordered by gap
print('plotting systematic targets...')
panelGrapher(systematicTargets,'systematic')

# 4. bootstrapping for WT downregulation
print()
print('bootstrapping for WT downregulation...')

# 4.1. storing original information
originalTrajectories=[]
for protein in listOfProteins:
    for genotype in genotypes:
        y=[]
        for condition in conditionNames:
            for replicate in replicates:
                y.append(proteome[genotype][replicate][condition][protein])
        originalTrajectories.append(y)

# 4.2. sampling n proteomes
hydra=multiprocessing.pool.Pool(numberOfThreads)
distribution=[None for i in range(iterations)]
distribution=hydra.map(bootstrapWTdown,distribution)
    
# 4.3. compute the p-value of finding systematic targets
largeElements=[element for element in distribution if element >= len(systematicTargets)]
pvalue=len(largeElements)/len(distribution)
print('%s samplings with more cases than found. p-value: %s'%(len(largeElements),pvalue))

# 4.4. making a histogram plot out of the distribution
label='WT down regulation'
histogramFigureFile=figureDir+'permutation.WT.down.pdf'
histogramBuilder(len(systematicTargets),distribution,label,histogramFigureFile)

# 5. bootstrapping for WT upregulation
print()
print('bootstrapping for WT upregulation...')
          
# 5.1. finding the number of cases that pass opposite pattern
allFlags=[]
for protein in listOfProteins:

    # selecting the trajectories for WT
    trajectoryWT=[]
    genotype='WT'
    for condition in conditionNames:
        for replicate in replicates:
            trajectoryWT.append(proteome[genotype][replicate][condition][protein])

    # selecting the trajectories for 744
    trajectory744=[]
    genotype='744'
    for condition in conditionNames:
        for replicate in replicates:
            trajectory744.append(proteome[genotype][replicate][condition][protein])

    # checking for pattern
    flag=checkWTup(sampleTimePoints,trajectoryWT,trajectory744)
    allFlags.append(flag)
numberOfUpregulated=sum(allFlags)
print('%s upregulated proteins found in WT.'%(numberOfUpregulated))
    
# 5.2. finding the distribution of cases from a shuffled bag
distribution=[None for element in range(iterations)]
distribution=hydra.map(bootstrapWTup,distribution)

# 5.3. computing the p-value
largeElements=[element for element in distribution if element >= numberOfUpregulated]
pvalue=len(largeElements)/len(distribution)
print('%s samplings with more cases than found. p-value: %s'%(len(largeElements),pvalue))

# 5.4. making a histogram plot out of the distribution
label='WT up regulation'
histogramFigureFile=figureDir+'permutation.WT.up.pdf'
histogramBuilder(numberOfUpregulated,distribution,label,histogramFigureFile)

# 6. bootstrapping for 744 downregulation
print()
print('bootstrapping for 744 downregulation...')

# 6.1. finding the number of cases for that pattern
allFlags=[]
for protein in listOfProteins:

    # selecting the trajectories for WT
    trajectoryWT=[]
    genotype='WT'
    for condition in conditionNames:
        for replicate in replicates:
            trajectoryWT.append(proteome[genotype][replicate][condition][protein])

    # selecting the trajectories for 744
    trajectory744=[]
    genotype='744'
    for condition in conditionNames:
        for replicate in replicates:
            trajectory744.append(proteome[genotype][replicate][condition][protein])

    # checking for pattern
    flag=check744down(sampleTimePoints,trajectoryWT,trajectory744)
    allFlags.append(flag)
numberOfDownregulated744=sum(allFlags)
print('%s downregulated proteins found in 744.'%(numberOfDownregulated744))

# 6.2. finding the distribution of cases from a shuffled bag
distribution=[None for element in range(iterations)]
distribution=hydra.map(bootstrap744down,distribution)

# 6.3. computing the p-value
largeElements=[element for element in distribution if element >= numberOfDownregulated744]
pvalue=len(largeElements)/len(distribution)
print('%s samplings with more cases than found. p-value: %s'%(len(largeElements),pvalue))

# 6.4. making a histogram plot out of the distribution
label='744 down regulation'
histogramFigureFile=figureDir+'permutation.744.down.pdf'
histogramBuilder(numberOfDownregulated744,distribution,label,histogramFigureFile)

# 7. final message
print()
print('... ALL DONE.')
