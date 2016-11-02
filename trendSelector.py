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
    this function reads the annotation file and creates a dictionary between systematic names and protein names
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

def bootstrapperDown744(tempo):

    passed=0
    for randomProtein in range(len(listOfProteins)):
        trajectoryWT=random.choice(originalTrajectories)
        trajectory744=random.choice(originalTrajectories)
        timePoints=[1,1,1,2,2,2,3,3,3,4,4,4]
        flag=downRegulationChecker744(timePoints,trajectoryWT,trajectory744)
        if flag == True:
            passed=passed+1

    return passed

def bootstrapperUp(tempo):

    passed=0
    for randomProtein in range(len(listOfProteins)):
        trajectoryWT=random.choice(originalTrajectories)
        trajectory744=random.choice(originalTrajectories)
        timePoints=[1,1,1,2,2,2,3,3,3,4,4,4]
        flag=upregulationChecker(timePoints,trajectoryWT,trajectory744)
        if flag == True:
            passed=passed+1

    return passed

def bootstrappingConditionChecker(tempo):

    '''
    this function checks the 4 conditions again
    '''

    # f.1. compute how many proteins pass the filter
    passed=0
    for randomProtein in range(len(listOfProteins)):

        trajectoryWT=random.choice(originalTrajectories)
        trajectory744=random.choice(originalTrajectories)

        # checking the 4 rules
        success=0

        # f.1.1. condition 1: t-test on ST1 vs SR3
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

        # f.1.2. condition 2: slope < 1 fold change in range (-1/3)
        x=[1,1,1,2,2,2,3,3,3,4,4,4]
        y=trajectoryWT
        slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
        if slope < -(1/3):
            success=success+1

        # f.1.3. condition 3: slope slope > 1 fold change in range (-1/3)
        y=trajectory744
        slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
        if slope > -(1/3):
            success=success+1

        # f.1.4. condition 4: mean WT < mean 744 at SR3
        x=trajectoryWT[-3:]
        y=trajectory744[-3:]
        if numpy.median(x) < numpy.median(y):
            success=success+1

        # counting number of proteins that passed
        if success == 4:
            passed=passed+1

    return passed
    
def classifier(protein,verbose):

    '''
    this function classifies transcripts depending on several conditions
    '''

    evaluations=[None,None,None,None]
    grades=[None,None,None,None]
    
    # 1. condition 1: t-test on ST1 vs SR3
    genotype='WT'
    condition='ST1'
    x=[];y=[]
    for replicate in replicates:
        x.append(proteome[genotype][replicate][condition][protein])
    condition='SR3'
    for replicate in replicates:
        y.append(proteome[genotype][replicate][condition][protein])

    tempo,pvalueA=scipy.stats.shapiro(x)
    tempo,pvalueB=scipy.stats.shapiro(y)
    if min(pvalueA,pvalueB) < 0.05:
        statistic,pvalue=scipy.stats.mannwhitneyu(x,y)
    else:
        statistic,pvalue=scipy.stats.ttest_ind(x,y)

    grades[0]=pvalue
        
    if pvalue < 0.05:
        evaluations[0]=1
    else:
        evaluations[0]=0

    # 2. condition 2: slope < 1 fold change in range (-1/3)
    genotype='WT'
    x=[];y=[]
    for condition in ['ST1','ST3']:
        for replicate in replicates:
            x.append(int(condition[-1]))
            y.append(proteome[genotype][replicate][condition][protein])
    for condition in ['SR1','SR3']:
        for replicate in replicates:
            x.append(int(condition[-1])+1)
            y.append(proteome[genotype][replicate][condition][protein])

    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
    grades[1]=slope
    if slope < -(1/3):
        evaluations[1]=1
    else:
        evaluations[1]=0

    # 3. condition 3: slope slope > 1 fold change in range (-1/3)
    genotype='744'
    x=[];y=[]
    for condition in ['ST1','ST3']:
        for replicate in replicates:
            x.append(int(condition[-1]))
            y.append(proteome[genotype][replicate][condition][protein])
    for condition in ['SR1','SR3']:
        for replicate in replicates:
            x.append(int(condition[-1])+1)
            y.append(proteome[genotype][replicate][condition][protein])

    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
    grades[2]=slope
    if slope > -(1/3):
        evaluations[2]=1
    else:
        evaluations[2]=0

    # 4. condition 4: mean WT < mean 744 at SR3
    evaluations[3]=0
    grades[3]=1
    x=[];y=[]
    condition='SR3'
    genotype='WT'
    for replicate in replicates:
        x.append(proteome[genotype][replicate][condition][protein])
    genotype='744'
    for replicate in replicates:
        y.append(proteome[genotype][replicate][condition][protein])
    gap=min(y)-max(x)

    if numpy.median(x) < numpy.median(y):
        evaluations[3]=1
        grades[3]=gap
   
    return gap,evaluations,grades

def downRegulationChecker744(timePoints,trajectoryWT,trajectory744):

    '''
    this function return a binary for downregulation in the 744
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
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(timePoints,trajectory744)
    if slope < -(1/3):
        success=success+1

    # f.3. condition 3: slope > -1 fold change in range (-1/3) for WT
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(timePoints,trajectoryWT)
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

def expressionReader():

    '''
    this function reads the data and returns the expression in a dictionary format
    '''

    # f.1. checking that proteins were detected in both ST and SR
    STproteins=[]
    with open(stDataFile,'r') as f:
        next(f)
        for line in f:
            organism=line.split('\t')[0]
            proteinName=line.split('\t')[1]
            if organism == 'DESVH':
                STproteins.append(proteinName)
    print('%s proteins detected in ST conditions...'%len(STproteins))
    
    SRproteins=[]
    with open(srDataFile,'r') as f:
        next(f)
        for line in f:
            organism=line.split('\t')[0]
            proteinName=line.split('\t')[1]
            if organism == 'DESVH':
                SRproteins.append(proteinName)
    print('%s proteins detected in SR conditions...'%len(STproteins))

    commonProteins=[]
    for element in STproteins:
        if element in SRproteins:
            commonProteins.append(element)
    uniqueProteins=list(set(commonProteins))
    print('%s proteins found in both conditions...'%len(uniqueProteins))

    # f.2. building the expression variable
    proteome={}
    proteome=fileReader(stDataFile,uniqueProteins,proteome)
    proteome=fileReader(srDataFile,uniqueProteins,proteome)
                            
    return proteome,uniqueProteins

def fileReader(inputFileName,uniqueProteins,proteome):

    '''
    this function reads input files
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

                        # appending data
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

    matplotlib.pyplot.hist(distribution,bins=theBins,normed=True,color='black')
    matplotlib.pyplot.axvline(x=observedValue,linewidth=2,color='r',ls='--')
    matplotlib.pyplot.xlabel('count')
    matplotlib.pyplot.ylabel('p(count)')
    matplotlib.pyplot.title(label)

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(histogramFigureFile)
    matplotlib.pyplot.clf()

    return None

def panelGrapher(targets,label):

    '''
    this function makes a panel figure with target genes
    '''
    # f.1. dealing with essential genes panel
    figureFileName=figureDir+label+'.pdf'
    if label == 'essential':
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
            #print(targets[i],i,indexRow,indexCol)
            singleFigureGrapher(targets[i],axes,indexRow,indexCol,label)

            if (i+1)%10 == 0:
                fig.tight_layout()
                matplotlib.pyplot.savefig(figureFileName)
                matplotlib.pyplot.clf()
        
        # completing the last figure
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
    
    # f.1. 744
    genotype='744'
    theColor='#F67770'
    # ST1 --> ST3 
    x=[];y=[]
    for condition in ['ST1','ST3']:
        for replicate in replicates:
            x.append(int(condition[-1]))
            y.append(proteome[genotype][replicate][condition][protein])
    for condition in ['SR1','SR3']:
        for replicate in replicates:
            x.append(int(condition[-1])+1)
            y.append(proteome[genotype][replicate][condition][protein])
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
    x=numpy.array(x); y=numpy.array(y)
    model=slope*x+intercept
    cx=[x[i] for i in range(0,len(x),csize)]
    cy=[list(y[i:i+csize]) for i in range(0,len(y),csize)]
    bp=axes[indexRow,indexCol].boxplot(cy,positions=cx,patch_artist=True)
    setBoxColors(bp,theColor)
    axes[indexRow,indexCol].plot(x,model,'-',color=theColor,lw=3,label='744')

    # f.2. WT
    genotype='WT'
    theColor='#1FBFC3'
    # ST1 --> ST3
    conditions=['ST1','ST3']
    x=[];y=[]
    for condition in ['ST1','ST3']:
        for replicate in replicates:
            x.append(int(condition[-1]))
            y.append(proteome[genotype][replicate][condition][protein])
    for condition in ['SR1','SR3']:
        for replicate in replicates:
            x.append(int(condition[-1])+1)
            y.append(proteome[genotype][replicate][condition][protein])
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
    x=numpy.array(x); y=numpy.array(y)
    model=slope*x+intercept
    cx=[x[i] for i in range(0,len(x),csize)]
    cy=[list(y[i:i+csize]) for i in range(0,len(y),csize)]
    bp=axes[indexRow,indexCol].boxplot(cy,positions=cx,patch_artist=True)
    setBoxColors(bp,theColor)
    axes[indexRow,indexCol].plot(x,model,'-',color=theColor,lw=3,label='WT')

    foldChange=numpy.mean(cy[-1])-numpy.mean(cy[0])
    tempo,pvalueA=scipy.stats.shapiro(cy[0])
    tempo,pvalueB=scipy.stats.shapiro(cy[-1])
    if min(pvalueA,pvalueB) < 0.05:
        statistic,pvalue=scipy.stats.mannwhitneyu(x,y)
    else:
        statistic,pvalue=scipy.stats.ttest_ind(x,y)
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
        axes[indexRow,indexCol].set_xticklabels(['ST1','ST3','SR1','SR3'])
    elif indexRow == 0 and label == 'systematic' and protein in ['Q3V892','Q72FN6']:
        axes[indexRow,indexCol].set_xticklabels(['ST1','ST3','SR1','SR3'])
    else:
        axes[indexRow,indexCol].set_xticklabels([])

    if indexRow == 0:
        axes[indexRow,indexCol].legend(loc=3)

    return None

def slopeCheck(protein,expectedDistribution):

    '''
    this function checks the slope difference for a protein
    '''
    
    timePoints=[1,1,1,2,2,2,3,3,3,4,4,4]
    trajectories=[]
    
    for genotype in genotypes:
        trajectory=[]
        for condition in ['ST1','SR1','ST3','SR3']:
            for replicate in replicates:
                trajectory.append(proteome[genotype][replicate][condition][protein])
        trajectories.append(trajectory)

    # compute the delta slope between the WT and 744 and append to distribution
    trajectoryA=trajectories[0]
    trajectoryB=trajectories[1]

    slopeA,tempo,tempo,tempo,tempo=scipy.stats.linregress(timePoints,trajectoryA)
    slopeB,tempo,tempo,tempo,tempo=scipy.stats.linregress(timePoints,trajectoryB)

    slope=slopeA-slopeB

    # compute the index and print
    largerElements = [element for element in expectedDistribution if element > slope]
    pvalue=1-(len(largerElements)/len(expectedDistribution))
    if pvalue < 0.05:
        tag='*'
    else:
        tag=''
    message='\t'.join([annotation[protein],"{0:.4f}".format(slope),"{0:.4f}".format(pvalue),tag])

    print(message)

    return None

def upregulationChecker(timePoints,trajectoryWT,trajectory744):

    '''
    this function return a binary for upregulation in the WT
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
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(timePoints,trajectoryWT)
    if slope > 1/3:
        success=success+1

    # f.3. condition 3: slope slope < 1 fold change in range (1/3) for 744
    slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(timePoints,trajectory744)
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

def visualTargetsReader():

    '''
    this function reads the visual targets
    '''

    visualTargets=[]
    with open(visualTargetsFile,'r') as f:
        for line in f:
            vector=line.split('\t')
            proteinName=[proteinName for proteinName, systematicName in annotation.items() if systematicName == vector[0]][0]
            visualTargets.append(proteinName)

    return visualTargets

### MAIN

# 0. user defined functions
stDataFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dv-mmp_20160830-L1vsL3.txt'
srDataFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dv-mmp_20160830-LS1vsLS3.txt'
figureDir='figures/'
annotationFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dvh-protein-annotations.txt'
visualTargetsFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/proteins-with-trend.txt'
essentialityFile='/Users/alomana/gDrive2/projects/enigma/data/proteomics/dvu_essential-genes.txt'

genotypes=['WT','744']
replicates=['1','2','3']
timePoints=['ST1','SR1','ST3','SR3']

iterations=10000 # should be 10k

numberOfThreads=4

# 1. read data
print()
print('reading  data...')

# 1.1. read identifiers
annotation=annotationReader()

# 1.2. read expression data
proteome,listOfProteins=expressionReader()
listOfProteins.sort()
print('%s proteins found.'%len(listOfProteins))

# 1.3. read visual targets
visualTargets=visualTargetsReader()
print('%s visual targets identified.'%len(visualTargets))

# 1.4. reading essentiality
essentialGenes=essentialityReader()

# 2. classifying proteins
print()
print('classifying proteins ...')

found={}
gaps={}
notifications={}
for protein in listOfProteins:
    gap,evaluations,grades=classifier(protein,verbose=False)
    notifications[protein]=[evaluations,grades]
    gaps[protein]=gap
    if sum(evaluations) == 4:
        found[protein]=gap
systematicTargets=sorted(found,key=found.__getitem__,reverse=True)
print('%s systematic targets identified.'%len(systematicTargets))

# 3. compare systematic vs visual targets
print()
print('comparing 2 sets of targets...')
# 3.1. are visual targets a subset of visual targets?
intersection=list(set(visualTargets) & set(systematicTargets))
if len(systematicTargets) == len(intersection):
    print('systematic targets IS a subset of visual targets.')
else:
    print('systematic targets IS NOT a subset of visual targets.')
print()

# visual targets backup by systematic analysis
intersection=list(set(visualTargets) & set(systematicTargets))
print('%s visual targets in systematic targets:'%(len(intersection)))
for element in intersection:
    print(annotation[element])
print()

# missing targets in visual
new=[element for element in systematicTargets if element not in visualTargets]
print('%s systematic targets not in visual targets:'%len(new))
for element in new:
    print(annotation[element])
print()

# define visual targets failed...
failed=[element for element in visualTargets if element not in systematicTargets]

# 4. plotting
print('plotting...')

# 4.1. plotting all essential genes ordered by gap
print('plotting essential panel...')
selectedTargets=['Q72WJ8','Q72CT6','Q72WF7','Q72A54','Q726R4','Q72AR0','Q72FA9','Q725Z7','Q727Q1','Q72AF9']
panelGrapher(selectedTargets,'essential')

# 4.2. plotting all genes ordered by gap
print('plotting systematic targets...')
selectedGaps={key:value for key,value in gaps.items() if key in systematicTargets}
sortedTargets=sorted(selectedGaps,key=selectedGaps.__getitem__,reverse=True)
panelGrapher(sortedTargets,'systematic')

# 5. bootstrapping
print()
print('bootstrapping for downregulation...')

# 5.1. storing original information
originalTrajectories=[]
for protein in listOfProteins:
    for genotype in genotypes:
        y=[]
        for condition in ['ST1','SR1','ST3','SR3']:
            for replicate in replicates:
                y.append(proteome[genotype][replicate][condition][protein])
        originalTrajectories.append(y)

# 5.2. for n proteomes
hydra=multiprocessing.pool.Pool(numberOfThreads)
distribution=[None for i in range(iterations)]
distribution=hydra.map(bootstrappingConditionChecker,distribution)
    
# 5.4 compute the p-value of finding systematic targets
distribution.sort()
print(distribution[:100])
largeElements=[element for element in distribution if element >= len(systematicTargets)]
pvalue=len(largeElements)/len(distribution)
print('%s samplings with more cases than found. P-value: %s'%(len(largeElements),pvalue))

label='WT down regulation'
histogramFigureFile=figureDir+'permutation.WT.downregulation.pdf'
histogramBuilder(len(systematicTargets),distribution,label,histogramFigureFile)

# 6. bootstrapping on the opposite pattern
print()
print('bootstrapping for upregulation...')
          
# 6.1. finding the number of cases that pass opposite pattern
allFlags=[]
for protein in listOfProteins:

    timePoints=[1,1,1,2,2,2,3,3,3,4,4,4]
    # selecting the trajectories for WT
    trajectoryWT=[]
    genotype='WT'
    for condition in ['ST1','SR1','ST3','SR3']:
        for replicate in replicates:
            trajectoryWT.append(proteome[genotype][replicate][condition][protein])

    # selecting the trajectories for 744
    trajectory744=[]
    genotype='744'
    for condition in ['ST1','SR1','ST3','SR3']:
        for replicate in replicates:
            trajectory744.append(proteome[genotype][replicate][condition][protein])

    # checking for pattern
    flag=upregulationChecker(timePoints,trajectoryWT,trajectory744)
    allFlags.append(flag)
numberOfUpregulated=sum(allFlags)
print('%s upregulated proteins found'%(numberOfUpregulated))
    
# 6.2. finding the distribution of cases from a shuffled bag
distribution=[None for element in range(iterations)]
distribution=hydra.map(bootstrapperUp,distribution)

# 6.3. computing the p-value
distribution.sort()
print(distribution[:100])
largeElements=[element for element in distribution if element >= numberOfUpregulated]
pvalue=len(largeElements)/len(distribution)
print('%s samplings with more cases than found. P-value: %s'%(len(largeElements),pvalue))

label='WT up regulation'
histogramFigureFile=figureDir+'permutation.WT.upregulation.pdf'
histogramBuilder(numberOfUpregulated,distribution,label,histogramFigureFile)

# 7. Check that the pattern of downregulated proteins in 744. Is this more than chance?
print()
print('bootstrapping for downregulation in 744...')

# 7.1. finding the number of cases for that pattern
allFlags=[]
for protein in listOfProteins:

    timePoints=[1,1,1,2,2,2,3,3,3,4,4,4]
    # selecting the trajectories for WT
    trajectoryWT=[]
    genotype='WT'
    for condition in ['ST1','SR1','ST3','SR3']:
        for replicate in replicates:
            trajectoryWT.append(proteome[genotype][replicate][condition][protein])

    # selecting the trajectories for 744
    trajectory744=[]
    genotype='744'
    for condition in ['ST1','SR1','ST3','SR3']:
        for replicate in replicates:
            trajectory744.append(proteome[genotype][replicate][condition][protein])

    # checking for pattern
    flag=downRegulationChecker744(timePoints,trajectoryWT,trajectory744)
    allFlags.append(flag)
numberOfDownregulated744=sum(allFlags)
print('%s downregulated 744 proteins found'%(numberOfDownregulated744))

# 7.2. finding the distribution of cases from a shuffled bag
distribution=[None for element in range(iterations)]
distribution=hydra.map(bootstrapperDown744,distribution)

# 7.3. computing the p-value
distribution.sort()
print(distribution[:100])
largeElements=[element for element in distribution if element >= numberOfDownregulated744]
pvalue=len(largeElements)/len(distribution)
print('%s samplings with more cases than found. P-value: %s'%(len(largeElements),pvalue))

label='744 down regulation'
histogramFigureFile=figureDir+'permutation.744.dowregulation.pdf'
histogramBuilder(numberOfDownregulated744,distribution,label,histogramFigureFile)

# 8. final message
print('... all done.')
