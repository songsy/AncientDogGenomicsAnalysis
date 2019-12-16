#!/usr/bin/env python
# -*- coding: ASCII -*-

###This program provides GLs from biallelic SNPs from bam files while taking into account post mortem damage as estimate by MapDamage.
###Otherwise the algorithm is the same as GATK Unified Genotype for diploid calls
###one file is created:
###An GL file noting the genotype likelihood in the bed file
###The bed file should have the required alleles
##aDNA_GenoCaller <indexed bamfile> <bed file> <reference genome> <5C-T mapdamage file> <3G-A mapdamage file>


from sys import argv
import pysam
import math
import numpy as np
import string
import numpy.ma as ma
import random
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import fmin
from scipy.optimize import fminbound
from scipy.optimize import curve_fit
import time
from random import randint



filenamein=argv[1]
filenameinB=argv[2]
ref_file=argv[3]
mapdamageCT=argv[4]
mapdamageGA=argv[5]

min_RD=1 #int(argv[3])
MQ=15 #float(argv[4])
BQ=15 #float(argv[5])
theta=0.001


if filenamein[-4:] == '.bam':
    filenameout=string.split(filenamein,'.bam')[0]
else:
    filenameout=filenamein[:]

    
if '/' in filenamein:
    filenameout=string.split(filenameout,'/')[-1]

#plotfile=filenameout+'.aDNA.snp.expMDfit_WB.pdf'
#MDfile=filenameout+'.aDNA.snp.fitted_model_params_WB'
filenameout1=filenameout+'.'+filenameinB+'.aDNA.snp.gls'


def phred2prob(x):
    return 10.0**(-x/10.0)

def prob2phred(x):
    return -10*math.log10(x)

#exponential function
def exp_fit(x,a,b,c):
    return a*np.exp(-x*b)+c

#stretched exponential (weibell) function
def weibull_fit(x,a,b,c):
    return a*np.exp(-(x**c)*b)

#genotype called that incorporates damage for aDNA
def geno_caller_10GT_aDNA(X):
    GL=np.zeros(10)   #all 10 possible genotypes and order = AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
    hap=np.zeros((len(X),4))  #all 4 haploid possibilities, A,C,G,T
    all_dic={}
    all_dic['A']=0
    all_dic['C']=1
    all_dic['G']=2
    all_dic['T']=3
    count=0
    
    for g in range(len(X)):
        if all_dic.has_key(X[g][0])==False:
            continue
        err=phred2prob(X[g][1])
        hap[g]=err/3.0
        if X[g][0]=='A':
            hap[g][0]=1-err
            hap[g][2]=((1-X[g][3])*(err/3.0))+(X[g][3]*(1-err))
        elif X[g][0]=='C':
            hap[g][1]= ((1-X[g][2])*(1-err))+(X[g][2]*(err/3.0))
        elif X[g][0]=='G':
             hap[g][2]= ((1-X[g][3])*(1-err))+(X[g][3]*(err/3.0))           
        elif X[g][0]=='T':
            hap[g][3]=1-err
            hap[g][1]=((1-X[g][2])*(err/3.0))+(X[g][2]*(1-err))                   

        GL[0]=GL[0]+math.log10(hap[g][0])                      
        GL[1]=GL[1]+math.log10((hap[g][0]+hap[g][1])/2)       
        GL[2]=GL[2]+math.log10((hap[g][0]+hap[g][2])/2)
        GL[3]=GL[3]+math.log10((hap[g][0]+hap[g][3])/2)
        GL[4]=GL[4]+math.log10(hap[g][1])
        GL[5]=GL[5]+math.log10((hap[g][1]+hap[g][2])/2)
        GL[6]=GL[6]+math.log10((hap[g][1]+hap[g][3])/2)
        GL[7]=GL[7]+math.log10(hap[g][2])
        GL[8]=GL[8]+math.log10((hap[g][2]+hap[g][3])/2)
        GL[9]=GL[9]+math.log10(hap[g][3])
        count+=1

    if count==0:
        GL.fill(-9)
    return GL






#open up mapdamage files and put results in an array
fileCT=open(mapdamageCT,'r')
CTdata=fileCT.read()
fileCT.close()
CTdata=string.split(CTdata,'\n')

if CTdata[-1]=='':
    del(CTdata[-1])


fileGA=open(mapdamageGA,'r')
GAdata=fileGA.read()
fileGA.close()
GAdata=string.split(GAdata,'\n')

if GAdata[-1]=='':
    del(GAdata[-1])

X_CT=[]
X_GA=[]

Y_CT=[]
Y_GA=[]
for g in range(1,len(CTdata)):
    k=string.split(CTdata[g],'\t')
    X_CT.append(float(k[0]))    
    Y_CT.append(float(k[1]))
    k=string.split(GAdata[g],'\t')
    X_GA.append(float(k[0])) 
    Y_GA.append(float(k[1]))


X_CT=np.asarray(X_CT)
Y_CT=np.asarray(Y_CT)
X_GA=np.asarray(X_GA)
Y_GA=np.asarray(Y_GA)


try:
    #Perform curve fitting of MapDamage results to an exponential function
    print '\nFitting the following weibull model to the MapDamage data : a*exp(-(x^c)*b)'
    fitParams_CT = curve_fit(weibull_fit, X_CT, Y_CT)
    CTa=fitParams_CT[0][0]
    CTb=fitParams_CT[0][1]
    CTc=fitParams_CT[0][2]

    fitParams_GA = curve_fit(weibull_fit, X_GA, Y_GA)
    GAa=fitParams_GA[0][0]
    GAb=fitParams_GA[0][1]
    GAc=fitParams_GA[0][2]

    out='Fit the following weibull model to the MapDamage data : a*exp(-(x^c)*b)\n\n'
    out=out+'Fitted coefficients for CT changes:\na\t'+str(CTa)+'\nb\t'+str(CTb)+'\nc\t'+str(CTc)+'\n'
    out=out+'\nFitted coefficients for GA changes:\na\t'+str(GAa)+'\nb\t'+str(GAb)+'\nc\t'+str(GAc)+'\n'
    print '\n'+out


##    #print coefficient inference to file
##    fileMD=open(MDfile,'w')
##    fileMD.write(out)
##    fileMD.close()
##
##    print '\nWrote fitted coefficients to '+MDfile
##
##    #set up plot area
##    rcParams['figure.figsize'] = 10, 6 
##    plt.ylabel('Substitution Frequency', fontsize = 16)
##    plt.xlabel('Read position', fontsize = 16)
##    plt.xlim(0,26)
##
##    #plot points and fitted curve
##    plt.plot(X_CT,Y_CT,'ro')
##    plt.plot(X_GA,Y_GA,'bo')
##    plt.plot(X_CT,weibull_fit(X_CT, fitParams_CT[0][0], fitParams_CT[0][1], fitParams_CT[0][2]),'r',ms=10,linewidth=2.0,label='C>T')
##    plt.plot(X_GA,weibull_fit(X_GA, fitParams_GA[0][0], fitParams_GA[0][1], fitParams_GA[0][2]),'b',ms=10,linewidth=2.0,label='G>A')
##    plt.legend()
##
##
##    # save plot to a file
##    plt.savefig(plotfile, bbox_inches=0, dpi=600)
##    plt.close()
##    print '\nPlotted MapDamage curve fit to '+plotfile

    #Make a reference list to determine how much to adjust a particular quality score given it's read position (maxed at 300 here)
    CT_decay=[]
    GA_decay=[]

    for g in range(300):
        CTpoint=weibull_fit(g+1,CTa,CTb,CTc)-theta
        GApoint=weibull_fit(g+1,GAa,GAb,GAc)-theta
        if CTpoint<0.0:
            CTpoint=0.0
        if GApoint<0.0:
            GApoint=0.0
        CT_decay.append(CTpoint)
        GA_decay.append(GApoint)
#in case weibull does not fit (occurs for low damage) use exponential fit
except:
    #Perform curve fitting of MapDamage results to an exponential function
    print '\nCould not fit weibull'
    print '\nFitting the following exponential model to the MapDamage data : a*exp(-x*b)+c'
    fitParams_CT = curve_fit(exp_fit, X_CT, Y_CT)
    CTa=fitParams_CT[0][0]
    CTb=fitParams_CT[0][1]
    CTc=fitParams_CT[0][2]

    fitParams_GA = curve_fit(exp_fit, X_GA, Y_GA)
    GAa=fitParams_GA[0][0]
    GAb=fitParams_GA[0][1]
    GAc=fitParams_GA[0][2]


    out='Fit the following exponential model to the MapDamage data : a*exp(-x*b)+c\n\n'
    out=out+'Fitted coefficients for CT changes:\na\t'+str(CTa)+'\nb\t'+str(CTb)+'\nc\t'+str(CTc)+'\n'
    out=out+'\nFitted coefficients for GA changes:\na\t'+str(GAa)+'\nb\t'+str(GAb)+'\nc\t'+str(GAc)+'\n'
    print '\n'+out


##    #print coefficient inference to file
##    fileMD=open(MDfile,'w')
##    fileMD.write(out)
##    fileMD.close()
##
##    print '\nWrote fitted coefficients to '+MDfile
##
##    #set up plot area
##    rcParams['figure.figsize'] = 10, 6 
##    plt.ylabel('Substitution Frequency', fontsize = 16)
##    plt.xlabel('Read position', fontsize = 16)
##    plt.xlim(0,26)
##
##    #plot points and fitted curve
##    plt.plot(X_CT,Y_CT,'ro')
##    plt.plot(X_GA,Y_GA,'bo')
##    plt.plot(X_CT,exp_fit(X_CT, fitParams_CT[0][0], fitParams_CT[0][1], fitParams_CT[0][2]),'r',ms=10,linewidth=2.0,label='C>T')
##    plt.plot(X_GA,exp_fit(X_GA, fitParams_GA[0][0], fitParams_GA[0][1], fitParams_GA[0][2]),'b',ms=10,linewidth=2.0,label='G>A')
##    plt.legend()
##
##
##    # save plot to a file
##    plt.savefig(plotfile, bbox_inches=0, dpi=600)
##    plt.close()
##    print '\nPlotted MapDamage curve fit to '+plotfile

    #Make a reference list to determine how much to adjust a particular quality score given it's read position (maxed at 300 here)
    CT_decay=[]
    GA_decay=[]

    for g in range(300):
        CT_decay.append(exp_fit(g+1,CTa,CTb,CTc)-theta)
        GA_decay.append(exp_fit(g+1,GAa,GAb,GAc)-theta)


###open up reference file
ref=pysam.FastaFile(ref_file)


###set up various look up dictionaries
all_dic={}
all_dic['A']=0
all_dic['C']=1
all_dic['G']=2
all_dic['T']=3
all_dic[0]='A'
all_dic[1]='C'
all_dic[2]='G'
all_dic[3]='T'
all_dic[-9]='.'

GL_dic={}

GL_dic['AA']=0
GL_dic['AC']=1
GL_dic['AG']=2
GL_dic['AT']=3
GL_dic['CC']=4
GL_dic['CG']=5
GL_dic['CT']=6
GL_dic['GG']=7
GL_dic['GT']=8
GL_dic['TT']=9

GL_dic[0]='AA'
GL_dic[1]='AC'
GL_dic[2]='AG'
GL_dic[3]='AT'
GL_dic[4]='CC'
GL_dic[5]='CG'
GL_dic[6]='CT'
GL_dic[7]='GG'
GL_dic[8]='GT'
GL_dic[9]='TT'
GL_dic[-9]='./.'

homs=[0,4,7,9]

alt_dic={}
alt_dic[0]=[0,0]
alt_dic[1]=[0,1]
alt_dic[2]=[0,2]
alt_dic[3]=[0,3]
alt_dic[4]=[1,1]
alt_dic[5]=[1,2]
alt_dic[6]=[1,3]
alt_dic[7]=[2,2]
alt_dic[8]=[2,3]
alt_dic[9]=[3,3]

tri_dic={}
tri_dic[1]='A,C'
tri_dic[2]='A,G'
tri_dic[3]='A,T'
tri_dic[5]='C,G'
tri_dic[6]='C,T'
tri_dic[8]='G,T'

rev_com={}
rev_com['A']='T'
rev_com['C']='G'
rev_com['G']='C'
rev_com['T']='A'


LL0_map={}
LL0_map[0]=[0]
LL0_map[1]=[4]
LL0_map[2]=[7]
LL0_map[3]=[9]
            
LL1_map={}
LL1_map[0]=[1,2,3]
LL1_map[1]=[1,5,6]
LL1_map[2]=[2,5,8]
LL1_map[3]=[3,6,8]
            
LL2_map={}
LL2_map[0]=[4,5,6,7,8,9]
LL2_map[1]=[0,2,3,7,8,9]
LL2_map[2]=[0,1,3,4,6,9]
LL2_map[3]=[0,1,2,4,5,7]


###make a list of regions to be interogated
file = open(filenameinB)
data=file.read()
data=string.split(data,'\n')
file.close()

if data[-1]=='':
    del(data[-1])

SNPdir={}
SNPll={}
SNPlist=[]

count=0

for g in range(len(data)):
    k=string.split(data[g])
    key=k[0]+'_'+k[1]+'_'+k[2]
    chromo=k[0]
    start=int(k[1])
    end=int(k[2])
    all1=k[3]
    all2=k[4]
    SNPlist.append(key)
    SNPdir[key]=k            




###open up bam file (must be indexed)
samfile = pysam.AlignmentFile(filenamein, "rb")
samp_name=samfile.header['RG'][0]['SM']
count=0


###set up output files
fileout=open(filenameout1,'w')
out='chr\tpos\tref\tgeno1\tgeno2\tgeno3\tall1\tall2\tDP_all1\tDP_all2\n'
fileout.write(out)



###Start running through regions in the bed
for gg in range(len(SNPlist)):
    print 'Calculating genotype likelihoods for locus '+str(gg+1)+' : '+SNPlist[gg]
    chromo=SNPdir[SNPlist[gg]][0]
    pos_start=int(SNPdir[SNPlist[gg]][1])
    pos_end=int(SNPdir[SNPlist[gg]][2])

    ###set up arrays to store results for a given bed entry (probably not to efficient for SNP data, more suited for long regions)
    seq_len=pos_end-pos_start
    GLs=np.zeros((seq_len,10),dtype='float32')  ##Genotype likelihoods
    RDs=np.zeros((seq_len,4),dtype='int32')     ##Read depth for each of the four bases
    REFs=np.zeros((seq_len),dtype='int32')      ##Reference allele index
    POS=np.zeros((seq_len),dtype='int32')       ##1-based position
    
    ###make a list of positions
    count=0
    for ggg in range(pos_start,pos_end):
        POS[count]=ggg+1
        count+=1

    ###prepopulate arrays with -9s
    REFs.fill(0)
    GLs.fill(1.0/10.0)

    ###Go base by base through each region in the bam file
    for pileupcolumn in samfile.pileup(chromo,pos_start,pos_end,truncate=True,stepper='all'):
        nucl=ref.fetch(chromo,pileupcolumn.pos,pileupcolumn.pos+1)  ##grab reference sequence for region
        nucl=nucl.upper()
        var_list=[]  ##a list to store read bases
        index=pileupcolumn.pos-pos_start  ##tells us what position we are in for the arrays given the nucleotide position
        REFs[index]=all_dic[nucl]  ##recrod the reference allele

        ###go read by read for each position
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:  ##ensure not an indel or duplicate
                if (pileupread.alignment.mapping_quality>=MQ) and (ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ):  ###Add base calls meeting the MQ and BQ filters
                    var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33,CT_decay[pileupread.query_position],GA_decay[len(pileupread.alignment.query_sequence)-1-pileupread.query_position]])

        ###rescale qualities that are greater than 40 to a max of 40.
        for ggg in range(len(var_list)):
            if var_list[ggg][1]>40:
                var_list[ggg][1]=40

        ###record read depth for each basetype
        if len(var_list)>0:
            all_list=list(zip(*var_list)[0])
            RDs[index]=[all_list.count('A'),all_list.count('C'),all_list.count('G'),all_list.count('T')]
        

        ###if minimum read depth is met, try calling the genotype 
        if len(var_list)>=min_RD:            
            GLs[index]=geno_caller_10GT_aDNA(var_list)  ##ancient DNA aware calculation of genotype likelihoods

    ###start writing calls to  vcf file                   
    for ggg in range(len(GLs)):
        out=chromo+'\t'+str(POS[ggg])+'\t'+all_dic[REFs[ggg]]+'\t'
        SNPdir[SNPlist[gg]]
        k=SNPdir[SNPlist[gg]]
        alls1=k[3]
        alls2=k[4]
        try:
            GL_val=np.array([GLs[ggg][GL_dic[alls1+alls1]],GLs[ggg][GL_dic[alls1+alls2]],GLs[ggg][GL_dic[alls2+alls2]]])
        except:
            GL_val=np.array([GLs[ggg][GL_dic[alls1+alls1]],GLs[ggg][GL_dic[alls2+alls1]],GLs[ggg][GL_dic[alls2+alls2]]])
        GL_val_sc=GL_val-np.max(GL_val)       
        scale_p=(10**GL_val_sc)/np.sum(10**GL_val_sc)
        out=out+str(round(scale_p[0],3))+'\t'+str(round(scale_p[1],3))+'\t'+str(round(scale_p[2],3))+'\t'+alls1+'\t'+alls2+'\t'+str(RDs[ggg][all_dic[alls1]])+'\t'+str(RDs[ggg][all_dic[alls2]])+'\n'
        fileout.write(out) ##write to emit all

fileout.close()
        


