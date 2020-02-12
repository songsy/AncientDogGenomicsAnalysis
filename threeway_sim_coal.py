##import time
import subprocess
from subprocess import Popen
import os
##import dircache
import string
from string import join
import math
##from time import strftime
import random
from random import randint
from random import uniform
from random import gauss
from random import gammavariate
from random import betavariate
from math import sqrt
from sys import argv
import numpy as np
import datetime
import ms

##########################################
def prob_coal(theta,t):
    prob=1-math.exp(-(2.0/theta)*t)
    return prob

#Analytical solution
theta1=4E-4
theta2=4E-4
tau1=2.5E-5
tau2=4.5E-5
tau3=7E-5
alpha=0.9
beta=1-alpha


prob_AB_a=alpha* (prob_coal(theta1,tau2-tau1) + ((1-prob_coal(theta1,tau2-tau1))/3))
prob_AC_a=alpha* ((1-prob_coal(theta1,tau2-tau1))/3)

prob_AB_b=beta* ((1-prob_coal(theta2,tau3-tau2))/3)
prob_AC_b=beta* (prob_coal(theta2,tau3-tau2) + ((1-prob_coal(theta2,tau3-tau2))/3))

exp_rat=(prob_AB_a+prob_AB_b)/(prob_AC_a+prob_AC_b)



##########################################
#Set up simulation characteristics common to all runs

nbseq=100000 #16434      #number of loci ## LRB length of bed file from GPhocs
seq_len=1000   #length of loci   ## LRB Number of bp in each segment.
n=[1,1,1,1]             #number of chromosomes, ## Eur,HXH, Ind, Asia
job=1        
##########################################
#Start individual simulation

para_out=[]     #list to store parameters and values

##########################################
#Parameter set up


####################################################################################    
####################################################################################
####################################################################################
####################################################################################
#Simulate invidual loci. Paste in specific ms set arguments here

##########################################
#Set up global arg (though 'Theta' will change)
ms_args = ['ms',str(sum(n)),'1','-t',str(theta1*seq_len),'-I',str(len(n))] ## ms TotalSampleSize 1 -t Theta -I NumberPops

#-----------Population structure of island model-----------#
for x in range(len(n)):
    ms_args.append(str(n[x])) ## Append samples sized per population

#-----------Migration pulse (ASI>HXH)--------------#            
ms_args.extend(['-es',str((tau1-(tau1*1E-12))/theta1),'2',str(alpha)]) ##


###-----------First population split1 (EUR/HXH)--------------#            
##ms_args.extend(['-ej',str(T1/(4*float(N1))),'2','1']) ##
#-----------First population split1 (EUR/HXH)--------------#            
ms_args.extend(['-ej',str((tau1)/theta1),'2','1']) ##
ms_args.extend(['-ej',str((tau1)/theta1),'5','4']) ##


#-----------Second population split2 (EUR/IND)--------------#            
ms_args.extend(['-ej',str((tau2)/theta1),'3','1']) ##

#-----------Population size change split2--------------#            
ms_args.extend(['-en',str((tau2)/theta1),'1',str(theta2/theta1)]) ## 4 was written in gibbon

#-----------Second population split3 (EUR/ASI)--------------#            
ms_args.extend(['-ej',str((tau3)/theta1),'4','1']) ##


reg_use=0

AB=0
AC=0

while reg_use<nbseq:

    #RUN SIMULATION
    ms.ms_main(len(ms_args), ms_args)
    ms.setSeed(randint(1,100000)*job,randint(1,100000)*job,randint(1,100000)*job)   
    ms.generateSamples()

    #Process msoutput           
    nbss=ms.getNumSites()
    
#	print nbss
    
    pos=[]           
    for x in range(nbss):
        pos.append(ms.getPosition(x,seq_len))           


    alleles=[]
    for x in range(nbss):
        loc=[]
        for m in range(sum(n)):
            loc.append(ms.getSite(m,x))
        alleles.append(loc)

    alleles_use=[]
    for g in range(len(alleles)):
        if alleles[g][0]=='1':
            if (alleles[g][1]=='1') and (alleles[g][2]=='0'):
                alleles_use.append(alleles[g])
            elif  (alleles[g][1]=='0') and (alleles[g][2]=='1'):
                alleles_use.append(alleles[g])

    if len(alleles_use)>0:
        random.shuffle(alleles_use)
        if string.join(alleles_use[0][:3],'')=='110':
            AB+=1
        elif string.join(alleles_use[0][:3],'')=='101':
            AC+=1

    reg_use=reg_use+1

print 'analytical\tsimulated'
print str(exp_rat)+'\t'+str(float(AB)/float(AC))
