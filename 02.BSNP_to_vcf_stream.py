#!/usr/bin/env python
# python 02.BSNP_to_vcf.py
# Shiya Song
# 27th October 2015

import math
import sys
import os
import gzip

name = sys.argv[1]
print '##fileformat=VCFv4.1\n##samtoolsVersion=0.1.17 (r973:277)\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">'
print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
print '##FORMAT=<ID=PL,Number=.,Type=Integer,Description="List of Phred-scaled genotype likelihoods, number of values is (#ALT+1)*(#ALT+2)/2">'
print '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' %(name)

field={"AA":9,"AC":10,"AG":11,"AT":12,"CC":13,"CG":14,"CT":15,"GG":16,"GT":17,"TT":18}
annotation={'R':'AG','Y':'CT','K':'GT','M':'AC','S':'CG','W':'AT','A':'AA','T':'TT','C':'CC','G':'GG','N':'NN'}
for line in sys.stdin:
	if line[0]=="#":
		continue
	else:
		line = line.strip().split("\t")
		chr = line[0]
		pos = int(line[1])
		ref = line[2]
		alt = line[4]
		DP = line[30]
		qual = line[5]
		if ref==alt:
			LL = 1-float(line[field[ref+alt]])
			try:
				GQ = int(round(-10*math.log(LL,10)))
			except ValueError:
				GQ= 1000
#				print >>f_out,"%s\t%i\t.\t%s\t.\t%s\t.\tDP=%s\tPL\t0" %(chr,pos,ref,qual,DP)
			print "%s\t%i\t.\t%s\t.\t%s\t.\tDP=%s\tGT:PL:GQ\t%s:%s:%s" %(chr,pos,ref,qual,DP,"0/0","0",GQ)
		else:
			if alt!="N" and ref!="N":
#					LL = 1-float(line[field[annotation[alt]]])
#					qual = round(-10*math.log(LL,10))
				if annotation[alt][0]!=annotation[alt][1]:
					if annotation[alt][0]==ref or annotation[alt][1]==ref:
						tag1="0/1"
						L1=float(line[field[annotation[ref]]])
						L2=float(line[field[annotation[alt]]])
						if annotation[alt][0]==ref:
							another = annotation[alt][1]
						else:
							another = annotation[alt][0]
						L3=float(line[field[annotation[another]]])
						GQ=1-L2
						try:
							tag3=int(round(-10*math.log(GQ,10)))
						except ValueError:
#								print pos,L1,L2,L3,GQ
							tag3=1000
						try:
							p1=-10*math.log(L1/max(L1,L2,L3),10)
						except ValueError:
							p1=1000
						try:
							p2=-10*math.log(L2/max(L1,L2,L3),10)
						except ValueError:
							p2=1000
						try:
							p3=-10*math.log(L3/max(L1,L2,L3),10)
						except ValueError:
							p3=1000
						tag2=str(int(round(p1)))+","+str(int(round(p2)))+","+str(int(round(p3)))
						print "%s\t%i\t.\t%s\t%s\t%s\t.\tDP=%s\tGT:PL:GQ\t%s:%s:%i" %(chr,pos,ref,another,qual,DP,tag1,tag2,tag3)
					else:
						tag1="1/2"
						LL=1-float(line[field[annotation[alt]]])
						try:
							GQ = int(round(-10*math.log(LL,10)))
						except ValueError:
							GQ= 1000
						print "%s\t%i\t.\t%s\t%s,%s\t%s\t.\tDP=%s\tGT:PL:GQ\t%s:%s:%i" %(chr,pos,ref,annotation[alt][0],annotation[alt][1],qual,DP,"1/2","0",GQ)
				else:
					if annotation[ref][0]!=annotation[ref][1]:
						ref = 'N'
						print "%s\t%i\t.\t%s\t%s\t%s\t.\tDP=%s\tGT:PL:GQ\t0:0:0" %(chr,pos,ref,alt,qual,DP)
						continue
					tag1="1/1"
					L1=float(line[field[annotation[ref]]])
					if ref<alt:
						allele=ref+alt
					else:
						allele=alt+ref
					L2=float(line[field[allele]])
					L3=float(line[field[annotation[alt]]])
#						print pos,L1,L2,L3
					GQ=1-L3
					try:
						tag3=int(round(-10*math.log(GQ,10)))
					except ValueError:
#							print pos,L1,L2,L3,GQ
						tag3=1000
					try:
						p1=-10*math.log(L1/max(L1,L2,L3),10)
					except ValueError:
						p1=1000
					try:
						p2=-10*math.log(L2/max(L1,L2,L3),10)
					except ValueError:
						p2=1000
					try:
						p3=-10*math.log(L3/max(L1,L2,L3),10)
					except ValueError:
						p3=1000
					tag2=str(int(round(p1)))+","+str(int(round(p2)))+","+str(int(round(p3)))
					print "%s\t%i\t.\t%s\t%s\t%s\t.\tDP=%s\tGT:PL:GQ\t%s:%s:%i" %(chr,pos,ref,alt,qual,DP,tag1,tag2,tag3)
			else:
				print "%s\t%i\t.\t%s\t%s\t%s\t.\tDP=%s\tGT:PL:GQ\t0:0:0" %(chr,pos,ref,alt,qual,DP)
#cmds ='rm '+m
#os.popen(cmds)