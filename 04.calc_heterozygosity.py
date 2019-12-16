#!/usr/bin/env python
# python 03.parse_interval_vcf.py
# Shiya Song
# 3rd June 2013
# Parse vcf file, mark indel(around 3bp) as N, mark those with DP<5 as N, mark simple repeats as N, record clustered SNP

import re
import genutils
import numpy as np
import math
import sys
from NGS_utils import *
from itertools import chain
import gzip
from numpy.random import binomial
import argparse

annotation={'AG':'R','GA':'R','TC':'Y','CT':'Y','GT':'K','TG':'K','AC':'M','CA':'M','GC':'S','CG':'S','AT':'W','TA':'W','XA':'N','XG':'N','XC':'N','XT':'N','AX':'N','GX':'N','CX':'N','TX':'N'}

def get_chrom_list():
	chrom_list=[]
	chromOrder={}
	f=open('/home/jmkidd/kidd-lab/genomes/canFam3.1/chromOrder.txt','r')
	for i,line in enumerate(f):
		chrom_fa = line.rstrip()
		chrom=chrom_fa.replace('.fa','')
		chrom_list.append(chrom)
		if chrom_fa=='chrM.fa':
			break
		chromOrder[chrom] = i
	return chrom_list,chromOrder

def read_trf_file(trf_file):
	trf = {}
	if trf_file[-2:]=='gz':
		f=gzip.open(trf_file,"r")
	else:
		f=open(trf_file,"r")
	for line in f:
		line=line.strip().split("\t")
		chr = line[0]
		if chr not in trf.keys():
			trf[chr]=[]
		trf[chr].append([line[0],int(line[1])+1,int(line[2])])
	return trf

def read_locus_file(file):              # locus interval is bed format, zero based, transform into 1-based
	locus={}
	f=open(file,"r")
	for line in f:
		line = line.strip().split("\t")
		chr = line[0]
		pos1 = line[1]
		pos2 = line[2]
		if chr not in locus.keys():
			locus[chr]=[]
		locus[chr].append([chr,int(pos1)+1,int(pos2)]) 
	return locus	
					
def calc_heterozygosity_window_based(vcf_file,trf,rmsk,chr,fout):
	if chr in trf.keys():
		trf_mask=set(chain(*(xrange(start, end+1) for chrom,start,end in trf[chr])))
	else:
		trf_mask=[]
	if chr in rmsk.keys():
		rmsk_mask=set(chain(*(xrange(start, end+1) for chrom,start,end in rmsk[chr])))
	else:
		rmsk_mask=[]
	pos_index=0
	first = True
	f = gzip.open(vcf_file,"r")
	prev_genotype="AA"
	prev_CpG=False
	prev_het=False
	for line in f:
		if line[0]=='#':
			continue
		line = line.strip().split("\t")
		info = line[7]
		chr = line[0]
		pos = int(line[1])
		m1 = re.search('DP=(\d+)',info)
		genotype = line[9].split(":")[0]
		if pos_index%200000==0:
			if first is True:
				first=False
			else:
#				window.append([chr,pos1,pos2,het,L,float(het)/L])
				info=[chr,pos1,pos2,het1,L1,float(het1)/L1,het2,L2,float(het2)/L2]
				print >>fout,"\t".join(map(str,info))
#				print "\t".join(map(str,info))
			pos1=pos
			het1,L1,het2,L2=(0,0,0,0)
			rev_genotype="AA"
			prev_CpG=False
			prev_het=False
		pos_index+=1
		pos2=pos
		if int(m1.group(1))<5:   # DP<5
			prev_genotype="AA"
			prev_CpG=False
			prev_het=False
			continue
		if pos in trf_mask:
			prev_genotype="AA"
			prev_CpG=False
			prev_het=False
			continue
		if pos in rmsk_mask:
			prev_genotype="AA"
			prev_CpG=False
			prev_het=False
			continue
		L1+=1
		L2+=1
		if genotype=="0/1":
			het1+=1
			cur_genotype=line[3]+line[4]
		elif genotype=="1/1":
			cur_genotype=line[4]+line[4]
		elif genotype=="0":
			cur_genotype=line[3]+line[3]
		if "C" in prev_genotype and "G" in cur_genotype:
			if prev_CpG:
				L2-=1
			else:
				L2-=2
			prev_CpG=True
			if prev_het:
				het2-=1
			prev_het=False
		else:
			if genotype=="0/1":
				het2+=1
				prev_het=True
			else:
				prev_het=False
		prev_genotype=cur_genotype

def parse_vcf(vcf_file,mask,locus,trf,fasta_file,hg19_seq):
	random = 0
	f = gzip.open(vcf_file,"r")
	total_length=0
	het = 0
	i=0
	string = []
	mask_index = 0
	for j in range(locus[i][2]-locus[i][1]+1):
		string.append(0)
	for j in range(len(trf)):
		if trf[j][1]>= locus[i][1] and trf[j][1]<=locus[i][2]:
			if trf[j][2]<=locus[i][2]:
				for k in range(trf[j][1]-locus[i][1],trf[j][2]-locus[i][1]+1):
					string[k]=1
			else:
				for k in range(trf[j][1]-locus[i][1],locus[i][2]-locus[i][1]+1):
					string[k]=1
		elif trf[j][2]>= locus[i][1] and trf[j][2]<=locus[i][2]:
			if trf[j][1]<=locus[i][1]:
				for k in range(0,trf[j][2]-locus[i][1]+1):
					string[k]=1
			else:
				for k in range(trf[j][1]-locus[i][1],trf[j][2]-locus[i][1]+1):
					string[k]=1
		elif trf[j][1]<=locus[i][1] and trf[j][2]>=locus[i][2]:
			for k in range(0,locus[i][2]-locus[i][1]+1):
				string[k]=1
		elif trf[j][1]>locus[i][2]:
			break

	for line in f:
		if line[0]=='#':
			continue
		line = line.strip().split("\t")
		info = line[7]
		chr = line[0]
		pos = int(line[1])
		m1 = re.search('DP=(\d+)',info)
		if i>=len(locus):
			break
		if pos>=locus[i][1] and pos<=locus[i][2]:
			if string[pos-locus[i][1]]==1:
				string[pos-locus[i][1]]="N"
			else:
				if DP<5:
					string[pos-locus[i][1]]="N"
				elif string[pos-locus[i][1]]!="N":
					if line[4]==".":					# reference allele
						string[pos-locus[i][1]]=line[3]
					else:	# non reference allele
						genotype = info.split(":")[0]
						allele = line[4]
						if genotype=="0/1":
							string[pos-locus[i][1]]=annotation[allele+line[3]]
							het +=1
						elif genotype=="1/1":
							string[pos-locus[i][1]]=allele
						else:
							print 'wrong',line
		elif pos>locus[i][2]:
			seq=""
			seq1=""
			seq2=""
			length = locus[i][2]-locus[i][1]+1
			for m in range(locus[i][1],locus[i][2]+1):
				if string[m-locus[i][1]]!=0 and string[m-locus[i][1]]!=1 and string[m-locus[i][1]]!=2:
					seq += string[m-locus[i][1]]
					seq1 +=string1[m-locus[i][1]]
					seq2 +=string2[m-locus[i][1]]
				elif string[m-locus[i][1]]==2:
					seq += hg19_sequence[m-locus[i][1]]
					seq1 += hg19_sequence[m-locus[i][1]]
					seq2 += hg19_sequence[m-locus[i][1]]
				else:
					seq += "N"
					seq1 += "N"
					seq2 += "N"
					length -=1
			print >>fasta_file[0],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
			print >>fasta_file[1],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
			print >>fasta_file[2],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
			print >>fasta_file[0],seq
			print >>fasta_file[1],seq1
			print >>fasta_file[2],seq2
			total_length += length
			while True:
				i +=1
				if i>=len(locus):
					break
#				print locus[i]
				hg19_sequence = hg19_seq['%s:%i-%i' %(locus[i][0],locus[i][1]-1,locus[i][2])]
				string = []
				string1 = []
				string2 = []
				for j in range(locus[i][2]-locus[i][1]+1):
					string.append(0)
					string1.append(0)
					string2.append(0)
				for j in range(len(trf)):
					if trf[j][1]>= locus[i][1] and trf[j][1]<=locus[i][2]:
						if trf[j][2]<=locus[i][2]:
							for k in range(trf[j][1]-locus[i][1],trf[j][2]-locus[i][1]+1):
								string[k]=1
								string1[k]=1
								string2[k]=1
						else:
							for k in range(trf[j][1]-locus[i][1],locus[i][2]-locus[i][1]+1):
								string[k]=1
								string1[k]=1
								string2[k]=1
					elif trf[j][2]>= locus[i][1] and trf[j][2]<=locus[i][2]:
						if trf[j][1]<=locus[i][1]:
							for k in range(0,trf[j][2]-locus[i][1]+1):
								string[k]=1
								string1[k]=1
								string2[k]=1
						else:
							for k in range(trf[j][1]-locus[i][1],trf[j][2]-locus[i][1]+1):
								string[k]=1
								string1[k]=1
								string2[k]=1
					elif trf[j][1]<=locus[i][1] and trf[j][2]>=locus[i][2]:
						for k in range(0,locus[i][2]-locus[i][1]+1):
							string[k]=1
							string1[k]=1
							string2[k]=1
					elif trf[j][1]>locus[i][2]:
						break
				while True:
					if locus[i][1]>mask[mask_index][2]:
						mask_index+=1
					else:
						break
				for j in range(locus[i][2]-locus[i][1]+1):
					if (locus[i][1]+j)>=mask[mask_index][1] and (locus[i][1]+j)<=mask[mask_index][2]:
						string[j]=2
						string1[j]=2
						string2[j]=2
					elif (locus[i][1]+j)>mask[mask_index][2]:
						while True:
							mask_index +=1
							if (locus[i][1]+j)<=mask[mask_index][2]:
								break
						if (locus[i][1]+j)>=mask[mask_index][1] and (locus[i][1]+j)<=mask[mask_index][2]:
							string[j]=2
							string1[j]=2
							string2[j]=2
				if pos<locus[i][1]:
					break
				elif pos>=locus[i][1] and pos<=locus[i][2]:
					if string[pos-locus[i][1]]==1:
						string[pos-locus[i][1]]="N"
						string1[pos-locus[i][1]]="N"
						string2[pos-locus[i][1]]="N"
					else:
						if DP<5:
							string[pos-locus[i][1]]="N"
							string1[pos-locus[i][1]]="N"
							string2[pos-locus[i][1]]="N"
						elif string[pos-locus[i][1]]!="N":
							genotype = info.split(":")[0]
							allele = line[4]
							if genotype=="0|1":
								string[pos-locus[i][1]]=annotation[allele+line[3]]
								string1[pos-locus[i][1]]=line[3]
								string2[pos-locus[i][1]]=allele
								het +=1
							elif genotype=="1|0":
								string[pos-locus[i][1]]=annotation[allele+line[3]]
								string2[pos-locus[i][1]]=line[3]
								string1[pos-locus[i][1]]=allele
								het +=1
							elif genotype=="1|1" or genotype=="1/1":
								string[pos-locus[i][1]]=allele
								string1[pos-locus[i][1]]=allele
								string2[pos-locus[i][1]]=allele
							elif genotype=="0/1" or genotype=="1/0":
								het +=1
								random+=1
								string[pos-locus[i][1]]=annotation[allele+line[3]]
								'''
								indicator = binomial(1,0.5)
								if indicator == 0:
									string1[pos-locus[i][1]]=line[3]
									string2[pos-locus[i][1]]=allele
								else:
									string2[pos-locus[i][1]]=line[3]
									string1[pos-locus[i][1]]=allele
								'''
								string1[pos-locus[i][1]]="N"
								string2[pos-locus[i][1]]="N"
							elif genotype=="0|0":
								string[pos-locus[i][1]]=line[3]
								string1[pos-locus[i][1]]=line[3]
								string2[pos-locus[i][1]]=line[3]
							else:
								string[pos-locus[i][1]]="N"
								string1[pos-locus[i][1]]="N"
								string2[pos-locus[i][1]]="N"
								print 'wrong2',line
					break
				elif pos>locus[i][2]:
					seq=""
					seq1=""
					seq2=""
					length = locus[i][2]-locus[i][1]+1
					for m in range(locus[i][1],locus[i][2]+1):
						if string[m-locus[i][1]]!=0 and string[m-locus[i][1]]!=1 and string[m-locus[i][1]]!=2:
							seq += string[m-locus[i][1]]
							seq1 +=string1[m-locus[i][1]]
							seq2 +=string2[m-locus[i][1]]
						elif string[m-locus[i][1]]==2:
							seq += hg19_sequence[m-locus[i][1]]
							seq1 += hg19_sequence[m-locus[i][1]]
							seq2 += hg19_sequence[m-locus[i][1]]
						else:
							seq += "N"
							seq1 += "N"
							seq2 += "N"
							length -=1
					print >>fasta_file[0],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
					print >>fasta_file[1],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
					print >>fasta_file[2],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
					print >>fasta_file[0],seq
					print >>fasta_file[1],seq1
					print >>fasta_file[2],seq2
					total_length += length
	if i<len(locus):
		seq=""
		seq1=""
		seq2=""
		length = locus[i][2]-locus[i][1]+1
		for m in range(locus[i][1],locus[i][2]+1):
			if string[m-locus[i][1]]!=0 and string[m-locus[i][1]]!=1 and string[m-locus[i][1]]!=2:
				seq += string[m-locus[i][1]]
				seq1 +=string1[m-locus[i][1]]
				seq2 +=string2[m-locus[i][1]]
			elif string[m-locus[i][1]]==2:
				seq += hg19_sequence[m-locus[i][1]]
				seq1 += hg19_sequence[m-locus[i][1]]
				seq2 += hg19_sequence[m-locus[i][1]]
			else:
				seq += "N"
				seq1 += "N"
				seq2 += "N"
				length -=1
		print >>fasta_file[0],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
		print >>fasta_file[1],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
		print >>fasta_file[2],">%s.%i-%i %s %i %i" %(locus[i][0],locus[i][1],locus[i][2],locus[i][3],len(seq),len(seq))
		print >>fasta_file[0],seq
		print >>fasta_file[1],seq1
		print >>fasta_file[2],seq2
		total_length += length
	return het,total_length,random

def read_seq(fasta,chr):
	f=open(fasta,'r')
	find = False
	seq = {}
	for line in f:
		if line[0]=='>':
			chrom = line.strip().split(':')[0][1:]
			if chrom==chr:
				find = True
				tag=line.strip()[1:]
		elif find is True:
			seq[tag]=line.strip()
			find is False
	return seq

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='vcf to fa')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--vcf_dir",default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/BSNP',dest='vcf_dir',help="vcf dir")
	parser.add_argument("--chr",dest='chr',help="chr")
	args = parser.parse_args()		
	chrom_list,chromOrder=get_chrom_list()
	trf_file = "/home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1_trf.bed"
	rmsk_file = "/home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1_rmsk.bed"
	locus_file = "/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/neutralLoci-geneFlank10k-sep30k-3.1.bed"
	locus=read_locus_file(locus_file)
	trf = read_trf_file(trf_file)
	rmsk = read_trf_file(rmsk_file)
	print 'finish reading files'
	if True:
		print args.chr
		vcf_file = "%s/%s_%s_BSNP.vcf.gz" %(args.vcf_dir,args.sample,args.chr)
		fout=open("%s/%s_%s_heterzygosity.txt" %(args.vcf_dir,args.sample,args.chr),"w")
		calc_heterozygosity_window_based(vcf_file,trf,rmsk,args.chr,fout)
		

			
		
	

 