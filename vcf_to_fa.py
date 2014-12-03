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
import gzip
from numpy.random import binomial
import argparse

annotation={'AG':'R','GA':'R','TC':'Y','CT':'Y','GT':'K','TG':'K','AC':'M','CA':'M','GC':'S','CG':'S','AT':'W','TA':'W','XA':'N','XG':'N','XC':'N','XT':'N','AX':'N','GX':'N','CX':'N','TX':'N'}

def read_trf_file(trf_file):
	trf = {}
	print trf_file[-2:]
	if trf_file[-2:]=='gz':
		f=gzip.open(trf_file,"r")
	else:
		f=open(trf_file,"r")
	for line in f:
		line=line.strip().split("\t")
		chr = line[0]
		if chr not in trf.keys():
			trf[chr]=[]
#		if chr != "chr1":
#			continue
		trf[chr].append([line[0],int(line[1])+1,int(line[2])])
	return trf

def read_locus_file(file):              # locus interval is bed format, zero based, transform into 1-based
	locus={}
	f=open(file,"r")
	for line in f:
		line = line.strip().split("\t")
#		if line[0] != "chr1":
#			break
		chr = line[0]
		if chr=="chr2a":
			chr = "chr2A"
		if chr=="chr2b":
			chr = "chr2B"
		info = "%s:%s-%s" %(line[0],line[1],line[2])
		col = line[3].split(":")
		chr = col[0]
		pos1 = col[1].split("-")[0]
		pos2 = col[1].split("-")[1]
		if chr not in locus.keys():
			locus[chr]=[]
		locus[chr].append([chr,int(pos1)+1,int(pos2),info])   # first human then chimp
#	print len(locus)
	return locus	
					
def parse_vcf(vcf_file,mask,locus,trf,fasta_file,hg19_seq):
	random = 0
	f = gzip.open(vcf_file,"r")
	total_length=0
	het = 0
	i=0
#	print locus[i]
	hg19_sequence = hg19_seq['%s:%i-%i' %(locus[i][0],locus[i][1]-1,locus[i][2])]
	string = []
	string1 = []
	string2 = []
	mask_index = 0
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
#		print locus[i][1]+j,mask[mask_index][1],mask[mask_index][2]
		if (locus[i][1]+j)>=mask[mask_index][1] and (locus[i][1]+j)<=mask[mask_index][2]:
#			print j,mask_index,string[j]
			string[j]=2
			string1[j]=2
			string2[j]=2
		elif (locus[i][1]+j)>mask[mask_index][2]:
			while True:
				mask_index +=1
				if (locus[i][1]+j)<=mask[mask_index][2]:
					break
			if (locus[i][1]+j)>=mask[mask_index][1] and (locus[i][1]+j)<=mask[mask_index][2]:
				print j,mask_index,string[j]
				string[j]=2
				string1[j]=2
				string2[j]=2
	for line in f:
		if line[0]=='#':
			continue
		line = line.strip().split("\t")
		info = line[9]
		chr = line[0]
		pos = int(line[1])
		try:
			DP = int(info.split(":")[2])
		except IndexError:
			DP = 10
		if i>=len(locus):
			break
		if pos>=locus[i][1] and pos<=locus[i][2]:
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
					if line[4]==".":					# reference allele
						string[pos-locus[i][1]]=line[3]
						string1[pos-locus[i][1]]=line[3]
						string2[pos-locus[i][1]]=line[3]
					else:	# non reference allele
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
							print 'wrong',line
							string[pos-locus[i][1]]="N"
							string1[pos-locus[i][1]]="N"
							string2[pos-locus[i][1]]="N"
							continue
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
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--vcf_dir",default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',dest='vcf_dir',help="vcf dir")
	args = parser.parse_args()		
	chromOrder=get_chromOrder("human")
	trf_file = "/home/jmkidd/kidd-lab-scratch/shiya-projects/G-Phocs/hg19.simplerepeats.bed"
	locus_file = "locus-chimp.bed"
	locus=read_locus_file(locus_file)
	trf = read_trf_file(trf_file)
	outfile=open("heterozygosity_%s" %(args.sample),"w")
	fasta_file = [open("locus-%s.fa" %(args.sample),"w"),open("locus-%s-hap1.fa" %(args.sample),"w"),open("locus-%s-hap2.fa" %(args.sample),"w")]
	tot_het = 0
	tot_l = 0
	tot_random = 0
	for chr in range(1,23):
		hg19_seq = read_seq('locus-hg19.fa','chr'+str(chr))
		locus['chr'+str(chr)].sort( key=lambda w: int(w[1]))    # locus are in hg19
#		print locus['chr'+str(chr)]
#		vcf_file = "%s%s/all_sites/%s.chr%s.moleculo.phased.vcf.gz" %(args.vcf_dir,args.sample,args.sample,chr)
#		mask_file = "%s%s/all_sites/%s_chr%s.mask.bed.gz" %(args.vcf_dir,args.sample,args.sample,chr)
#		vcf_file = "%s%s/gVCF_calls/%s.chr%s.prism.v2.phased.vcf.gz" %(args.vcf_dir,args.sample,args.sample,chr)
		vcf_file = "%s%s/gVCF_calls/%s.%s.prism.phased.vcf.gz" %(args.vcf_dir,args.sample,args.sample,chr)
		mask_file = "%s%s/gVCF_calls/%s_%s.mask.bed.gz" %(args.vcf_dir,args.sample,args.sample,chr)	
#		vcf_file = "/share/jmkidd/songsy/complete-genomics/MKK/%s.chr%s.phased.vcf.gz" %(args.sample,chr)
#		mask_file = "/share/jmkidd/songsy/complete-genomics/MKK/%s_chr%s.mask.bed.gz" %(args.sample,chr)
		mask=read_trf_file(mask_file)
		het,total_length,random=parse_vcf(vcf_file,mask[str(chr)],locus['chr'+str(chr)],trf['chr'+str(chr)],fasta_file,hg19_seq)      # sys.argv[3] True if I made them after BSNP output, false if otherwise
		print >>outfile,chr,het,random,total_length,float(het)/total_length
		tot_het+=het
		tot_l+=total_length
		tot_random+=random
	print >>outfile,'tot',tot_het,tot_random,tot_l,float(tot_het)/tot_l
		

			
		
	

 