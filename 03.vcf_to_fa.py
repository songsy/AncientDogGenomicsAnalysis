#!/usr/bin/env python
# python 03.parse_interval_vcf.py
# Shiya Song
# 3rd June 2013
# Parse vcf file, mark indel(around 3bp) as N, mark those with DP<5 as N, mark simple repeats as N, record clustered SNP

import re
import numpy as np
import math
import sys
import gzip
from numpy.random import binomial
import argparse

annotation={'0':'N','AG':'R','GA':'R','TC':'Y','CT':'Y','GT':'K','TG':'K','AC':'M','CA':'M','GC':'S','CG':'S','AT':'W','TA':'W','XA':'N','XG':'N','XC':'N','XT':'N','AX':'N','GX':'N','CX':'N','TX':'N'}
DP_FILTER=5
MQ_FILTER=20
GQ_FILTER=30


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

def read_locus_file_v2(file):              # locus interval is bed format, zero based, transform into 1-based
	locus=[]
	f=open(file,"r")
	for line in f:
		line = line.strip().split("\t")
		chr = line[0]
		pos1 = line[1]
		pos2 = line[2]
		locus.append([chr,int(pos1)+1,int(pos2)]) 
	return locus


def parse_vcf(vcf_file,locus,fasta_file):
	f = gzip.open(vcf_file,"r")
	i=0
	total_length=0
	het = 0
	string = []
	for j in range(locus[i][2]-locus[i][1]+1):
		string.append(0)
	for line in f:
		if line[0]=='#':
			continue
		line = line.strip().split("\t")
		info = line[7]
		chr = line[0]
		pos = int(line[1])
		m1 = re.search('DP=(\d+)',info)
		geno = line[9]
		MQ = int(line[5])
		genotype = geno.split(":")[0]
		GQ = int(geno.split(":")[2])
		if i>=len(locus):
			break
		if pos>=locus[i][1] and pos<=locus[i][2]:
			if string[pos-locus[i][1]]==1:
				string[pos-locus[i][1]]="N"
			else:
				if int(m1.group(1))<DP_FILTER or MQ<MQ_FILTER or GQ<GQ_FILTER:
					string[pos-locus[i][1]]="N"
				elif string[pos-locus[i][1]]!="N":
					if line[4]==".":					# reference allele
						string[pos-locus[i][1]]=line[3]
					else:	# non reference allele
						allele = line[4].split(",")
						if genotype=="0/1":
							string[pos-locus[i][1]]=annotation[allele[0]+line[3]]
							het +=1
						elif genotype=="1/1":
							string[pos-locus[i][1]]=allele[0]
						elif genotype=="1/2":
							string[pos-locus[i][1]]=annotation[allele[0]+allele[1]]
						else:
							string[pos-locus[i][1]]="N"
		elif pos>locus[i][2]:
			seq=""
			length = locus[i][2]-locus[i][1]+1
			for m in range(locus[i][1],locus[i][2]+1):
				if string[m-locus[i][1]]!=0 and string[m-locus[i][1]]!=1 and string[m-locus[i][1]]!=2:
					seq += string[m-locus[i][1]]
				else:
					seq += "N"
					length -=1
			print >>fasta_file[0],">%s:%i-%i" %(locus[i][0],locus[i][1]-1,locus[i][2])
			print >>fasta_file[0],seq
			total_length += length
			while True:
				i +=1
				if i>=len(locus):
					break
#				print locus[i],i,len(locus)
				string = []
				for j in range(locus[i][2]-locus[i][1]+1):
					string.append(0)
				if chr==locus[i][0] and pos>=locus[i][1] and pos<=locus[i][2]:
					if string[pos-locus[i][1]]==1:
						string[pos-locus[i][1]]="N"
					else:
						if int(m1.group(1))<DP_FILTER or MQ<MQ_FILTER or GQ<GQ_FILTER:
							string[pos-locus[i][1]]="N"
						elif string[pos-locus[i][1]]!="N":
							if line[4]==".":					# reference allele
								string[pos-locus[i][1]]=line[3]
							else:	# non reference allele
								allele = line[4].split(",")
								if genotype=="0/1":
									string[pos-locus[i][1]]=annotation[allele[0]+line[3]]
									het +=1
								elif genotype=="1/1":
									string[pos-locus[i][1]]=allele[0]
								elif genotype=="1/2":
									string[pos-locus[i][1]]=annotation[allele[0]+allele[1]]
								else:
									string[pos-locus[i][1]]="N"
					break
				else:
					seq=""
					length = locus[i][2]-locus[i][1]+1
					for m in range(locus[i][1],locus[i][2]+1):
						if string[m-locus[i][1]]!=0 and string[m-locus[i][1]]!=1:
							seq += string[m-locus[i][1]]
						else:
							seq += "N"
							length -=1
					print >>fasta_file[0],">%s:%i-%i" %(locus[i][0],locus[i][1]-1,locus[i][2])
					print >>fasta_file[0],seq
					total_length += length
	if i<len(locus):
		seq=""
		length = locus[i][2]-locus[i][1]+1
		for m in range(locus[i][1],locus[i][2]+1):
			if string[m-locus[i][1]]!=0 and string[m-locus[i][1]]!=1:
				seq += string[m-locus[i][1]]
			else:
				seq += "N"
				length -=1
		print >>fasta_file[0],">%s:%i-%i" %(locus[i][0],locus[i][1]-1,locus[i][2])
		print >>fasta_file[0],seq
		total_length += length

def parse_vcf_v2(vcf_file,locus,fasta_file):
	f = gzip.open(vcf_file,"r")
	i=0
	total_length=0
	het = 0
	prev_chr="None"
	for line in f:
		if line[0]=='#':
			continue
		line = line.strip().split("\t")
		info = line[7]
		chr = line[0]
		pos = int(line[1])
		m1 = re.search('DP=(\d+)',info)
		geno = line[9]
		ref = line[3]
		allele = line[4].split(",")
		MQ = line[5]
		if MQ==".":
			MQ=100
		else:
			MQ=float(MQ)
		genotype = geno.split(":")[0]
		GQ = int(geno.split(":")[2])
		if chr!=prev_chr:
			index=0
			if chr not in locus.keys():
				break
			loci=locus[chr][index]
			print_indicator = False
#			print loci
			string1 = []
			string2 = []
			for j in range(loci[2]-loci[1]+1):
				string1.append(0)
				string2.append(0)
			prev_chr=chr
		if pos>=loci[1] and pos<=loci[2]:
			if string1[pos-loci[1]]==1:
				string1[pos-loci[1]]="N"
			else:
				if int(m1.group(1))<DP_FILTER or GQ<GQ_FILTER or MQ<MQ_FILTER:
					string1[pos-loci[1]]="N"
				elif string1[pos-loci[1]]!="N":
					if ref=="N" or line[4]=="N":
						string1[pos-loci[1]]="N"
					elif line[4]==".":					# reference allele
						string1[pos-loci[1]]=line[3]
					else:	# non reference allele
						allele = line[4].split(",")
						if genotype=="0/1":
							string1[pos-loci[1]]=annotation[allele[0]+line[3]]
#							print pos,string[pos-loci[1]]
							het +=1
						elif genotype=="1/1":
							string1[pos-loci[1]]=allele[0]
#							print pos,string[pos-loci[1]]
						elif genotype=="1/2":
							string1[pos-loci[1]]=annotation[allele[0]+allele[1]]
						else:
							string1[pos-loci[1]]="N"
			if ref=="C" and allele[0]=="T" or ref=="G" and allele[0]=="A" or ref=="T" and allele[0]=="C" or ref=="A" and allele[0]=="G" :
				string2[pos-loci[1]]="N"
			else:
				string2[pos-loci[1]]=string1[pos-loci[1]]
		elif pos>loci[2]:
			if print_indicator is False:
				seq=""
				seq2=""
				length = loci[2]-loci[1]+1
				for m in range(loci[1],loci[2]+1):
					if string1[m-loci[1]]!=0:
						seq += string1[m-loci[1]]
						seq2 += string2[m-loci[1]]
					else:
						seq += "N"
						seq2 += "N"
						length -=1
		#			print seq
				print >>fasta_file[0],">%s:%i-%i" %(loci[0],loci[1]-1,loci[2])
				print >>fasta_file[0],seq
				print >>fasta_file[1],">%s:%i-%i" %(loci[0],loci[1]-1,loci[2])
				print >>fasta_file[1],seq2
				print_indicator = True
				total_length += length
			while True:
				index +=1
				if index>=len(locus[chr]):
					break
				loci=locus[chr][index]
				print_indicator = False
				string1 = []
				string2 = []
				for j in range(loci[2]-loci[1]+1):
					string1.append(0)
					string2.append(0)
				if pos>=loci[1] and pos<=loci[2]:
					if string1[pos-loci[1]]==1:
						string1[pos-loci[1]]="N"
					else:
						if int(m1.group(1))<DP_FILTER or GQ<GQ_FILTER or MQ<MQ_FILTER:
							string1[pos-loci[1]]="N"
						elif string1[pos-loci[1]]!="N":
							if ref=="N" or line[4]=="N":
								string1[pos-loci[1]]="N"
							elif line[4]==".":					# reference allele
								string1[pos-loci[1]]=line[3]
							else:	# non reference allele
								allele = line[4].split(",")
								if genotype=="0/1":
									string1[pos-loci[1]]=annotation[allele[0]+line[3]]
		#							print pos,string[pos-loci[1]]
									het +=1
								elif genotype=="1/1":
									string1[pos-loci[1]]=allele[0]
		#							print pos,string[pos-loci[1]]
								elif genotype=="1/2":
									string1[pos-loci[1]]=annotation[allele[0]+allele[1]]
								else:
									string1[pos-loci[1]]="N"
					if ref=="C" and allele[0]=="T" or ref=="G" and allele[0]=="A" or ref=="T" and allele[0]=="C" or ref=="A" and allele[0]=="G" :
						string2[pos-loci[1]]="N"
					else:
						string2[pos-loci[1]]=string1[pos-loci[1]]
					break
				elif pos<loci[1]:
					break
				else:
					if print_indicator is False:
						seq=""
						seq2=""
						length = loci[2]-loci[1]+1
						for m in range(loci[1],loci[2]+1):
							if string1[m-loci[1]]!=0:
								seq += string1[m-loci[1]]
								seq2 += string2[m-loci[1]]
							else:
								seq += "N"
								seq2 += "N"
								length -=1
		#					print seq
						print >>fasta_file[0],">%s:%i-%i" %(loci[0],loci[1]-1,loci[2])
						print >>fasta_file[0],seq
						print >>fasta_file[1],">%s:%i-%i" %(loci[0],loci[1]-1,loci[2])
						print >>fasta_file[1],seq2
						print_indicator = True
						total_length += length
	if print_indicator is False:
		seq=""
		seq2=""
		length = loci[2]-loci[1]+1
		for m in range(loci[1],loci[2]+1):
			if string1[m-loci[1]]!=0:
				seq += string1[m-loci[1]]
				seq2 += string2[m-loci[1]]
			else:
				seq += "N"
				seq2 += "N"
				length -=1
	#		print seq
		print >>fasta_file[0],">%s:%i-%i" %(loci[0],loci[1]-1,loci[2])
		print >>fasta_file[0],seq
		print >>fasta_file[1],">%s:%i-%i" %(loci[0],loci[1]-1,loci[2])
		print >>fasta_file[1],seq2
		total_length += length
	print het,total_length,float(het)/total_length
	return het,total_length

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
	parser.add_argument("--vcf_dir",default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/GPhoCS/HC',dest='vcf_dir',help="vcf dir")
	args = parser.parse_args()		
	locus_file = "/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/neutralLoci-geneFlank10k-sep30k.bed"
	locus=read_locus_file(locus_file)
	string="MQ%s_DP%s_GQ%s" %(MQ_FILTER,DP_FILTER,GQ_FILTER)
	fasta_file = [open("%s_%s.fa" %(args.sample,string),"w"),open("%s_%s_rm_damage.fa" %(args.sample,string),"w")]
	vcf_file = "%s/%s_HC_GPhoCS.vcf.gz" %(args.vcf_dir,args.sample)
	het,total_length=parse_vcf_v2(vcf_file,locus,fasta_file)
	outfile=open("heterozygosity_canfam3_%s.txt" %(string),"a")
	print >>outfile,args.sample,het,total_length,float(het)/total_length

			
		
	

 
