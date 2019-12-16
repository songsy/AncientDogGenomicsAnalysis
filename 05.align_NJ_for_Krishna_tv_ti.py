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
from itertools import chain
import pysam
import random

annotation={'AG':'R','GA':'R','TC':'Y','CT':'Y','GT':'K','TG':'K','AC':'M','CA':'M','GC':'S','CG':'S','AT':'W','TA':'W','XA':'N','XG':'N','XC':'N','XT':'N','AX':'N','GX':'N','CX':'N','TX':'N'}

def merge_seq(seq1,seq2,ti,tv,Length,ref_seq):
	seq=""
	assert len(seq1)==len(ref_seq)
	for i in range(len(seq1)):
		if seq1[i]!="N" and seq2[i]!="N":
			if seq1[i]==seq2[i]:
				seq+=seq1[i]
				if ref_seq[i]!=seq1[i]:
					if annotation[seq1[i]+ref_seq[i]] in ['M','S','W','K']:
						tv[0]+=1
					elif annotation[seq1[i]+ref_seq[i]] in ['R','Y']:
						ti[0]+=1
			else:
				seq+=annotation[seq1[i]+seq2[i]]
				if annotation[seq1[i]+seq2[i]] in ['M','S','W','K']:
					tv[1]+=1
				elif annotation[seq1[i]+seq2[i]] in ['R','Y']:
					ti[1]+=1
			Length+=1
		else:
			seq+="N"
	return seq,ti,tv,Length

def delta_function(a,b):
	if a==b:
		return 1
	else:
		return 0

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
	
def read_GhoCSseq_GQ_20(file,locus,ref,SEQ):
	f=open(file,'r')
	ti=[0,0] # homo,het
	tv=[0,0] # homo,het
	Length=0
	for index,line in enumerate(f):
		line=line.rstrip().split('\t')
		info=line[0]
		seq=line[1]
		tag=info.split(":")[2]
		chrom,start,end=tag.split("_")
		ref_seq=ref.fetch(chrom,int(start),int(end))
		ref_seq=ref_seq.upper()
		if index%2==0:
			seq1=line[1]
		else:
			seq2=line[1]
			SEQ[tag],ti,tv,Length=merge_seq(seq1,seq2,ti,tv,Length,ref_seq)
	return SEQ,ti,tv,Length

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='vcf to fa')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--vcf_dir",default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/GPhoCS/input_3.1_v2',dest='vcf_dir',help="vcf dir")
	args = parser.parse_args()		
	locus_file = "/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/neutralLoci-geneFlank10k-sep30k-3.1.bed"
	locus=read_locus_file(locus_file)
	SEQ={}

	ref_file='/home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1-withUn/canFam3.1.withUn.fa'
	ref=pysam.FastaFile(ref_file)
	
	f=open('sampleList.txt','r')
	for line in f:
		a=line.rstrip().split('\t')
		name=a[1]
		file='%s/%s.GhoCSseq_GQ30_BQ15_MQ15_DP7' %(args.vcf_dir,name)
		SEQ={}
		SEQ,ti,tv,Length=read_GhoCSseq_GQ_20(file,locus,ref,SEQ)
		print name,",".join(map(str,ti)),",".join(map(str,tv)),float(tv[0])/ti[0],float(tv[1])/ti[1]
		
	for value in [30]:
#		for name in ['Kirshbaum','Herxheim']:
#		for name in ['Herxheim']:
		for name in ['Kirshbaum','Herxheim','NewGrange']:
			SEQ[name]={}
			file='%s/%s_unrecal.GhoCSseq_aDNA_GQ%s_BQ15_MQ15_DP7' %(args.vcf_dir,name,value)
			SEQ={}
			SEQ,ti,tv,Length=read_GhoCSseq_GQ_20(file,locus,ref,SEQ)
			print name,",".join(map(str,ti)),",".join(map(str,tv)),float(tv[0])/ti[0],float(tv[1])/ti[1]