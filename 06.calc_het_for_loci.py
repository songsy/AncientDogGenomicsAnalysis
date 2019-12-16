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
import gzip
from numpy.random import binomial
import argparse
from itertools import chain
import glob
import pysam

def merge_seq(seq1,seq2,ref_seq,het,Length,mut_type):
	assert len(ref_seq)==len(seq1)
	for i in range(len(seq1)):
		if seq1[i]!="N" and seq2[i]!="N":
			if seq1[i]==seq2[i]:
				Length+=1
			else:
				het+=1
			if seq1[i]!=ref_seq[i]:
				if ref_seq[i]+seq1[i] not in mut_type.keys():
					mut_type[ref_seq[i]+seq1[i]]=0
				mut_type[ref_seq[i]+seq1[i]]+=1
			if seq2[i]!=ref_seq[i]:
				if ref_seq[i]+seq2[i] not in mut_type.keys():
					mut_type[ref_seq[i]+seq2[i]]=0
				mut_type[ref_seq[i]+seq2[i]]+=1
	return het,Length,mut_type

def calc_het(file,ref):
	mut_type={}
	f=open(file,'r')
	het = 0
	Length = 0
	for index,line in enumerate(f):
		line=line.rstrip().split('\t')
		info=line[0]
		seq=line[1]
		tag=info.split(:)[2]
		chrom,start,end=tag.split(_)
		ref_seq=ref.fetch(chrom,int(start),int(end))
		ref_seq=ref_seq.upper()
		if index%2==0:
			seq1=line[1]
		else:
			seq2=line[1]
			het,Length,mut_type=merge_seq(seq1,seq2,ref_seq,het,Length,mut_type)
	return het,Length,mut_type

def calc_het_from_vcf():
	transition={}
	het ={}
	f=gzip.open('CTC_BSNP.vcf.gz','r')
	for line in f:
		if line[0]=='#':
			continue
		line = line.rstrip().split('\t')
		geno=line[9].split(':')[0]
		if geno=="1/1" or geno=="0/1":
			if len(line[4])>1:
				 continue
			ref_alt=line[3]+line[4]
			DP=int(line[7].split('=')[1])
			GQ=int(line[9].split(':')[2])
			if ref_alt not in transition.keys():
				transition[ref_alt]=[]
			transition[ref_alt].append([DP,GQ])
			if geno=="0/1":
				if ref_alt not in het.keys():
					het[ref_alt]=[]
				het[ref_alt].append([DP,GQ])
	for i in transition.keys():
		print i,len(transition[i])
	for i in het.keys():
		print i,len(het[i])
	return het,transition

if __name__==__main__:
	parser = argparse.ArgumentParser(description='vcf to fa')
	parser.add_argument(--sample, dest='sample',help=sample name)
	parser.add_argument(--vcf_dir,default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/GPhoCS/summary_stats',dest='vcf_dir',help=vcf dir)
	args = parser.parse_args()

	ref_file='/home/jmkidd/kidd-lab/genomes/canFam3/canFam3-cat/canFam3.fa'
	ref=pysam.FastaFile(ref_file)
	for file in glob.glob('*GhoCSseq_aDNArefbias*'):	
		if file!='bam_list_aDNA.GhoCSseq_GQ_20' and file!='CanMap3.GhoCSseq':
			het,Length,mut_type=calc_het(file,ref)
			a=[mut_type[i] for i in sorted(mut_type.keys())]
			if Length==0:
				continue
			print mut_type
			print file,"\t",het,"\t",Length,"\t",float(het)/Length,"\t","\t".join(map(str,a))
