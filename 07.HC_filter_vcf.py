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
import pysam

def parse_gvcf(vcf_file,ref_seq):
	f = gzip.open(vcf_file,"r")
	i=0
	total_length=0
	het = 0
	prev_chr="None"
	for line in f:
		if line[0]=='#':
			print line.rstrip()
			continue
		line = line.rstrip().split('\t')
		geno_info=line[9].split(':')
		format=line[8].split(":")
		format={m:j for j,m in enumerate(format)}
		if 'GQ' in format.keys():
			GQ=int(geno_info[format['GQ']])
		else:
			GQ=100
		if 'DP' in format.keys():
			DP=int(geno_info[format['DP']])
		else:
			DP = 0
		if 'PL' in format.keys():
			PL=geno_info[format['PL']]
		else:
			PL = 0
		info=line[7]
		m1 = re.search('MQ=(\d+)',info)
		if m1!=None:
			MQ=int(m1.group(1))
		else:
			MQ=0

		chr = line[0]
		pos = int(line[1])
		genotype = geno_info[0]
		ref = line[3]
		altAllele=line[4]
		if ref=='N' or DP==0:
			continue
		if altAllele=="<NON_REF>":	 # reference block
			end_pos_match = re.search("END=(\d+)",info)
			if end_pos_match:
				end_pos = int(end_pos_match.group(1))
			else:
#				print 'End not found:',line
				end_pos = pos
			for each_pos in range(pos,end_pos+1):
				ref_allele=ref_seq.fetch(chr,int(each_pos-1),int(each_pos))
				ref_allele=ref_allele.upper()
				line[1]=each_pos
				line[3]=ref_allele
				line[4]='.'
				line[7]="DP=%s" %(DP)
				line[8]="GT:PL:GQ"
				line[9]="%s:%s:%s" %(genotype,PL,GQ)
				print "\t".join(map(str,line))
		elif re.match("^[ACTGactg]$", ref) and re.match("^[ACTGactg],<NON_REF>$", altAllele):		# snps, exclude indels
			dp_match = re.search("DP=(\d+)", info)
			if not dp_match:
				continue
			DP = int(dp_match.group(1))
			line[4]=line[4].replace(',<NON_REF>','')
			line[7]="DP=%s" %(DP)
			line[8]="GT:PL:GQ"
			line[9]="%s:%s:%s" %(genotype,PL,GQ)
			print "\t".join(line)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='vcf to fa')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--vcf_dir",default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/GPhoCS/BSNP',dest='vcf_dir',help="vcf dir")
	args = parser.parse_args()		
	
	ref_file='/home/jmkidd/kidd-lab/genomes/canFam3/canFam3-cat/canFam3.fa'
	ref_seq=pysam.FastaFile(ref_file)
	parse_gvcf("%s_GPhoCS_gvcf.gz" %(args.sample),ref_seq)

			
		
	

 
