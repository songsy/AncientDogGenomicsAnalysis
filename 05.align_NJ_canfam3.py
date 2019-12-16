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

annotation={'AG':'R','GA':'R','TC':'Y','CT':'Y','GT':'K','TG':'K','AC':'M','CA':'M','GC':'S','CG':'S','AT':'W','TA':'W','XA':'N','XG':'N','XC':'N','XT':'N','AX':'N','GX':'N','CX':'N','TX':'N'}

def merge_seq(seq1,seq2):
	seq=""
	for i in range(len(seq1)):
		if seq1[i]!="N" and seq2[i]!="N":
			if seq1[i]==seq2[i]:
				seq+=seq1[i]
			else:
				seq+=annotation[seq1[i]+seq2[i]]
		else:
			seq+="N"

	return seq

def delta_function(a,b):
	if a==b:
		return 1
	else:
		return 0

def read_GhoCSseq_GQ_20(file,locus,SEQ):
	f=open(file,'r')
	for index,line in enumerate(f):
		line=line.rstrip().split('\t')
		info=line[0]
		seq=line[1]
		tag=info.split(":")[2]

		if index%2==0:
			seq1=line[1]
		else:
			seq2=line[1]
			SEQ[tag]=merge_seq(seq1,seq2)
	return SEQ

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

def read_fasta(file,locus,SEQ,ref):
	annotation={'R':'GA','Y':'CT','K':'TG','M':'CA','S':'CG','W':'TA','A':'AA','T':'TT','C':'CC','G':'GG','N':'NN'}
	f=open(file,"r")
	mut_type={}
	het = 0
	Length = 0
	for line in f:
		line=line.strip()
		if line[0]==">":
			tag="_".join(line[1:].split(":"))
			tag=tag.replace("-","_")
			chrom,start,end=tag.split("_")
			ref_seq=ref.fetch(chrom,int(start),int(end))
			ref_seq=ref_seq.upper()
		else:
			SEQ[tag]=line
			seq1=[annotation[line[m]][0] for m in range(len(line))]
			seq2=[annotation[line[m]][1] for m in range(len(line))]
			het,Length,mut_type=merge_seq_v2(seq1,seq2,ref_seq,het,Length,mut_type)
	return SEQ,het,Length,mut_type

def merge_seq_v2(seq1,seq2,ref_seq,het,Length,mut_type):
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

def remove_CG(record):
	annotation={'R':'GA','Y':'CT','K':'TG','M':'CA','S':'CG','W':'TA','A':'AA','T':'TT','C':'CC','G':'GG','N':'NN'}
	pos = []
	for i in record:
		for j in range(len(record[i])-1):
			if "C" in annotation[record[i][j]] and "G" in annotation[record[i][j+1]]:
				pos.append(j)
				pos.append(j+1)
	for i in record:
		seq = ""
		for j in range(len(record[i])):
			if j in pos:
				seq+="N"
			else:
				seq +=record[i][j]
		record[i]=seq
	return record

def build_record(SEQ,loci_name,length):
	record={}
	for sample in SEQ.keys():
		try:
			record[sample]=SEQ[sample][loci_name]
		except KeyError:
			record[sample]="N"*length
	return record

def comparison(seq1,seq2):  # E8.1
	annotation={'R':'GA','Y':'CT','K':'TG','M':'CA','S':'CG','W':'TA','A':'AA','T':'TT','C':'CC','G':'GG'}
	assert len(seq1)==len(seq2)
	dist = 0
	length = 0
	for i in range(len(seq1)):
		if seq1[i]!="N" and seq2[i]!="N":
			a=annotation[seq1[i]]
			b=annotation[seq2[i]]
			dist+= 1- float(0.5)*max(delta_function(a[0],b[0])+delta_function(a[1],b[1]),delta_function(a[0],b[1])+delta_function(a[1],b[0]))
			length+=1
	return dist,length

def comparison_v2(seq1,seq2):  #E8.2
	annotation={'R':'GA','Y':'CT','K':'TG','M':'CA','S':'CG','W':'TA','A':'AA','T':'TT','C':'CC','G':'GG'}
	assert len(seq1)==len(seq2)
	dist = 0
	length = 0
	for i in range(len(seq1)):
		if seq1[i]!="N" and seq2[i]!="N":
			a=annotation[seq1[i]]
			b=annotation[seq2[i]]
			dist+= 1- float(0.25)*(delta_function(a[0],b[0])+delta_function(a[1],b[1])+delta_function(a[0],b[1])+delta_function(a[1],b[0]))
			length+=1
	return dist,length

def aln_to_fasta(SEQ,locus):
	sample_name=sorted(SEQ.keys())
	matrix=np.zeros((len(sample_name),len(sample_name)))
	size=np.zeros((len(sample_name),len(sample_name)))
	f=open('Gphocs_input_HC_canfam3_MQ20_DP5_GQ30.txt',"w")
#	f=open('Gphocs_input_canfam3_GQ30_BQ5_MQ15_DP5_unrecal_aDNA_varyGQ.txt',"w")
	remove_locus=[]
	print >>f,"26250"
	f.write("\n")
	f_matrix=open("dist_matrix_HC_canfam3_MQ20_DP5_GQ30.txt","w")
#	f_matrix=open("dist_matrix_canfam3_GQ30_BQ5_MQ15_DP5_unrecal_aDNA_varyGQ.txt","w")
	for chr in sorted(locus.keys()):
		print chr
		for loci in locus[chr]:
			zero=False
			pos1=loci[1]
			pos2=loci[2]
			loci_name="%s_%s_%s" %(chr,pos1-1,pos2)
			record=build_record(SEQ,loci_name,pos2-pos1+1)
			record=remove_CG(record)

			for j in range(len(sample_name)):
				for k in range(j+1,len(sample_name)):
					dist,length=comparison(record[sample_name[j]],record[sample_name[k]])
					if length==0:
						zero = True
					matrix[j][k]+=dist
					size[j][k]+=length
					matrix[k][j]+=dist
					size[k][j]+=length
			if zero is False:
				print >>f,"%s.%s-%s\t%i\t%s" %(chr,pos1-1,pos2,len(sample_name),pos2-pos1+1)
				for j in sample_name:
					print >>f,"%s\t%s" %(j,record[j])
				f.write("\n")
			else:
				remove_locus.append(loci)
				print "remove locus",loci
	print >>f_matrix,"dist_matrix:"
	print >>f_matrix,"\t".join(sample_name)
	for i in range(len(sample_name)):
		print >>f_matrix,"\t".join(map(str,matrix[i]))
	print >>f_matrix,"length_matrix:"
	print >>f_matrix,"\t".join(sample_name)
	for i in range(len(sample_name)):
		print >>f_matrix,"\t".join(map(str,size[i]))
	print >>f_matrix,"distance_matrix"
	print >>f_matrix,"\t".join(sample_name)
	for i in range(len(sample_name)):
		print >>f_matrix,"\t".join(map(str,[x/y for x,y in zip(matrix[i],size[i])]))
	return remove_locus

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='vcf to fa')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--vcf_dir",default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/GPhoCS/HC',dest='vcf_dir',help="vcf dir")
	args = parser.parse_args()		
	locus_file = "/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/neutralLoci-geneFlank10k-sep30k.bed"
	locus=read_locus_file(locus_file)
	SEQ={}
# READ SEQ DATA:
	ref_file='/home/jmkidd/kidd-lab/genomes/canFam3/canFam3-cat/canFam3.fa'
	ref=pysam.FastaFile(ref_file)

	'''
	for name in ['Basenji','Dingo','IsraeliWolf','CroatianWolf','ChineseWolf','GoldenJackal']:
		SEQ[name]={}
		file='%s/%s.GhoCSseq_GQ30_BQ5_MQ15_DP5' %(args.vcf_dir,name)
		SEQ[name]=read_GhoCSseq_GQ_20(file,locus,SEQ[name])
		print name,len(SEQ[name]),sorted(SEQ[name].keys())[0]
	
	file='locus-Boxer.fa'
	SEQ['Boxer']={}
	SEQ['Boxer'],het,Length,mut_type=read_fasta(file,locus,SEQ['Boxer'],ref)
	
	for value in [10,20,30,40]:
		for name in ['Kirshbaum','Herxheim']:
			SEQ[name+"_%s" %(value)]={}
			file='%s/%s_unrecal.GhoCSseq_aDNA_GQ%s_BQ5_MQ15_DP5' %(args.vcf_dir,name,value)
			SEQ[name+"_%s" %(value)]=read_GhoCSseq_GQ_20(file,locus,SEQ[name+"_%s" %(value)])
			print name,len(SEQ[name]),sorted(SEQ[name].keys())[0]

	remove_locus=aln_to_fasta(SEQ,locus)
	print 'removed_locus',len(remove_locus)
	'''
	fout=open("heterozygosity_canfam3_MQ20_DP5_GQ30_detail.txt","a")
	for name in ['Basenji','Dingo','IsraeliWolf','CroatianWolf','ChineseWolf','GoldenJackal','CTC','HxH']:
		SEQ[name]={}
		file='%s/%s_MQ20_DP5_GQ30.fa' %(args.vcf_dir,name)
		SEQ[name],het,Length,mut_type=read_fasta(file,locus,SEQ[name],ref)
		print name,len(SEQ[name]),sorted(SEQ[name].keys())[0]
		a=[mut_type[i] for i in sorted(mut_type.keys())]
		print >>fout,name,"\t",het,"\t",Length,"\t",float(het)/Length,"\t","\t".join(map(str,a))
	file='locus-Boxer.fa'
	SEQ['Boxer']={}
	SEQ['Boxer'],het,Length,mut_type=read_fasta(file,locus,SEQ['Boxer'],ref)
	remove_locus=aln_to_fasta(SEQ,locus)
	print 'removed_locus',len(remove_locus)

