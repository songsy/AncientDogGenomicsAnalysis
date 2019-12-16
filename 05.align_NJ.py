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

def read_GhoCSseq_GQ_20(file,locus,trf,rmsk,SEQ):
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

def read_fasta(file,locus,trf,rmsk,SEQ):
	f=open(file,"r")
	for line in f:
		line=line.strip()
		if line[0]==">":
			tag="_".join(line[1:].split(":"))
			tag=tag.replace("-","_")
		else:
			SEQ[tag]=line
	return SEQ

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

def comparison(seq1,seq2):
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

def remove_repeat(record,trf_mask,rmsk_mask,chr,pos1,pos2):
	remove_list=[]
	for pos in range(pos1,pos2+1):
		if pos in trf_mask or pos in rmsk_mask:
			remove_list.append(pos-pos1)
	if len(remove_list)==0:
		return record
	else:
		print chr,pos1,pos2,len(remove_list)
	for i in record:
		seq = ""
		for j in range(len(record[i])):
			if j in remove_list:
				seq+="N"
			else:
				seq +=record[i][j]
		record[i]=seq
	return record

def aln_to_fasta(SEQ,locus,trf,rmsk):
	sample_name=sorted(SEQ.keys())
	matrix=np.zeros((len(sample_name),len(sample_name)))
	size=np.zeros((len(sample_name),len(sample_name)))
	f=open('Gphocs_input_v3.txt',"w")
	remove_locus=[]
	print >>f,"26250"
	f.write("\n")
	f_matrix=open("dist_matrix_v3.txt","w")
	for chr in sorted(locus.keys()):
		chr_id='chr'+str(int(chr[3:]))
		if chr_id in trf.keys():
			trf_mask=set(chain(*(xrange(start, end+1) for chrom,start,end in trf[chr_id])))
		else:
			trf_mask=[]
		if chr_id in rmsk.keys():
			rmsk_mask=set(chain(*(xrange(start, end+1) for chrom,start,end in rmsk[chr_id])))
		else:
			rmsk_mask=[]
		print chr,len(trf_mask),len(rmsk_mask)
		for loci in locus[chr]:
			zero=False
			pos1=loci[1]
			pos2=loci[2]
			loci_name="%s_%s_%s" %(chr,pos1-1,pos2)
			record=build_record(SEQ,loci_name,pos2-pos1+1)
#			record=remove_repeat(record,trf_mask,rmsk_mask,chr,pos1,pos2)
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


if __name__=="__main__":
	parser = argparse.ArgumentParser(description='vcf to fa')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--vcf_dir",default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/GPhoCS/BSNP/',dest='vcf_dir',help="vcf dir")
	args = parser.parse_args()		
	trf_file = "/home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1_trf.bed"
	rmsk_file = "/home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1_rmsk.bed"
	locus_file = "/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/neutralLoci-geneFlank10k-sep30k-3.1.bed"
	locus=read_locus_file(locus_file)
	trf = read_trf_file(trf_file)
	rmsk = read_trf_file(rmsk_file)
#	fasta_file = [open("locus-%s.fa" %(args.sample),"w"),open("locus-%s-rm-damage.fa" %(args.sample),"w")]

	SEQ={}
# READ SEQ DATA:
	'''
	for name in ['Basenji','Dingo','IsraeliWolf','CroatianWolf','ChineseWolf','GoldenJackal','Kirshbaum']:
		SEQ[name]={}
		file='%s.GhoCSseq_GQ_20' %(name)
		SEQ[name]=read_GhoCSseq_GQ_20(file,locus,trf,rmsk,SEQ[name])
		print name,len(SEQ[name]),sorted(SEQ[name].keys())[0]
	file='Boxer_ref.fa'
	SEQ['Boxer']={}
	read_fasta(file,locus,trf,rmsk,SEQ['Boxer'])

	for file in glob.glob('locus-*-rm-damage.fa'):
		name=file.split('-')[1]
		SEQ[name]={}
		file='locus-%s.fa' %(name)
		SEQ[name]=read_fasta(file,locus,trf,rmsk,SEQ[name])
		print name,len(SEQ[name]),sorted(SEQ[name].keys())[0]
	'''
	for name in ['Basenji','Dingo','IsraeliWolf','CroatianWolf','ChineseWolf','GoldenJackal','CTC','HxH']:
		SEQ[name]={}
		file='%s_MQ0_DP0_GQ20.fa' %(name)
		SEQ[name]=read_fasta(file,locus,trf,rmsk,SEQ[name])
		print name,len(SEQ[name]),sorted(SEQ[name].keys())[0]

	file='locus-Boxer.fa'
	SEQ['Boxer']={}
	read_fasta(file,locus,SEQ['Boxer'])
	remove_locus=aln_to_fasta(SEQ,locus,trf,rmsk)
	print 'removed_locus',len(remove_locus)
