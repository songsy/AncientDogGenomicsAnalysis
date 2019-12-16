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

def merge_seq(seq1,seq2,het,Length):
	seq=""
	for i in range(len(seq1)):
		if seq1[i]!="N" and seq2[i]!="N":
			if seq1[i]==seq2[i]:
				seq+=seq1[i]
			else:
				seq+=annotation[seq1[i]+seq2[i]]
				het+=1
			Length+=1
		else:
			seq+="N"

	return seq,het,Length

def delta_function(a,b):
	if a==b:
		return 1
	else:
		return 0

def read_GhoCSseq_GQ_20(file,locus,SEQ):
	f=open(file,'r')
	het=0
	Length=0
	for index,line in enumerate(f):
		line=line.rstrip().split('\t')
		info=line[0]
		seq=line[1]
		tag=info.split(":")[2]

		if index%2==0:
			seq1=line[1]
		else:
			seq2=line[1]
			SEQ[tag],het,Length=merge_seq(seq1,seq2,het,Length)
	return SEQ,het,Length

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

def calc_singleton(record):
	singleton=np.array([0 for i in sorted(record.keys())])
	dict={'A':0,'C':1,'G':2,'T':3,'N':4}
	annotation={'R':'GA','Y':'CT','K':'TG','M':'CA','S':'CG','W':'TA','A':'AA','T':'TT','C':'CC','G':'GG','N':'NN'}
	for j in range(len(record[record.keys()[0]])):
		allele=""
		for i in sorted(record.keys()):
			allele+=annotation[record[i][j]]
		count=np.array([allele.count("A"),allele.count("C"),allele.count("G"),allele.count("T"),allele.count("N")])
		if count[4]>0:
			continue
		for index,i in enumerate(sorted(record.keys())):
			single= False
			a=annotation[record[i][j]]
			count1=np.array([a.count("A"),a.count("C"),a.count("G"),a.count("T"),a.count("N")])
			count2=count-count1
			for k in a:
				assert count2[dict[k]]>=0
				if count2[dict[k]]==0:
					single = True
			if single:
				singleton[index]+=1
	return singleton
			
def remove_damage(record,ref_seq):
	annotation={'R':'GA','Y':'CT','K':'TG','M':'CA','S':'CG','W':'TA','A':'AA','T':'TT','C':'CC','G':'GG','N':'NN'}
	pos = []
	for i in record:
		assert len(record[i])==ref_seq
		for j in range(len(record[i])):
			if (ref_seq=="C" and "T" in annotation[record[i][j]]) or (ref_seq=="T" and "C" in annotation[record[i][j]]):
				pos.append(j)
			elif (ref_seq=="G" and "A" in annotation[record[i][j]]) or (ref_seq=="A" and "G" in annotation[record[i][j]]):
				pos.append(j)
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
#	f=open('Gphocs_input_canfam3_MQ0_DP0_GQ20.txt',"w")
	f=open('Gphocs_input_canfam3.1_subset2_GQ30_BQ15_MQ15_DP7_unrecal_aDNA.txt',"w")
	remove_locus=[]
	print >>f,"26250"
	f.write("\n")
#	f_matrix=open("dist_matrix_canfam3_MQ0_DP0_GQ20.txt","w")
	f_matrix=open("dist_matrix_canfam3.1_subset2_GQ30_BQ15_MQ15_DP7_unrecal_aDNA.txt","w")
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

def aln_to_fasta_v2(SEQ,locus):  # contain the result for removing potential damaging bases
	sample_name=sorted(SEQ.keys())
	matrix=np.zeros((len(sample_name),len(sample_name)))
	size=np.zeros((len(sample_name),len(sample_name)))
#	f=open('Gphocs_input_canfam3_MQ0_DP0_GQ20.txt',"w")
#	f=open('Gphocs_input_canfam3.1_all_sample_GQ30_BQ5_MQ15_DP5_unrecal_aDNA.txt',"w")
	f=open('Gphocs_input_canfam3.1_all_sample_GQ30_BQ5_MQ15_DP5_unrecal_aDNA_rm_damage.txt',"w")
	remove_locus=[]
	print >>f,"26250"
	f.write("\n")
#	f_matrix=open("dist_matrix_canfam3_MQ0_DP0_GQ20.txt","w")
	f_matrix=open("dist_matrix_canfam3.1_all_sample_GQ30_BQ5_MQ15_DP5_unrecal_aDNA_rm_damage.txt","w")
	for chr in sorted(locus.keys()):
		print chr
		for loci in locus[chr]:
			zero=False
			pos1=loci[1]
			pos2=loci[2]
			loci_name="%s_%s_%s" %(chr,pos1-1,pos2)
			ref_seq=ref.fetch(chr,pos1-1,pos2)
			ref_seq=ref_seq.upper()
			record=build_record(SEQ,loci_name,pos2-pos1+1)
			record=remove_CG(record)
			record=remove_damage(record,ref_seq)
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

def calc_ratio(record):
	singleton=np.array([0 for i in sorted(record.keys())])
	dict={'A':0,'C':1,'G':2,'T':3,'N':4}
	annotation={'R':'GA','Y':'CT','K':'TG','M':'CA','S':'CG','W':'TA','A':'AA','T':'TT','C':'CC','G':'GG','N':'NN'}
	AB=0
	AC=0
	for j in range(len(record[record.keys()[0]])):
		if record["GoldenJackal"][j]=='N' or record["VillageIndian2"][j]=='N' or record["VillageEuropePT"][j]=='N' or record["NewGrange"][j]=='N':
			continue
		jackal=annotation[record["GoldenJackal"][j]]
		Indian=annotation[record["VillageIndian2"][j]]
		Europe=annotation[record["VillageEuropePT"][j]]
		HXH=annotation[record["NewGrange"][j]]
		index=[random.randint(0,1) for k in range(4)]
		if jackal[index[0]]==Indian[index[1]] and HXH[index[3]]==Europe[index[2]] and HXH[index[3]]!=Indian[index[1]]:
			AB+=1
			print 'AB',jackal[index[0]],Indian[index[1]],Europe[index[2]],HXH[index[3]]
		elif jackal[index[0]]==HXH[index[3]] and Indian[index[1]]==Europe[index[2]] and HXH[index[3]]!=Indian[index[1]]:
			AC+=1
			print 'AC',jackal[index[0]],Indian[index[1]],Europe[index[2]],HXH[index[3]]
	return AB,AC

def aln_to_fasta_v3(SEQ,locus):  # calc singletons
	sample_name=sorted(SEQ.keys())
	print sample_name
	matrix=np.zeros((len(sample_name),len(sample_name)))
	size=np.zeros((len(sample_name),len(sample_name)))
#	f=open('Gphocs_input_canfam3_MQ0_DP0_GQ20.txt',"w")
#	f=open('Gphocs_input_canfam3.1_subset2_GQ30_BQ15_MQ15_DP7_unrecal_aDNA.txt',"w")
	remove_locus=[]
#	print >>f,"26250"
#	f.write("\n")
#	f_matrix=open("dist_matrix_canfam3_MQ0_DP0_GQ20.txt","w")
#	f_matrix=open("dist_matrix_canfam3.1_subset2_GQ30_BQ15_MQ15_DP7_unrecal_aDNA.txt","w")
	Singleton=np.array([0 for i in sample_name])
	AB=0
	AC=0
	for chr in sorted(locus.keys()):
		print chr
		for loci in locus[chr]:
			zero=False
			pos1=loci[1]
			pos2=loci[2]
			loci_name="%s_%s_%s" %(chr,pos1-1,pos2)
			record=build_record(SEQ,loci_name,pos2-pos1+1)
			record=remove_CG(record)
			singleton=calc_singleton(record)
#			ab,ac=calc_ratio(record)
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
#				print >>f,"%s.%s-%s\t%i\t%s" %(chr,pos1-1,pos2,len(sample_name),pos2-pos1+1)
#				for j in sample_name:
#					print >>f,"%s\t%s" %(j,record[j])
#				f.write("\n")
				Singleton=Singleton+singleton
#				print Singleton
#				AB+=ab
#				AC+=ac
			else:
				remove_locus.append(loci)
				print "remove locus",loci
	'''
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
	'''
	print sample_name
	print Singleton
#	print AB,AC,float(AB)/AC
	return remove_locus
	
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='vcf to fa')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--vcf_dir",default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/GPhoCS/input_3.1_v2',dest='vcf_dir',help="vcf dir")
	args = parser.parse_args()		
	locus_file = "/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/neutralLoci-geneFlank10k-sep30k-3.1.bed"
	locus=read_locus_file(locus_file)
	SEQ={}

#	print 12
	ref_file='/home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1-withUn/canFam3.1.withUn.fa'
	ref=pysam.FastaFile(ref_file)

	f=open('sampleList_subset2.txt','r')
	for line in f:
		a=line.rstrip().split('\t')
		name=a[1]
#	for name in ['zoey','Basenji','Dingo','IsraeliWolf','CroatianWolf','ChineseWolf','GoldenJackal']:
		SEQ[name]={}
		file='%s/%s.GhoCSseq_GQ30_BQ15_MQ15_DP7' %(args.vcf_dir,a[0])
		SEQ[name],het,Length=read_GhoCSseq_GQ_20(file,locus,SEQ[name])
#		print name,len(SEQ[name]),sorted(SEQ[name].keys())[0]
		print name,het,Length,float(het)/Length
	

	for value in [30]:
#		for name in ['Kirshbaum','Herxheim']:
#		for name in ['Herxheim']:
		for name in ['NewGrange']:
			SEQ[name]={}
			file='%s/%s_unrecal.GhoCSseq_aDNA_GQ%s_BQ15_MQ15_DP7' %(args.vcf_dir,name,value)
			SEQ[name],het,Length=read_GhoCSseq_GQ_20(file,locus,SEQ[name])
#			print name,len(SEQ[name+"_%s" %(value)]),sorted(SEQ[name+"_%s" %(value)].keys())[0]
			print name,het,Length,float(het)/Length

	remove_locus=aln_to_fasta_v3(SEQ,locus) 
	print 'removed_locus',len(remove_locus)
