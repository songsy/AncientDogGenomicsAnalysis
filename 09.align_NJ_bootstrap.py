#!/usr/bin/env python
# python 03.parse_interval_vcf.py
# Shiya Song
# 3rd June 2013
# Parse vcf file, mark indel(around 3bp) as N, mark those with DP<5 as N, mark simple repeats as N, record clustered SNP

import re
import numpy as np
import gzip
from numpy.random import binomial
import argparse,time,itertools,multiprocessing,math,glob,os,subprocess,pickle,sys
from itertools import chain
import pysam

annotation={'AG':'R','GA':'R','TC':'Y','CT':'Y','GT':'K','TG':'K','AC':'M','CA':'M','GC':'S','CG':'S','AT':'W','TA':'W','XA':'N','XG':'N','XC':'N','XT':'N','AX':'N','GX':'N','CX':'N','TX':'N'}

def get_sample_name(file):
	sample_list=[]
	f=open(file,'r')
	for line in f:
		name=line.rstrip()
		sample_list.append(name)
	return sample_list

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

def read_GhoCSseq_GQ_20(file,locus,SEQ,i):   
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
			merged_seq,het,Length=merge_seq(seq1,seq2,het,Length)
			assert len(merged_seq)==len(SEQ[tag][i])
			SEQ[tag][i]=merged_seq
	return SEQ,het,Length

def read_locus_file(file,SEQ,sample_list):              # locus interval is bed format, zero based, transform into 1-based
	locus=[]
	f=open(file,"r")
	for line in f:
		line = line.strip().split("\t")
		chr = line[0]
		pos1 = int(line[1])
		pos2 = int(line[2])
		sequence=["N"*(pos2-pos1) for i in sample_list]
		SEQ["_".join(line)]=sequence
		locus.append([chr,int(pos1)+1,int(pos2),"_".join(line)]) 
	return locus,SEQ

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
	for seq in record:
		for j in range(len(seq)-1):
			if "C" in annotation[seq[j]] and "G" in annotation[seq[j+1]]:
				pos.append(j)
				pos.append(j+1)
	for i in range(len(record)):
		old_seq=record[i]
		seq = ""
		for j in range(len(old_seq)):
			if j in pos:
				seq+="N"
			else:
				seq +=old_seq[j]
		record[i]=seq
	return record

def build_record(SEQ,length):
	record={}
	for sample in SEQ.keys():
		try:
			record[sample]=SEQ[sample]
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
	f=open('Gphocs_input_canfam3.1_all_sample_GQ30_BQ5_MQ15_DP5_unrecal_aDNA.txt',"w")
	remove_locus=[]
	print >>f,"26250"
	f.write("\n")
#	f_matrix=open("dist_matrix_canfam3_MQ0_DP0_GQ20.txt","w")
	f_matrix=open("dist_matrix_canfam3.1_all_sample_GQ30_BQ5_MQ15_DP5_unrecal_aDNA.txt","w")
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

def calc_distance_for_each_loci(locus,SEQ,sample_name,index):
	matrix=np.zeros((len(sample_name),len(sample_name)))
	size=np.zeros((len(sample_name),len(sample_name)))
	status=1
	pos1=int(locus[index][1])
	pos2=int(locus[index][2])
	tag=locus[index][3]
	chr=locus[index][0]
	ref_seq=ref.fetch(chr,pos1-1,pos2)
	ref_seq=ref_seq.upper()
	record=SEQ[tag]
	record=remove_CG(record)
#	record=remove_damage(record,ref_seq)
	for j in range(len(sample_name)):
		for k in range(j+1,len(sample_name)):
			dist,length=comparison(record[j],record[k])
			if length==0:
				status=0
			matrix[j][k]=dist
			size[j][k]=length
			matrix[k][j]=dist
			size[k][j]=length
	return (matrix,size,status)


def calc_distance_for_each_loci_wrapper(args):
	args2 = args[0] + (args[1],)
	stats=calc_distance_for_each_loci(*args2)
	return stats

def calc_NJ(dist,size,status,sample_name):
	f_matrix=open("dist_matrix_canfam3.1_all_rm_damage.txt","w")
	matrix=np.zeros((len(sample_name),len(sample_name)))
	length=np.zeros((len(sample_name),len(sample_name)))
	for i in range(len(status)):
		if status[i]==1:
			matrix+=dist[i]
			length+=size[i]
	print >>f_matrix,"dist_matrix:"
	print >>f_matrix,"\t".join(sample_name)
	for i in range(len(sample_name)):
		print >>f_matrix,"\t".join(map(str,matrix[i]))
	print >>f_matrix,"length_matrix:"
	print >>f_matrix,"\t".join(sample_name)
	for i in range(len(sample_name)):
		print >>f_matrix,"\t".join(map(str,length[i]))
	print >>f_matrix,"distance_matrix"
	print >>f_matrix,"\t".join(sample_name)
	for i in range(len(sample_name)):
		print >>f_matrix,"\t".join(map(str,[x/y for x,y in zip(matrix[i],length[i])]))

def remove_damage(record,ref_seq):
	annotation={'R':'GA','Y':'CT','K':'TG','M':'CA','S':'CG','W':'TA','A':'AA','T':'TT','C':'CC','G':'GG','N':'NN'}
	pos = []
	for seq in record:
		for j in range(len(seq)):
			if (ref_seq[j]=="C" and "T" in annotation[seq[j]]) or (ref_seq[j]=="T" and "C" in annotation[seq[j]]):
				pos.append(j)
			elif (ref_seq[j]=="G" and "A" in annotation[seq[j]]) or (ref_seq[j]=="A" and "G" in annotation[seq[j]]):
				pos.append(j)
	if len(pos)!=0:
		print pos
	for i in range(len(record)):
		old_seq=record[i]
		seq = ""
		for j in range(len(old_seq)):
			if j in pos:
				seq+="N"
			else:
				seq +=old_seq[j]
		record[i]=seq
	return record

def aln_to_fasta_bootstrap(dist,size,status,sample_name):
	Nbootstrap=100
	values = np.array(status)
	searchval = 1
	usable=np.where(values == searchval)[0]
	for index in range(Nbootstrap):
		print 'bootstrap:',index
		f_matrix=open("dist_matrix_canfam3.1_bootstrap_%s_rm_damage.txt" %(index),"w")
		choose=np.random.choice(usable,size=len(usable))
		matrix=np.zeros((len(sample_name),len(sample_name)))
		length=np.zeros((len(sample_name),len(sample_name)))
		for i in choose:
			matrix+=dist[i]
			length+=size[i]
		print >>f_matrix,"dist_matrix:"
		print >>f_matrix,"\t".join(sample_name)
		for i in range(len(sample_name)):
			print >>f_matrix,"\t".join(map(str,matrix[i]))
		print >>f_matrix,"length_matrix:"
		print >>f_matrix,"\t".join(sample_name)
		for i in range(len(sample_name)):
			print >>f_matrix,"\t".join(map(str,length[i]))
		print >>f_matrix,"distance_matrix"
		print >>f_matrix,"\t".join(sample_name)
		for i in range(len(sample_name)):
			print >>f_matrix,"\t".join(map(str,[x/y for x,y in zip(matrix[i],length[i])]))

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='vcf to fa')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--vcf_dir",default='/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/GPhoCS/input_3.1',dest='vcf_dir',help="vcf dir")
	args = parser.parse_args()

	print 'Start'
	# get sample name
	sample_list=get_sample_name('%s/sampleList.txt' %(args.vcf_dir))
	print 'sample list:',sample_list
	# Initialize SEQ and locus, SEQ[tag]=['','',''], a list of sequence in the order of sample_list
	SEQ={}  # SEQ a dictionary, key is the tag of loci, value is ["N"*length for each sample in sample_list]
	locus_file = "/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/neutralLoci-geneFlank10k-sep30k-3.1.bed"
	locus,SEQ=read_locus_file(locus_file,SEQ,sample_list)  

#	print 12
	ref_file='/home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1-withUn/canFam3.1.withUn.fa'
	ref=pysam.FastaFile(ref_file)

	for i in range(len(sample_list)):
		name=sample_list[i]
		if name in ['Kirshbaum','Herxheim']:
			file='%s/%s_unrecal.GhoCSseq_aDNA_GQ30_BQ5_MQ15_DP5' %(args.vcf_dir,name)
		else:
			file='%s/%s.GhoCSseq_GQ30_BQ5_MQ15_DP5' %(args.vcf_dir,name)
		SEQ,het,Length=read_GhoCSseq_GQ_20(file,locus,SEQ,i)
		print name,het,Length,float(het)/Length

#	m=calc_distance_for_each_loci(locus,SEQ,sample_list,0)
#	print m[0],m[1],m[2]

	pool = multiprocessing.Pool(processes=20)
	result1=pool.map(calc_distance_for_each_loci_wrapper,itertools.izip(itertools.repeat((locus,SEQ,sample_list)),range(len(locus))))
	pool.close()
	pool.join()

	dist=[i[0] for i in result1]
	size=[i[1] for i in result1]
	status=[i[2] for i in result1]
	print sum(status)
	calc_NJ(dist,size,status,sample_list)
	print 'finish NJ tree'
	aln_to_fasta_bootstrap(dist,size,status,sample_list)

