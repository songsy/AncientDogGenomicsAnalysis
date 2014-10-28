#!/usr/bin/env python
# python 05.locus-to-muscle.py
# Shiya Song
# 7rd June 2013

import glob
import os
import sys

def import_sample_name(a):
	sample_name=a+["chimp"] 
	return sample_name

def read_locus_file(file):              # locus interval is bed format, zero based, transform into 1-based
	locus={}
	dict={}
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
		dict[info]=line[3]	#map chimp to human
#	print len(locus)
	return locus,dict

def blat_to_strand(file,dict,locus):
	f=open(file,"r")
	strand={}
	first = True
	last_remove=""
	for line in f:
		if line[0].isdigit() is False:
			continue
		line=line.strip().split("\t")
		name=line[9]
		chr = line[13]
		pos1 = int(line[15])
		pos2 = int(line[16])
		if first is True:
			tag = dict[name]       # human coordinate
			chrom=tag.split(":")[0]
			p1=int(tag.split(":")[1].split("-")[0])
			p2=int(tag.split(":")[1].split("-")[1])
			if chrom==chr and (abs(p1-pos1)<=100 or abs(p2-pos2)<=100):
				strand[name]=[int(line[0]),line[8],chr,pos1,pos2]
			first = False
			last_name=name
			continue
		if last_name==name:
			if chrom==chr and (abs(p1-pos1)<=100 or abs(p2-pos2)<=100):
				if strand.has_key(name) is False:
					strand[name]=[int(line[0]),line[8],chr,pos1,pos2]
				elif int(line[0])>strand[name][0]:
#					print "overlap",last_name
					strand[name]=[int(line[0]),line[8],chr,pos1,pos2]
			elif int(line[0])>900:
				print "multi hit", last_name
				if name!=last_remove:
					locus.remove([chrom,p1,p2,name])
					last_remove= name
		else:
			if strand.has_key(last_name) is False:
				print "no hit",last_name
				if [chrom,p1,p2,last_name] in locus:
					locus.remove([chrom,p1,p2,last_name])
					print "removed"
			tag = dict[name]       # human coordinate
			chrom=tag.split(":")[0]
			p1=int(tag.split(":")[1].split("-")[0])
			p2=int(tag.split(":")[1].split("-")[1])
			if chrom==chr and (abs(p1-pos1)<=100 or abs(p2-pos2)<=100):
				strand[name]=[int(line[0]),line[8],chr,pos1,pos2]
			elif int(line[0])>900:
				print "multi hit", last_name
				if name!=last_remove:
					locus.remove([chrom,p1,p2,name])
					last_remove= name
			last_name=name
	if strand.has_key(last_name) is False:
		print "no hit",last_name
		locus.remove([chrom,p1,p2,last_name])
	return strand,locus
		
def get_coordinate(sample_name):
	dict={}
	i="Victoria"
	file = i + "_within-locus.fa"
	with open(file,"r") as f1:
		for line in f1:
			if line.startswith('>'):
				tag=line.strip()[1:].split(" ")[0].replace(".",":")
				tag_hg = line.strip()[1:].split(" ")[1]
				dict[tag_hg]=tag
	return dict

def print_locus(file,locus):
	f=open(file,"w")
	for i in locus:
		print >>f,"\t".join(map(str,i))
	
def fasta_to_muscle(sample_name,strand,dict,muscle_dir):
	for i in sample_name:
		file = i + "_within-locus.fa"
#		file = i + "_test.fa"
		with open(file,"r") as f1:
			for line in f1:
				if line.startswith('>'):
					tag=line.strip()[1:].split(" ")[0].replace(".",":")
					tag_hg = line.strip()[1:].split(" ")[1]
					length=line.strip()[1:].split(" ")[2]
					f=open(muscle_dir+"/"+tag+".fa","a")
					print >>f,">%s:%s %s" %(i,tag,length)
				else:
					print >>f,line.strip()
					f.close()
		f1.close()
	seq=""
	with open("locus-hg19-new.fa","r") as f2:
		for line in f2:
			if line.startswith('>'):
				if seq!="":
					if strand.has_key(tag_hg) and strand[tag_hg][1]=="-":
						seq=reverse_comp(seq)
					print >>f, seq
					f.close()
					cmd = "muscle -in %s/%s.fa -out %s/%s.aln" %(muscle_dir,tag,muscle_dir,tag)
					os.popen(cmd)
				seq = ""
				tag_hg=line.strip()[1:]
				tag=dict[tag_hg]
				f=open(muscle_dir+"/"+ tag+".fa","a")
				print >>f,">human:%s" %(tag_hg)
			else:
				seq+= line.strip()	
	f2.close()

def reverse_comp(seq):
	reverse={"A":"T","T":"A","C":"G","G":"C"}
	new_seq=""
	for i in seq[::-1]:
		new_seq+=reverse[i.upper()]
	return new_seq
	
def read_in_seq(aln_file):
	record={}
	f=open(aln_file,"r")
	for line in f:
		line=line.strip()
		if line[0]==">":
			sample=line[1:].split(" ")[0].split(":")[0]
			record[sample]=""
		else:
			record[sample]+=line
	return record

def remove_gap(record):
	for i in record:
		record[i]=record[i].replace("-","N")
	return record

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
		length=len(record[i])
	return record,length

def delta_function(a,b):
	if a==b:
		return 1
	else:
		return 0

def calc_heterozygosity(record):
	num = {}
	for i in record:
		num[i] = float(0)
		for j in range(len(record[i])):
			if record[i][j] in ["R","Y","K","M","S","W"]:
				num[i] += 1
	return num
				
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
			
def aln_to_fasta(outfile,sample_name,locus,muscle_dir):
	f=open(outfile,"w")
	remove_locus=[]
	print >>f,"26250"
	f.write("\n")
#	f_heter = open("locus_heterozygosity_2E4W1R","w")
#	print >>f_heter,"chrom\tpos1\tpos2\thet_Victoria\thet_9732\thet_X00108\thet_KB3784\tdiv_Victoria\tdiv_9732\tdiv_X00108\tdiv_KB3784"
	f_matrix=open("dist_matrix_" + outfile.replace("_neutral.txt",""),"w")
#	sample_name.append("gorGor3")
	sample_name.append("human")
	matrix=[]              # set up the matrix to record the distance
	size = []				# record the total length of each comparison
	heter={}
	total_length = 0
	div = {}
	for i in range(len(sample_name)):
		matrix.append([])
		size.append([])
		heter[sample_name[i]]=0
		for j in range(len(sample_name)):
			if i==j:
				matrix[i].append(0)
				size[i].append(1)
			else:
				matrix[i].append(float(0))
				size[i].append(0)
	for i in locus:
		zero = False
		if i[0]=="chrX" or i[3][:4]=="chrX":
			print i
			continue
		aln_file = "%s/%s:%s-%s.aln" %(muscle_dir,i[0],i[1],i[2])
		record=read_in_seq(aln_file)
		record=remove_gap(record)
		record,seq_length=remove_CG(record)
		seq_length=len(record["Victoria"])
		for j in range(len(sample_name)):
			for k in range(j+1,len(sample_name)):
				dist,length=comparison(record[sample_name[j]],record[sample_name[k]])
				if length==0:
					zero = True
					print "remove locus",i
				if sample_name[j]=="Victoria" and sample_name[k]=="human" and length!=0:
					print "locus dist",i,dist,length,dist/length
				if sample_name[k]=="human" and sample_name[j]!="human" and length!=0:
					div[sample_name[j]]=dist/length
				matrix[j][k]+=dist
				size[j][k]+=length
				matrix[k][j]+=dist
				size[k][j]+=length
		if zero is False:
			print >>f,"%s.%s-%s\t%i\t%s" %(i[0],i[1],i[2],len(sample_name),seq_length)
			for j in sample_name:
				print >>f,"%s\t%s" %(j,record[j])
			f.write("\n")
			num = calc_heterozygosity(record)
			print num
#			print "\t".join(map(str,[num[x]/seq_length for x in num]))
#			print >>f_heter, "%s\t%s\t%s\t" %(i[0],i[1],i[2]),"\t".join(map(str,[num[m]/seq_length for m in sample_name]))
			for i in sample_name:
				heter[i]+= num[i]
			total_length+= seq_length
		else:
			remove_locus.append(i)
		print aln_file,"finished"
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
	print >>f_matrix, "Heterozygosity"
	print >>f_matrix,"\t".join(map(str,[heter[m] for m in sample_name])),"\t".join(map(str,[heter[m]/total_length for m in sample_name]))
	return remove_locus

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--sample", dest='sample',nargs='+',help="sample name")
	parser.add_argument("--muscle_dir",dest='muscle_dir',help="muscle dir")
	parser.add_argument("--outfile", dest='outfile',help="output file")
	parser.add_argument("--locusfile", dest='locus_file',help="locus bed file")

	if os.path.isdir(args.muscle_dir) is False:
		os.popen("mkdir "+muscle_dir)

	sample_name=import_sample_name(a)
	print sample_name
	locus_file = "locus-chimp.bed"
	locus,dict=read_locus_file(locus_file)
#	dict = get_coordinate(sample_name)
	strand,locus=blat_to_strand("locus_chimp_hg19.psl",dict,locus)
"""
file = "locus_final.bed"
print_locus(file,locus)
for i in strand:
	print i,"\t","\t".join(map(str,strand[i]))
for i in strand:
	if strand[i][1]=="-":
		print "-",i,strand[i]
	if strand[i][0]<=900:
		print i,strand[i]
"""	
	fasta_to_muscle(sample_name,strand,dict,muscle_dir)
	remove_locus=aln_to_fasta(args.outfile,sample_name,locus,muscle_dir)

	for i in remove_locus:
		locus.remove(i)
	print_locus(args.locus_file,locus)

						
						
					
			