#!/usr/bin/env python
# python 01.bcftools_cmds.py
# Shiya Song
# 26th October 2015
import sys
import glob
import os
import argparse

def bcftools_cmds():
	cmds = "/home/jmkidd/kidd-lab/progs/samtools-0.1.16/samtools pileup -csf fasta/CHROM BAM|"
	cmds += "BSNP -i - -o file1 -on file2 -os file3 -ka 1.0 -p0 0.001 -mq 1 -mp 0 -th 0.85 -tm 0 -si 1 -ns 1 -ig 1 -v 1 -pb 1 -st 0"
	for chr in chrom_list:
		cmds1 = cmds
		cmds1 = cmds1.replace("BAM","../BAM/%s_%s.bam" %(args.name,chr))
		cmds1 = cmds1.replace("CHROM",chr+".fa")
		cmds1 = cmds1.replace("file1","%s_%s.ALL.SNP_out" %(args.name,chr))
		cmds1 = cmds1.replace("file2","%s_%s.ALL.AUX_out" %(args.name,chr))
		cmds1 = cmds1.replace("file3","%s_%s.SUM_out" %(args.name,chr))
		print cmds1

def BSNP_cmds():
	cmds = "/home/jmkidd/kidd-lab/progs/samtools-0.1.16/samtools pileup -csf /home/jmkidd/kidd-lab/genomes/canFam3/canFam3-cat/canFam3.fa BAM|"
	cmds += "BSNP -i - -o file1 -on file2 -os file3 -ka 1.0 -p0 0.001 -mq 1 -mp 0 -th 0.85 -tm 0 -si 1 -ns 1 -ig 1 -v 1 -pb 1 -st 0"
	for bam in ['Basenji','Dingo','IsraeliWolf','CroatianWolf','ChineseWolf','GoldenJackal','CTC']:
		cmds1 = cmds
		cmds1 = cmds1.replace("BAM","/mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/Freedman_bams/%s_GPhoCs_merged.bam" %(bam))
		cmds1 = cmds1.replace("file1","%s.ALL.SNP_out" %(bam))
		cmds1 = cmds1.replace("file2","%s.ALL.AUX_out" %(bam))
		cmds1 = cmds1.replace("file3","%s.SUM_out" %(bam))
		print cmds1

def GATK_cmds():
	for bam in ['Basenji','Dingo','IsraeliWolf','CroatianWolf','ChineseWolf','GoldenJackal','CTC','HxH']:
		cmds = "java -Xmx4g -jar /home/jmkidd/kidd-lab/progs/GATK-3.4/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/jmkidd/kidd-lab/genomes/canFam3/canFam3-cat/canFam3.fa --emitRefConfidence GVCF -variant_index_type LINEAR -variant_index_parameter 128000 "
		cmds+= "-o %s_GPhoCS_gvcf.gz -I /mnt/EXT/Kidd-scratch/shiya-projects/ancient_dog/Krishna/Freedman_bams/%s_GPhoCs_merged.bam" %(bam,bam)
		print cmds

def get_chrom_list():
	chrom_list=[]
	f=open('/home/jmkidd/kidd-lab/genomes/canFam3.1/chromOrder.txt','r')
	for line in f:
		chrom_fa = line.rstrip()
		chrom=chrom_fa.replace('.fa','')
		chrom_list.append(chrom)
		if chrom_fa=='chrM.fa':
			break
	return chrom_list

def split_bam():
	for chr in chrom_list:
		cmds='samtools view -b %s %s > %s_%s.bam' %(args.bam,chr,args.name,chr)
		print cmds

def BSNP_to_vcf_cmds():
	for chr in chrom_list:
		cmds='python /home/songsy/script/G-phocs/02.BSNP_to_vcf.py %s %s' %(args.name,chr)
		print cmds

def calc_het_cmds():
	for chr in chrom_list:
		cmds='python /home/songsy/script/G-phocs/04.calc_heterozygosity.py --sample %s --chr %s' %(args.name,chr)
		print cmds

def liftover_cmds():
	for name in ['Basenji','Dingo','IsraeliWolf','CroatianWolf','ChineseWolf','GoldenJackal','CTC']:
		cmds='java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g /home/jmkidd/kidd-lab/progs/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar '
		cmds+='-T LiftoverVariants -R /home/jmkidd/kidd-lab/genomes/canFam3/canFam3-cat/canFam3.fa -V %s_BSNP.vcf.gz -chain /home/jmkidd/kidd-lab/genomes/lift-over/canFam3.0TocanFam3.1.over.chain ' %(name)
   		cmds+='-dict /home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1-cat/canFam3.1.fa.dict -o %s_3.1_BSNP.vcf ' %(name)
  		print cmds
  
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write bcftools commands')
	parser.add_argument("--bam", dest='bam',help="bam file")
	parser.add_argument("--name", dest='name',help="name prSefix")  
	args = parser.parse_args()

	chrom_list=get_chrom_list()
#	split_bam()
#	bcftools_cmds()
#	BSNP_to_vcf_cmds()
#	calc_het_cmds()
#	BSNP_cmds()
#	liftover_cmds()
#	GATK_cmds()
