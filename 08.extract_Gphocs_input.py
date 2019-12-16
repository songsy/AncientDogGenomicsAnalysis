import os

file = "Gphocs_input_canfam3_GQ30_BQ5_MQ15_DP5_unrecal_aDNA_varyGQ.txt"
#sample_name=["Basenji","Boxer","ChineseWolf","CroatianWolf","Dingo","GoldenJackal","Herxheim_30","IsraeliWolf","Kirshbaum_30"]
sample_name=["Basenji","Kirshbaum_30","ChineseWolf","CroatianWolf","Dingo","GoldenJackal","IsraeliWolf"]
group = [sample_name]

for item in group:
	f = open(file,"r")
	print item
	f_out = open("Gphocs_input_canfam3_GQ30_CTC_neutral.txt","w")
	for line in f:
		line = line.strip()
		if line=="":
			f_out.write("\n")
		elif line[0].isdigit() is True:
			print >>f_out,line
		elif line[0]=="c":
			col = line.split("\t")
			print >>f_out,"%s\t%i\t%s" %(col[0],len(item),col[2])
		else:
			col = line.split("\t")
			if col[0] in item:
				line=line.replace("_30","")
				print >>f_out,line