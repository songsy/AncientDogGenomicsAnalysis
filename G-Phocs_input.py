import os
import glob

BSNP_file = "H7C1_neutral.txt"
sample_name=["NA12878","NA21302","NA19240","HG02799","HG03108","HG03428","HGDP01029"]
group=[]
#group.append(["NA12878","NA21302","NA19240","HGDP01029"])
#group.append(["NA12878.1","NA12878.2","NA21302.1","NA21302.2","NA19240.1","NA19240.2","HGDP01029"])
#group.append(["NA12878","NA19240","HGDP01029"])
group.append(["NA12878","NA21302","HG02799","HGDP01029"])
group.append(["NA12878","NA21302","HG03108","HGDP01029"])
group.append(["NA12878","NA21302","HG03428","HGDP01029"])
for item in group:
	f = open(BSNP_file,"r")
	print item
	if item[2]=='HG02799':
		file_name = "Eu_MKK_GWD_San"
	elif item[2]=='HG03108':
		file_name = "Eu_MKK_ESN_San"
	elif item[2]=='HG03428':
		file_name = "Eu_MKK_MSL_San"
	print file_name
	f_out = open(file_name+"/"+file_name+"_neutral.txt","w")
	for line in f:
		line = line.strip()
		if line=="":
			f_out.write("\n")
		elif line[0].isdigit() is True:
			print >>f_out,line
		elif line[:3]=="chr":
			col = line.split("\t")
			print >>f_out,"%s\t%i\t%s" %(col[0],len(item)+1,col[2])
		else:
			col = line.split("\t")
			if col[0] in item or col[0]=="chimp":
				print >>f_out,line
