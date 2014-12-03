import argparse

def delta_function(a,b):
	if a==b:
		return 1
	else:
		return 0
				
def comparison(seq1,seq2):
	assert len(seq1)==len(seq2)
	dist = 0
	length = 0
	a=seq1
	b=seq2
	dist+= 1- float(0.5)*max(delta_function(a[0],b[0])+delta_function(a[1],b[1]),delta_function(a[0],b[1])+delta_function(a[1],b[0]))
	return dist

def calc_distance(file,distance_matrix,total_length):
	f=open(file,'r')
	for line in f:
		line=line.strip().split('\t')
		total_length+=int(line[2])
		seq = line[3].split(',')[0]
		for i in range(len(distance_matrix)):
			for j in range(i,len(distance_matrix)):
				seq1=seq[i*2:i*2+2]
				seq2=seq[j*2:j*2+2]
				distance_matrix[i][j]+=comparison(seq1,seq2)
				distance_matrix[j][i]+=comparison(seq1,seq2)
	return distance_matrix,total_length

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='calculate distance matrix')
	parser.add_argument("--sample", dest='sample',nargs='+',help="sample name")
	args = parser.parse_args()		
	f_matrix=open('distance_matrix_wgs.txt','w')
	sample_name=[]
	for i in range(len(args.sample)):
		sample_name.append(args.sample[i])
#		sample_name.append(args.sample[i]+'.1')
#		sample_name.append(args.sample[i]+'.2')
	matrix=[]              # set up the matrix to record the distance
	total_length = 0
	for i in range(len(sample_name)):
		matrix.append([])
		for j in range(len(sample_name)):
			matrix[i].append(float(0))
	for i in range(1,23):
		file = '%s_all_pop.txt' %(i)
		matrix,total_length=calc_distance(file,matrix,total_length)
		print i
	print >>f_matrix,"dist_matrix:"
	print >>f_matrix,"\t".join(sample_name)
	for i in range(len(sample_name)):
		print >>f_matrix,"\t".join(map(str,matrix[i]))
	print >>f_matrix,"total_length:",total_length
	print >>f_matrix,"distance_matrix"
	print >>f_matrix,"\t".join(sample_name)
	for i in range(len(sample_name)):
		print >>f_matrix,"\t".join(map(str,[x/total_length for x in matrix[i]]))

