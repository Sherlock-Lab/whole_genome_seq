import sys


#sys.argv[1] = window size
#sys.argv[2] = index tag
# first need to create pile-up file with samtools samtools mpileup -s AAGAGGCA-CTAAGCCT.dedup.realigned.bam > mpileup.txt

#converts ascii quality to numeric and gets average for given sequence
def get_avgmp(mapq):
	sum_mq = 0
	for letter in mapq:
		sum_mq = sum_mq + ord(letter)
	avg_mq = float(sum_mq)/float(len(mapq))-33.0
	return avg_mq

index_tag = sys.argv[2]
input_path = (index_tag+"_data/"+index_tag+"_mpileup.txt")
output_path = (index_tag+"_data/"+index_tag+"_cov1kb.csv")

#print input_path
#print output_path

input = open(input_path,'r')
output = open(output_path,'w')
output.write("chrom,win.start,reads.gc,reads.mapq,counts\n")

start = 0
window = int(sys.argv[1])	
end = start + window

chrom_convert = {"I":"chr1","II":"chr2","III":"chr3","IV":"chr4","V":"chr5","VI":"chr6","VII":"chr7","VIII":"chr8","IX":"chr9","X":"chr10","XI":"chr11","XII":"chr12","XIII":"chr13","XIV":"chr14","XV":"chr15","XVI":"chr16","Mito":"mt"}
curr_chr = "chr1"

#get beginning parameters
line = input.readline()
line = line.strip().split()
chr = chrom_convert[line[0]]
pos = int(line[1])
cov_list = [int(line[3])]
mapq_list = [get_avgmp(line[6])]

print (chr)
print (pos)
print (cov_list)
print (mapq_list)

#loop through lines in file
for line in input:
	
	#get information from line
	line = line.strip().split()
	chr = chrom_convert[line[0]]
	pos = int(line[1])
	cov = int(line[3])
	if len(line)==7:
		mp = get_avgmp(line[6])
	
	#if we aren't past the end of the window, and aren't on a new chromosome: add to lists
	if pos < end and curr_chr == chr:
		cov_list.append(cov)
		if len(line)==7:
			mapq_list.append(mp)
	
	#if we are in the new window, or on a new chromosome
	else:
		#output stats for previous window
		avg_window_cov = sum(cov_list)/len(cov_list)
		avg_window_mp = sum(mapq_list)/len(mapq_list)
		#if pos < end:
			#end = pos
		#print(curr_chr+","+str(start)+",0.38,"+str(avg_window_mp)+","+str(avg_window_cov))
		output.write(curr_chr+","+str(start)+",0.38,"+str(avg_window_mp)+","+str(avg_window_cov)+"\n")
		#print len(cov_list)
		#print len(mapq_list)
		

		#re-set parameters for new window
		#if we are on a new chromosome, reset current end of window and chromosome accordingly
		if curr_chr != chr:
			start = 0
			end = window
			curr_chr = chr
		
		#if we are continuing on same chromosome, reset current end of window accordningly
		else:
			start = end
			end = end + window
		
		#reset the two lists for the new window with the information for this position
		cov_list = [cov]
		if len(line)==7:
			mapq_list = [mp]
		else:
			mapq_list = [0]
		
		
			
	
input.close()
output.close()