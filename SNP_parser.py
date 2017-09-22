import os.path
import sys

#define working directory from which this program is being run
wkdir = sys.argv[1]
#wkdir = '~/Desktop/SherlockLab/revision_work/sequence_analysis/'

index_key = sys.argv[2]

vcfFileSuffix = 'gatk.merged.harsh.filtered.ann.final.vcf.txt'
outFileSuffix = '_list_filtered.txt'
outFileStringent = '_list_filtered_stringent.txt'

all_index_tags_dic  = {}
index_file = open (index_key, 'r')
sample_count = -1
for line in index_file:
	(strain, I7, seq7, I5, seq5) = line.rstrip('\n').split('\t') # split the line
	sample_count += 1
	all_index_tags_dic[strain] = seq7 + '-' + seq5 + '_S' + str(sample_count)
index_file.close()

os.system('mkdir mutations')

all_mut_info_file = open('mutations/all_mut_info.txt', 'w')
mutation_dic = {}
mut_info_dic = {}
mut_id_info = {}
for strain in all_index_tags_dic:
	index_tag = strain + '-' + all_index_tags_dic[strain]
	print (index_tag)	

	#define the folder for data output for this clone
	data_dir = wkdir + index_tag + "_data/"
	
	#define files names of this particular index tag
	vcfDir = data_dir + index_tag + vcfFileSuffix
	
	vcfFile = open(vcfDir, 'r')

	for line in vcfFile:
		#-F CHROM -F POS -F REF -F ALT -F FILTER -F SNPEFF_GENE_NAME -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_EFFECT -F QUAL -F DP -GF AD -GF GQ -GF PL
		(ref_chr, pos, ref_base, mut, filter, gene, func_class, codon_chg, aa_chg, effect, qual, filtered_depth, allele_depth1, GQ1, PL1, allele_depth2, GQ2, PL2) = line.rstrip().split('\t')
	
		mut_ID = ref_chr + '_' + pos + '_' + ref_base + '_' + mut
		if allele_depth1 == 'NA':
			allele_depth = allele_depth2
			GQ = GQ2
			PL = PL2
		else:
			allele_depth = allele_depth1
			GQ = GQ1
			PL = PL1



		mut_info = str.join('\t', (ref_chr, pos, ref_base, mut, filter, gene, func_class, codon_chg, aa_chg, effect, qual, filtered_depth, allele_depth, GQ, PL))
		mutation_strain_id = strain + '_' + mut_ID
		mut_id_info[mutation_strain_id] = mut_info
		info = strain + '\t' + mut_info + '\n'
		all_mut_info_file.write(info)
			
		if mut_ID in mutation_dic:
	
			mutation_dic[mut_ID].append(strain)
	
		else:
			mutation_dic[mut_ID] = [strain]
			mut_info_dic[mut_ID] = mut_info
			
			
	vcfFile.close()
			
all_mut_info_file.close()	


threshold_100 = mysamples
threshold_80 = mysamples * 0.8
threshold_50 = mysamples * 0.5
threshold_20 = mysamples * 0.2
threshold_10 = mysamples * 0.1
threshold_1 = 1


threshold_50_lst = []

outFileDic = {'all_file' : open('mutations/all_SNP.txt', 'w'), 'threshold_100_file': open('mutations/SNP_100.txt', 'w'), 'threshold_80_file': open('mutations/SNP_more_80.txt', 'w'), 'threshold_50_file': open('mutations/SNP_more_50.txt', 'w'), 'threshold_20_file': open('mutations/SNP_more_20.txt', 'w'), 'threshold_10_file': open('mutations/SNP_more_10.txt', 'w'), 'threshold_10_1_file': open('mutations/SNP_10_1.txt', 'w'), 'threshold_1_file': open('mutations/SNP_1.txt', 'w')}

print("Making all threshold files")
for mut in mutation_dic:
	mut_occurance = len(mutation_dic[mut])
	
	x = mut + '\t' + str(mut_occurance) + '\n'
	outfile = outFileDic['all_file']
	outfile.write(x)

	if mut_occurance == threshold_100:
		w100 = mut + '\t' + str(mut_occurance) + '\t' + mut_info_dic[mut] + '\n'
		outfile = outFileDic['threshold_100_file']
		outfile.write(w100)
		threshold_50_lst.append(mut)
	
	elif mut_occurance >= threshold_80 and mut_occurance < threshold_100 :
		w80 = mut + '\t' + str(mut_occurance) + '\t' + mut_info_dic[mut] + '\n'
		outfile = outFileDic['threshold_80_file']
		outfile.write(w80)
		threshold_50_lst.append(mut)
		
	elif mut_occurance >= threshold_50 and mut_occurance < threshold_80:
		w50 = mut + '\t' + str(mut_occurance) + '\t' + mut_info_dic[mut] + '\n'
		outfile = outFileDic['threshold_50_file']
		outfile.write(w50)
		threshold_50_lst.append(mut)
		
	elif mut_occurance >= threshold_20 and mut_occurance < threshold_50:
		for strain_ref in mutation_dic[mut]:
			mutation_ref_id = strain_ref + '_' + mut
			mut_id_info[mutation_ref_id]

			w20 = mut + '\t' + strain_ref + '\t' + mut_id_info[mutation_ref_id]  + '\t' + str(mutation_dic[mut]) + '\n'
			outfile = outFileDic['threshold_20_file']
			outfile.write(w20)
		
	elif mut_occurance >= threshold_10 and mut_occurance < threshold_20:
		for strain_ref in mutation_dic[mut]:
			mutation_ref_id = strain_ref + '_' + mut
			mut_id_info[mutation_ref_id]
			
			w10 = mut + '\t' + strain_ref + '\t' + mut_id_info[mutation_ref_id]  + '\t' + str(mutation_dic[mut]) + '\n'
			outfile = outFileDic['threshold_10_file']
			outfile.write(w10)
		
	elif mut_occurance > 1 and mut_occurance < threshold_10:
		for strain_ref in mutation_dic[mut]:
			mutation_ref_id = strain_ref + '_' + mut
			mut_id_info[mutation_ref_id]
			
			w10_1 = mut + '\t' + strain_ref + '\t' + mut_id_info[mutation_ref_id]  + '\t' + str(mutation_dic[mut]) + '\n'
			outfile = outFileDic['threshold_10_1_file']
			outfile.write(w10_1)
		
	elif mut_occurance == 1:
		w1 = mut + '\t' + str(mutation_dic[mut]) + '\t' + mut_info_dic[mut] + '\n'
		outfile = outFileDic['threshold_1_file']
		outfile.write(w1)
		
	else:
		print("Neither of those " + mut, mut_occurance, mut_info_dic[mut])
		
for key_file in outFileDic:
	outFileDic[key_file].close()	
		

os.system('sort mutations/all_SNP.txt -n -r -k2,2 > mutations/all_SNP_sorted.txt')
os.system('rm mutations/all_SNP.txt')


error_region_file = open('Regions_erroneous.txt', 'r')
err_reg_dic = {}
for err in error_region_file:
	
	chr, s_pos, e_pos = err.rstrip('\n').split('\t')
	posistion = s_pos + '_' + e_pos
	if chr not in err_reg_dic:
		err_reg_dic[chr] = [posistion]
	else:
		err_reg_dic[chr].append(posistion)

all_Stringent = wkdir + "mutations/all_SNP_stringent.txt"
stringent_all = open(all_Stringent, 'w')
out_Stringent = wkdir + "mutations/all_out_stringent.txt"
out_stringent = open(out_Stringent, 'w')
print("Removing frequent mutation and making new filtered list of mutation ")
for strain in all_index_tags_dic:
	index_tag = strain + '-' + all_index_tags_dic[strain]
	print (index_tag)	

	#define the folder for data output for this clone
	data_dir = wkdir + index_tag + "_data/"
	
	#define files names of this particular index tag
	vcfDir = data_dir + index_tag + vcfFileSuffix
	outDir = data_dir + index_tag + outFileSuffix
	outStringent = data_dir + index_tag + outFileStringent
	
	vcfFile = open(vcfDir, 'r')
	outFile = open(outDir, 'w')
	stringentFile = open(outStringent, 'w')

	for line in vcfFile:
		#-F CHROM -F POS -F REF -F ALT -F FILTER -F SNPEFF_GENE_NAME -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE  -F SNPEFF_EFFECT -F QUAL -F DP -GF AD -GF GQ -GF PL
		(ref_chr, pos, ref_base, mut, filter, gene, func_class, codon_chg, aa_chg, effect, qual, filtered_depth, allele_depth1, GQ1, PL1, allele_depth2, GQ2, PL2) = line.rstrip().split('\t')
	
		mut_id = ref_chr + '_' + pos + '_' + ref_base + '_' + mut
		if allele_depth1 == 'NA':
			allele_depth = allele_depth2
			GQ = GQ2
			PL = PL2
		else:
			allele_depth = allele_depth1
			GQ = GQ1
			PL = PL1

		if mut_id not in threshold_50_lst :
			kw = str.join('\t', (ref_chr, pos, ref_base, mut, filter, gene, func_class, codon_chg, aa_chg, effect, qual, filtered_depth, allele_depth, GQ, PL))
			outFile.write(kw + '\n')
			
			ref_A = allele_depth.split(',')[0]
			mut_A = allele_depth.split(',')[1]
			threshold_mut = (float(mut_A) + float(ref_A)) * 0.8
			erronous_region = False
			
			if ref_chr in err_reg_dic:
				pos_err_lst = err_reg_dic[ref_chr]

				for pos_err in pos_err_lst:
					pos_err_s, pos_err_e = pos_err.split('_')
					if int(pos) >= int(pos_err_s) and int(pos) <= int(pos_err_e):
						erronous_region = True
						break
			
			if float(qual) > 100 and 5 < int(filtered_depth) and float(mut_A) > threshold_mut and  float(mut_A) > 3 and erronous_region == False:
				sw = str.join('\t', (ref_chr, pos, ref_base, mut, filter, gene, func_class, codon_chg, aa_chg, qual, filtered_depth, allele_depth, GQ, PL))
				stringentFile.write(sw + '\n')
				all_sw = str.join('\t', (index_tag, ref_chr, pos, ref_base, mut, filter, gene, func_class, codon_chg, aa_chg, effect, qual, filtered_depth, allele_depth, GQ, PL))
				stringent_all.write(all_sw + '\n')
			else:
				out_sw = str.join('\t', (index_tag, ref_chr, pos, ref_base, mut, filter, gene, func_class, codon_chg, aa_chg, effect, qual, filtered_depth, allele_depth, GQ, PL))
				out_stringent.write(out_sw + '\n')
	
	vcfFile.close()
	outFile.close()
	stringentFile.close()


stringent_all.close()


os.system('sort mutations/all_SNP_stringent.txt -r -k7,7 -k3,3 > mutations/all_SNP_stringent_gene.txt')
os.system('sort mutations/all_SNP_stringent.txt -r -k1,1 > mutations/all_SNP_stringent_strain.txt')
os.system('rm mutations/all_SNP_stringent.txt')










