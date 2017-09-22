#this should all be run from the directory containing the seq_pipeline.py program
#this directory should also contain a folder of fastq files from the sequencer

#along with the index_key.txt file the cov_plot.R file, coverageLH.py file and
#if Filter_SNP is True the SNP_parser file with the Regions_erroneous.txt

import os

#define which steps need to be carried out
mkOutputDir = True
trim = True
mapping = True
coverage = True
SV = False
variantCall = True
annotation = True
deletionCheck = False # does not work
Filter_SNP = True

#define working directory from which this program is being run
wkdir = 'YOUR_WORKING_DIRECTORY_HERE'

#define directory where all the fastq data files are and file prefix and suffix
dataDir = wkdir+'fastq_files/'

r1FileSuffix = '_L001_R1_001.fastq'
r2FileSuffix = '_L001_R2_001.fastq'


all_index_tags_dic  = {}
index_file = open ('index_key.txt', 'r') # Define your filename that contain the sample
sample_count = 0
for line in index_file:
	(strain, I7, seq7, I5, seq5) = line.rstrip('\n').split('\t') # split the line
	sample_count += 1
	all_index_tags_dic[strain] = seq7 + '-' + seq5 + '_S' + str(sample_count) # Change according to the file name you get from the facility to match the fastq files name.
index_file.close()


beginningTrim = str(10)
qualityCutoff = str(30)
lengthThresh = str(12)
read1Adaptor = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
read2Adaptor = 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'

#define bwa and parameters
RefIndex = '/Volumes/Promise_Pegasus/Lucas/whole_genome_seq/Reference_index/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'
#The reference is already build !!!!!!
# Best to be in the Reference_index directory and build those there.
# bwa index fasta_file.fa
# samtools faidx fasta_file.fa
# java -jar /Volumes/Promise_Pegasus/picard.jar CreateSequenceDictionary R=fasta_file.fa O=path/to/fasta_file.dict

#define picard directory
picard = '/Volumes/Promise_Pegasus/picard.jar'

#define GATK directory
gatk = '/Volumes/Promise_Pegasus/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
#RefIndex = '/Volumes/Promise_Pegasus/Lucas/whole_genome_seq/Reference_gatk/S288C_reference_sequence_R64-2-1_20150113.fsa'

#define snpEff directory
snpEff = '/Volumes/Promise_Pegasus/snpEff/snpEff.jar'
snpEff_config = '/Volumes/Promise_Pegasus/snpEff/snpEff.config'

for strain in all_index_tags_dic:
	index_tag = strain + '-' + all_index_tags_dic[strain] # make sure this match the file name of the fastq files
	print (index_tag)
	#TAG SPECIFIC PARAMETERS
	#define the index tag data to be analyzed
	#index_tag = 'TAAGGCGA-TAGATCGC'

	#define files names of this particular index tag
	fastqR1 = dataDir + index_tag + r1FileSuffix
	fastqR2 = dataDir + index_tag + r2FileSuffix

	#make a folder for data output for this clone
	outputDir = wkdir + index_tag + "_data/"
	outDir = index_tag + "_data/"

	#define read group info needed for correct SAM files
	readGroupInfo = '@RG\\tID:group1\\tSM:'+ index_tag + '\\tPL:illumina\\tLB:12_08_2014\\tPU:'+index_tag

	if mkOutputDir == True:
		os.system("mkdir " + index_tag + "_data")
		#print("mkdir " + index_tag + "_data")

	if trim == True:

		#Trimming reads for quality, adaptor sequence and minimum length
		print("First round of trimming")
		os.system("cutAdapt --minimum-length " + lengthThresh +  " -q " + qualityCutoff + " -u " + beginningTrim + " -a "+ read1Adaptor + " -o "+ outputDir + "tmp.1.fastq -p " + outputDir + "tmp.2.fastq " +  fastqR1 + " " + fastqR2)
		#print("cutAdapt --minimum-length " + lengthThresh +  " -q " + qualityCutoff + " -u " + beginningTrim + " -a "+ read1Adaptor + " -o "+ outputDir + "tmp.1.fastq -p " + outputDir + "tmp.2.fastq " +  fastqR1 + " " + fastqR2)
	
		print("Second round of trimming")
		os.system("cutadapt --minimum-length " + lengthThresh +  " -q " + qualityCutoff + " -u " + beginningTrim + " -a " + read2Adaptor + " -o " + outputDir + index_tag + "_trimmedR2.fastq -p " + outputDir + index_tag + "_trimmedR1.fastq " + outputDir + "tmp.2.fastq " + outputDir + "tmp.1.fastq" )
		#print("cutadapt --minimum-length " + lengthThresh +  " -q " + qualityCutoff + " -u " + beginningTrim + " -a " + read2Adaptor + " -o " + outputDir + index_tag + "_trimmedR2.fastq -p " + outputDir + index_tag + "_trimmedR1.fastq " + outputDir + "tmp.2.fastq " + outputDir + "tmp.1.fastq")
		
		print("Removing temporary files")
		os.system("rm " + outputDir + "tmp.1.fastq")
		os.system("rm " + outputDir + "tmp.2.fastq")
		#print("rm " + outputDir + "tmp.1.fastq")
		#print("rm " + outputDir + "tmp.2.fastq")

		#Format from cutadapt documentation:
		#cutadapt -q 10 -a ADAPTER_FWD --minimum-length 20 -o tmp.1.fastq -p tmp.2.fastq reads.1.fastq reads.2.fastq
		#cutadapt -q 15 -a ADAPTER_REV --minimum-length 20 -o trimmed.2.fastq -p trimmed.1.fastq tmp.2.fastq tmp.1.fastq
	
	if mapping == True:
		# map with BWA -- right now all the parameters are default
		print("Mapping with BWA")
		os.system('bwa mem -M -R "' + readGroupInfo + '" ' + RefIndex + ' ' + outputDir + index_tag + '_trimmedR1.fastq ' + outputDir + index_tag + '_trimmedR2.fastq > ' + outputDir + index_tag + '.sam')
		#print('bwa mem -M -R "' + readGroupInfo + '" ' + RefIndex + ' ' + outputDir + index_tag + '_trimmedR1.fastq ' + outputDir + index_tag + '_trimmedR2.fastq > ' + outputDir + index_tag + '.sam')
		
		#Sort SAM file (saves as bam file)
		print("Sorting SAM file")
		os.system('java -Xmx2g -jar ' + picard + ' SortSam INPUT=' + outputDir + index_tag + '.sam OUTPUT=' + outputDir + index_tag + '.bam SORT_ORDER=coordinate')
		#print('java -Xmx2g -jar ' + picard + ' SortSam INPUT=' + outputDir + index_tag + '.sam OUTPUT=' + outputDir + index_tag + '.bam SORT_ORDER=coordinate')
		

		#Mark duplicates 
		print("Marking duplicates")
		os.system('java -Xmx2g -jar ' + picard + ' MarkDuplicates ASSUME_SORTED=true REMOVE_DUPLICATES=true INPUT= ' + outputDir + index_tag + '.bam OUTPUT= ' + outputDir + index_tag + '.dedup.bam METRICS_FILE= '+ outputDir + index_tag + '.dedup.metrics' )
		#print('java -Xmx2g -jar ' + picard + ' MarkDuplicates ASSUME_SORTED=true REMOVE_DUPLICATES=true INPUT= ' + outputDir + index_tag + '.bam OUTPUT= ' + outputDir + index_tag + '.dedup.bam METRICS_FILE= '+ outputDir + index_tag + '.dedup.metrics' )
		
		#Build BAM index
		print("Building BAM index")
		os.system('java -Xmx2g -jar ' + picard + ' BuildBamIndex INPUT= ' + outputDir + index_tag + '.dedup.bam')
		#print('java -Xmx2g -jar ' + picard + ' BuildBamIndex INPUT= ' + outputDir + index_tag + '.dedup.bam' )
		
		
		#Realignments around indels
		print("Realigning around indels")
		os.system("java -jar -Xmx2g " + gatk + " -T RealignerTargetCreator -R " + RefIndex + " -I " + outputDir + index_tag + '.dedup.bam -o ' + outputDir + index_tag + "_target_intervals.list" )
		#print("java -jar -Xmx2g " + gatk + " -T RealignerTargetCreator -R " + RefIndex + " -I " + outputDir + index_tag + '.dedup.bam -o ' + outputDir + index_tag + "_target_intervals.list")
		os.system("java -jar -Xmx2g " + gatk + " -T IndelRealigner -R " + RefIndex + " -I " + outputDir + index_tag + '.dedup.bam -targetIntervals ' + outputDir + index_tag + "_target_intervals.list -o " + outputDir + index_tag + '.dedup.realigned.bam' )
		#print("java -jar -Xmx2g " + gatk + " -T IndelRealigner -R " + RefIndex + " -I " + outputDir + index_tag + '.dedup.bam -targetIntervals ' + outputDir + index_tag + "_target_intervals.list -o " + outputDir + index_tag + '.dedup.realigned.bam' )
		
		
		#Delete unecessary files
		print("Deleting unneeded files")
		os.system("rm " + outputDir + index_tag + ".sam")
		os.system("rm " + outputDir + index_tag + ".bam")
		os.system("rm " + outputDir + index_tag + ".dedup.bam")
		#os.system("rm " + outputDir + index_tag + ".coverage_stats.sample_interval_statistics")
		#os.system("rm " + outputDir + index_tag + ".coverage_stats.sample_interval_summary")
		os.system("rm " + outputDir + index_tag + "_target_intervals.list")
		os.system("rm " + outputDir + index_tag + ".dedup.bai")
		os.system("rm " + outputDir + index_tag + ".dedup.metrics")
		
		#print("rm " + outputDir + index_tag + ".sam")
		#print("rm " + outputDir + index_tag + ".bam")
		#print("rm " + outputDir + index_tag + ".dedup.bam")
		#print("rm " + outputDir + index_tag + ".coverage_stats.sample_interval_statistics")
		#print("rm " + outputDir + index_tag + ".coverage_stats.sample_interval_summary")
		#print("rm " + outputDir + index_tag + "_target_intervals.list")
		#print("rm " + outputDir + index_tag + ".dedup.bai")
		#print("rm " + outputDir + index_tag + ".dedup.metrics")
		
	if coverage == True:	
		#make pile up file
		print("making pileup")
		#samtools mpileup -s AAGAGGCA-CTAAGCCT.dedup.realigned.bam > mpileup.txt
		os.system("samtools mpileup -s " + outputDir + index_tag + '.dedup.realigned.bam > ' + outputDir + index_tag +"_mpileup.txt")
		#print("samtools mpileup -s " + outputDir + index_tag + '.dedup.realigned.bam > ' + outputDir + index_tag +"_mpileup.txt")


		#make coverage file
		print("making coverage file")
		#python coverage.py index_tag 10000 coverage_10kb.csv
		os.system("python coverageLH.py " + index_tag + " 1000 ")

		#remove pile up file
		print("removing pileup")
		os.system("rm " + outputDir + index_tag +"_mpileup.txt")

		#make plots
		print("plotting")
		#Rscript cov_plot.r 
		print("Rscript cov_plot.r "+ wkdir + ' ' + index_tag + ' 40')

		os.system("Rscript cov_plot.r " + wkdir + ' ' + index_tag + ' 40')
		# here we give to the R script 3 arguments, the working directory the name of the sample and the y axis limit that you can change according to your expected coverage
		
	if SV == True:

		print('Calling SV')

		os.system('speedseq realign '+ RefIndex + ' ' + outputDir + index_tag + '.dedup.realigned.bam')
		#print('speedseq realign '+ RefIndex + ' ' + outputDir + index_tag + '.dedup.realigned.bam')

		os.system('speedseq sv -B ' + index_tag + '.realign.bam -S ' + index_tag + '.realign.splitters.bam -D ' + index_tag + '.realign.discordants.bam -R '+ RefIndex + ' -o '+ index_tag)
		#print('speedseq sv -B ' + index_tag + '.realign.bam -S ' + index_tag + '.realign.splitters.bam -D ' + + index_tag + '.realign.discordants.bam -R '+ RefIndex + ' -o ' + index_tag)

		os.system("java -jar -Xmx2g " + gatk + " -T VariantsToTable -raw -R " + RefIndex + " -F CHROM -F POS -F REF -F ALT -F SVTYPE -F SVLEN -GF SU -GF DP -V " + index_tag + '.vcf -o ' + index_tag + '_SV.txt')

		InFile = open(index_tag + '_SV.txt', 'r')
		OutFile =open(index_tag + 'filtered_SV.txt','w')

		for line in InFile:
			CHROM,POS,REF,ALT,SVTYPE,SVLEN,sample.SU,sample.DP = line.rstrip().split()
			if SVTYPE != "BND" and int(sample.DP) > 5:
				wSV = str.join('\t',(CHROM,POS,REF,ALT,SVTYPE,SVLEN,sample.SU,sample.DP))
				OutFile.write(wSV + '\n')
		InFile.close()
		OutFile.close()

	if variantCall == True:
		print("Calling variants with GATK's Unified Genotyper")
		os.system("java -jar -Xmx2g " + gatk + " -T UnifiedGenotyper -R " + RefIndex + " -I " + outputDir + index_tag + '.dedup.realigned.bam -ploidy 2 --genotype_likelihoods_model BOTH -stand_call_conf 30 -stand_emit_conf 10 -o ' + outputDir + index_tag + '_gatk.vcf')
		#print("java -jar -Xmx2g " + gatk + " -T UnifiedGenotyper -R " + RefIndex + " -I " + outputDir + index_tag + '.dedup.realigned.bam -ploidy 2 --genotype_likelihoods_model BOTH -stand_call_conf 30 -stand_emit_conf 10 -o ' + outputDir + index_tag + '_gatk.vcf')
		
		print("Filtering using Sandeep's threshold")
		
		
		#GATK recommended filtering
		#brute force filtering
		#Not sure what they all mean
		outfile1 = "gatk.snps.vcf"
		outfile2 = "gatk.snps.filtered.vcf"
		outfile3 = "gatk.indels.vcf"
		outfile4 = "gatk.indels.filtered.vcf"
		outfile5 = "gatk.merged.filtered.vcf"
		outfile6 = "gatk.merged.harsh.filtered.vcf"
		
		print("filering SNP")
		os.system("java -Xmx2g -jar  " + gatk + " -T SelectVariants -R " + RefIndex + " --variant " + outputDir + index_tag + '_gatk.vcf -o ' + outputDir + index_tag +  outfile1 + " -selectType SNP --excludeFiltered")
		#print("java -Xmx2g -jar  " + gatk + " -T SelectVariants -R " + RefIndex + " --variant " + outputDir + index_tag + '_gatk.vcf -o ' + outputDir + index_tag +  outfile1 + " -selectType SNP --excludeFiltered")

		os.system("java -Xmx2g -jar  " + gatk + " -T VariantFiltration -R " + RefIndex + " --variant " + outputDir + index_tag + outfile1 + ' -o ' + outputDir + index_tag +  outfile2 + ' --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "snp_filter"')
		#print("java -Xmx2g -jar  " + gatk + " -T VariantFiltration -R " + RefIndex + " --variant " + outputDir + index_tag + outfile1 + ' -o ' + outputDir + index_tag +  outfile2 + ' --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "snp_filter"')

		print("filtering INDELS")
		os.system("java -Xmx2g -jar  " + gatk + " -T SelectVariants -R " + RefIndex + " --variant " + outputDir + index_tag + '_gatk.vcf -o ' + outputDir + index_tag +  outfile3 + " -selectType INDEL --excludeFiltered")
		#print("java -Xmx2g -jar  " + gatk + " -T SelectVariants -R " + RefIndex + " --variant " + outputDir + index_tag + '_gatk.vcf -o ' + outputDir + index_tag +  outfile3 + " -selectType INDEL --excludeFiltered")

		os.system("java -Xmx2g -jar  " + gatk + " -T VariantFiltration -R " + RefIndex + " --variant " + outputDir + index_tag + outfile3 + ' -o ' + outputDir + index_tag +  outfile4 + ' --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "indel_filter" ')
		#print("java -Xmx2g -jar  " + gatk + " -T VariantFiltration -R " + RefIndex + " --variant " + outputDir + index_tag + outfile3 + ' -o ' + outputDir + index_tag +  outfile4 + ' --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "indel_filter" ')
		
		print("merging files")
		os.system("java -Xmx2g -jar  " + gatk + " -T CombineVariants -R " + RefIndex + " --variant " + outputDir + index_tag + outfile2 + " --variant " + outputDir + index_tag + outfile4 + ' -o ' + outputDir + index_tag +  outfile5 + " -genotypeMergeOptions UNIQUIFY ")# -priority " + outputDir + index_tag + outfile2 + ',' + outputDir + index_tag +  outfile4 ) # + " --assumeIdenticalSamples ")
		#print("java -Xmx2g -jar  " + gatk + " -T -T CombineVariants -R " + RefIndex + " --variant " + outputDir + index_tag + outfile2 + " --variant " + outputDir + index_tag + outfile4 + ' -o ' + outputDir + index_tag +  outfile5 + " --assumeIdenticalSamples ")

		print("final filtering of the merged file")
		os.system("java -Xmx2g -jar  " + gatk + " -T VariantFiltration -R " + RefIndex + ' --variant ' + outputDir + index_tag +  outfile5 + ' -o ' + outputDir + index_tag +  outfile6 + ' --filterExpression "QD < 10.0 || FS > 20.0 || MQ < 50.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || AN > 10 || AF<0.25" --filterName "harsh_filter" ' )
		#print("java -Xmx2g -jar  " + gatk + " -T VariantFiltration -R " + RefIndex + ' --variant ' + outputDir + index_tag +  outfile5 + ' -o ' + outputDir + index_tag +  outfile6 + ' --filterExpression "QD < 10.0 || FS > 20.0 || MQ < 50.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || AF<0.25" --filterName "harsh_filter" ' )


		os.system("rm " + outputDir + index_tag + outfile1)
		os.system("rm " + outputDir + index_tag + outfile1 + '.idx')
		os.system("rm " + outputDir + index_tag + outfile2)
		os.system("rm " + outputDir + index_tag + outfile2 + '.idx')
		os.system("rm " + outputDir + index_tag + outfile3)
		os.system("rm " + outputDir + index_tag + outfile3 + '.idx')
		os.system("rm " + outputDir + index_tag + outfile4)
		os.system("rm " + outputDir + index_tag + outfile4 + '.idx')
		os.system("rm " + outputDir + index_tag + outfile5)
		os.system("rm " + outputDir + index_tag + outfile5 + '.idx')

	if annotation ==True:
	
		outfile7 = "gatk.merged.harsh.filtered.ann.vcf"
		outfile8 = "gatk.merged.harsh.filtered.ann.final.vcf"
		
		print("Use snpEff to annotate a vcf file \n")
		os.system('java -Xmx2g -jar ' + snpEff + ' -c ' + snpEff_config + ' -v R64-1-1.75 -o gatk '+ outputDir + index_tag + outfile6 + ' > '+ outputDir + index_tag + outfile7 + ' -s ' + outputDir + index_tag + '_snpEff_summary.html')
		#print('java -Xmx2g -jar ' + snpEff + ' -c ' + snpEff_config + ' -v R64-1-1.75 -o gatk '+ outputDir + index_tag + outfile6 + ' > '+ outputDir + index_tag + outfile7)
		
		print("concatenate the vcf files using gatk tools \n")
		os.system("java -jar -Xmx2g " + gatk + " -T VariantAnnotator -R " + RefIndex + " -A SnpEff --variant " + outputDir + index_tag + outfile6 + " --snpEffFile " + outputDir + index_tag + outfile7 + " -L " + outputDir + index_tag + outfile6 + " -o " + outputDir + index_tag + outfile8)
		#print("java -jar -Xmx2g " + gatk + " -T VariantAnnotator -R " + RefIndex + " -A SnpEff --variant " + outputDir + index_tag + outfile6 + " --snpEffFile " + outputDir + index_tag + outfile7 + " -L " + outputDir + index_tag + outfile6 + " -o " + outputDir + index_tag + outfile8)

		print("From filtered and annotated VCF file to table")
		os.system("java -jar -Xmx2g " + gatk + " -T VariantsToTable -raw -R " + RefIndex + " -F CHROM -F POS -F REF -F ALT -F FILTER -F SNPEFF_GENE_NAME -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_EFFECT -F QUAL -F DP -GF AD -GF GQ -GF PL -V " + outputDir + index_tag + outfile8 +' -o ' + outputDir + index_tag + outfile8 +'.txt --allowMissingData' )
		#print("java -jar -Xmx2g " + gatk + " -T VariantsToTable -raw -R " + RefIndex + " -F CHROM -F POS -F REF -F ALT -F FILTER -F SNPEFF_GENE_NAME -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_EFFECT -F QUAL -F DP -GF AD -GF GQ -GF PL -V " + outputDir + index_tag + outfile8 +' -o ' + outputDir + index_tag + outfile8 +'.txt --allowMissingData' )
		
		os.system("rm " + outputDir + index_tag + outfile6 + '.idx')
		os.system("rm " + outputDir + index_tag + outfile6)
		os.system("rm " + outputDir + index_tag + outfile7 + '.idx')
		os.system("rm " + outputDir + index_tag + outfile7)

	if deletionCheck == True: # does not work 
		print("Checking coverage at gene deletions")
		#make pile up file
		#samtools mpileup -s AAGAGGCA-CTAAGCCT.dedup.realigned.bam > mpileup.txt
		os.system("samtools mpileup -s " + outputDir + index_tag + '.dedup.realigned.bam > ' + outputDir + index_tag +"_mpileup.txt")

		#call python script
		os.system("python deletion_check.py " + index_tag + " >> deletion_check.txt")

		#remove pile up file
		os.system("rm " + outputDir + index_tag +"_mpileup.txt")

if Filter_SNP == True : # use SNP_parser.py to filter the SNP group, consider the samples are originally the same compare to the Ref.

	#in this SNP_parser file, the same loop creates the same all_index_tags_dic
	# If you modified this loop here you have to modify it in the SNP_parser file.
	os.system('python SNP_parser.py ' + wkdir + ' index_key.txt')
