#this should all be run from the directory containing the seq_pipeline.py program
#this directory should also contain a folder of fastq files from the sequencer

import os

#define which steps need to be carried out
mkOutputDir = False
trim = False
mapping = True
coverage = False
CNV = False # does not work 
variantCall = True
annotation = True
deletionCheck = False # does not work 

#define working directory from which this program is being run
wkdir = '/Volumes/Promise\ Pegasus/Lucas/whole_genome_seq/FLC_clones/Replicat_2F'
#wkdir = '~/Desktop/SherlockLab/revision_work/sequence_analysis/'

#define directory where all the fastq data files are and file prefix and suffix
dataDir = wkdir+'fastq_files/'
#fastqPrefix = '141126_PINKERTON_0343_BC4J1PACXX_L7_'
r1FileSuffix = 'L004_R1_001.fastq'
r2FileSuffix = 'L004_R2_001.fastq'

#all_index_tags_dic  = {'P1-G12':'GACGTCGA-GAGCCTTA_S84', 'P1-C6':'CTGCGCAT-CCTAGAGT_S30'}
all_index_tags_dic  = {}
index_file = open ('index_key_NextSeq_KGS_F.txt', 'r')
sample_count = -1
for line in index_file:
	(strain, I7, seq7, I5, seq5) = line.rstrip('\n').split('\t') # split the line
	sample_count += 1
	if strain.startswith('P'):
		all_index_tags_dic[strain] = seq7 + '-' + seq5 + '_S' + str(sample_count)
index_file.close()

#define where cutadapt program is so it can be called and its parameters
#cutAdapt = '/Users/miajaffe/whole_genome/software/cutadapt-1.7.1/bin/cutadapt'
#cutAdapt = '/Users/jaffe/seq_software/cutadapt-1.7.1/bin/cutadapt'
beginningTrim = str(10)
qualityCutoff = str(30)
lengthThresh = str(12)
read1Adaptor = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
read2Adaptor = 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'

#define bwa and parameters
RefIndex = '/Volumes/Promise\ Pegasus/Lucas/whole_genome_seq/Reference_index/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa'
#The reference is already build !!!!!!
# Best to be in the Reference_index directory and build those there.
# bwa index fasta_file.fa
# samtools faidx fasta_file.fa
# java -jar /Volumes/Promise\ Pegasus/picard.jar CreateSequenceDictionary R=fasta_file.fa O=path/to/fasta_file.dict

#define picard directory
picard = '/Volumes/Promise\ Pegasus/picard.jar'

#define GATK directory
gatk = '/Volumes/Promise\ Pegasus/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
#RefIndex = '/Volumes/Promise\ Pegasus/Lucas/whole_genome_seq/Reference_gatk/S288C_reference_sequence_R64-2-1_20150113.fsa'

#define snpEff directory
snpEff = '/Volumes/Promise\ Pegasus/snpEff/snpEff.jar'
snpEff_config = '/Volumes/Promise\ Pegasus/snpEff/snpEff.config'

for strain in all_index_tags_dic:
	index_tag = strain + '-' + all_index_tags_dic[strain]
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
		#python coverage.py 10000 coverage_10kb.csv
		os.system("python coverageLH.py " + index_tag + " 1000 ")

		#remove pile up file
		print("removing pileup")
		os.system("rm " + outputDir + index_tag +"_mpileup.txt")

		#make plots
		print("plotting")
		#Rscript cov_plot.r 
		os.system("Rscript cov_plot.r " + index_tag)
		print("Rscript cov_plot.r " + index_tag)


	if variantCall == True:
		print("Calling variants with GATK's Unified Genotyper")
		os.system("java -jar -Xmx2g " + gatk + " -T UnifiedGenotyper -R " + RefIndex + " -I " + outputDir + index_tag + '.dedup.realigned.bam -ploidy 2 --genotype_likelihoods_model BOTH -stand_call_conf 30 -stand_emit_conf 10 -o ' + outputDir + index_tag + '_gatk.vcf')
		#print("java -jar -Xmx2g " + gatk + " -T UnifiedGenotyper -R " + RefIndex + " -I " + outputDir + index_tag + '.dedup.realigned.bam -ploidy 2 --genotype_likelihoods_model BOTH -stand_call_conf 30 -stand_emit_conf 10 -o ' + outputDir + index_tag + '_gatk.vcf')
		
		print("Filtering using Sandeep's threshold")
		
		
		#GATK recommended filtering
		#brute force filtering
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
		
