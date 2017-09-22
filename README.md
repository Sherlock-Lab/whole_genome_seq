# whole_genome_seq

This script is made to analyze whole genome sequencing data using a simple commande line.

You need to set up your working director, your file containing the index of the different sample and your fastq files names in the seq_pipleine file.
Make sure you match the fastq files names with the suffix because this is variable.

This script will:

  1 Trim using cutadapt 
  2 map using bwa to produce in the end the realigned .bam file
  3 make coverage plot using python and R script (make sure to adjust the yaxis limit according to your expected coverage)
  4 Call variant following GATK recommandation
  5 Annontate Variant using SNPeff (recommanded by GATK)to produce in the end a table containing all the mutations and informations.
  6 Once the previous setps done for all samples the script will merged all the final table file and parse the mutations to filter according to frequency of the mutations in all the samples, coverage, mapping quality. The threshold can be modified in the SNP_parser.py file. This produce files that you will still have to manually curate to figure out the real mutations.
  
  
