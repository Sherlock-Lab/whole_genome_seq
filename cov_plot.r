require(ggplot2)
#dat=read.csv("coverage_10kb.csv", header=TRUE, sep=',')
#dat$chrom=factor(x=dat$chrom, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","mt"))
#p = ggplot(dat)+geom_point(aes(x=win.start, y=counts, color=chrom))+ylim(0,200)+facet_wrap(~chrom, ncol=4, scale="free")+theme_bw()+theme(legend.position="none")
#ggsave(filename='test.jpg', plot=p)
setwd("/Volumes/Promise\ Pegasus/Lucas/whole_genome_seq/CLM_clones/")
#p = ggplot(dat)+geom_point(aes(x=win.start, y=counts, color=chrom))+facet_wrap(~chrom, ncol=4, scale="free")+theme_bw()+theme(legend.position="none")
#ggsave(filename='test.jpg', plot=p, width=11,height=8.5, dpi=500)
args = commandArgs(TRUE)
index_tag = args[1]

#file name: /Users/jaffe/Desktop/SherlockLab/whole_genome_seq/working_pipeline/AAGAGGCA-CTAAGCCT_data/AAGAGGCA-CTAAGCCT_cov1kb.csv

dat2=read.csv(paste(index_tag,"_data/",index_tag,"_cov1kb.csv",sep=''), header=TRUE, sep=',')
dat2$chrom=factor(x=dat2$chrom, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","mt"))
p = ggplot(dat2)+geom_point(aes(x=win.start, y=counts, color=chrom))+ylim(0,40)+facet_wrap(~chrom, ncol=4, scale="free")+theme_bw()+theme(legend.position="none")
#p = ggplot(dat2)+geom_point(aes(x=win.start, y=counts, color=chrom))+facet_wrap(~chrom, ncol=4, scale="free")+theme_bw()+theme(legend.position="none")
ggsave(filename=paste(index_tag,"_data/",index_tag,'cov_plot.jpg',sep=''), plot=p, width=11,height=8.5, dpi=500)