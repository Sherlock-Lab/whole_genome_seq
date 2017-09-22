require(ggplot2)

args = commandArgs(TRUE)

working_dir = args[1]

setwd(working_dir)

index_tag = args[2]
y_axis_limit = args[3]

dat2=read.csv(paste(index_tag,"_data/",index_tag,"_cov1kb.csv",sep=''), header=TRUE, sep=',')
dat2$chrom=factor(x=dat2$chrom, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","mt"))
p = ggplot(dat2)+geom_point(aes(x=win.start, y=counts, color=chrom))+ylim(0, as.numeric(y_axis_limit))+facet_wrap(~chrom, ncol=4, scale="free")+theme_bw()+theme(legend.position="none")

ggsave(filename=paste(index_tag,"_data/",index_tag,'cov_plot.jpg',sep=''), plot=p, width=11,height=8.5, dpi=500)