{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 args=commandArgs(trailingOnly=T)\
variant=args[1]\
hg38=args[2]\
\
library(data.table)\
library(coloc)\
library(ggplot2)\
library(cowplot)\
\
print(variant)\
hg38<-as.numeric(hg38)\
\
#Import file that has 2Mb window around meta-analysis peak snp\
data<-fread(paste("/scratch16/abattle4/rebecca/final/eqtlgen/snplist/", variant, sep=""), header=F)\
colnames(data)<-c("chr_hg38pos", "hg19pos", "ref", "alt", "eaf", "b", "se", "beta_95L", "beta_95U", "z", "pvalue", "-logpvalue", "q_statistic", "-logq_statistic", "i2", "n_studies", "n_samples", "effects", "chr", "hg38pos")\
data<-subset(data, data$hg38pos > (hg38-1000000) & data$hg38pos < (hg38+1000000))\
#we will use chr_hg19pos as a snp identifier to compare the meta-analysis and eqtlgen\
data$chr<-gsub("chr", "", data$chr)\
data$chr_hg19pos<-paste(data$chr, data$hg19pos, sep="_")\
##Convert EAF to MAF\
noflip<-subset(data, data$eaf < 0.5)\
noflip$maf<-noflip$eaf\
noflip$beta<-noflip$b\
flip<-subset(data, data$eaf > 0.5)\
flip$maf<-1-flip$eaf\
flip$beta<-flip$b*-1\
data<-rbind(noflip, flip)\
\
maf<-data[,c(21, 22)]\
\
#Import eQTLGen data for region and extract gene list\
genelist<-fread(paste("/scratch16/abattle4/rebecca/final/eqtlgen/eqtlgen_snplist/", variant, sep=""), header=F)\
colnames(genelist)<-c("pvalue", "rsid", "chr", "hg19pos", "alt", "ref", "z", "gene", "hgnc", "gene_chr", "gene_pos", "NrCohorts", "NrSamples", "FDR", "BonferroniP")\
genelist<-(unique(genelist$gene))\
\
for (gene in genelist)\{\
  print(gene)\
  eqtl<-fread(cmd=paste("zgrep -w ", gene, " /scratch16/abattle4/lab_data/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"))\
  colnames(eqtl)<-c("pvalue", "rsid", "chr", "hg19pos", "alt", "ref", "z", "gene", "hgnc", "gene_chr", "gene_pos", "NrCohorts", "NrSamples", "FDR", "BonferroniP")\
  eqtl$chr_hg19pos<-paste(eqtl$chr, eqtl$hg19pos, sep="_")\
  eqtl<-merge(eqtl, maf, by="chr_hg19pos")\
  \
  res<-coloc.abf(dataset1=list(pvalues=data$pvalue, MAF=data$maf, N=data$n_samples, snp=data$chr_hg19pos, sdY=0.93, type="quant"), dataset2=list(pvalues=eqtl$pvalue, MAF=eqtl$maf, N=eqtl$NrSamples, snp=eqtl$chr_hg19pos, sdY=0.99, type="quant"), p12=1e-6)\
  system(paste("echo", variant, gene, res[[1]][[2]], res[[1]][[3]], res[[1]][[4]], res[[1]][[5]], res[[1]][[6]], " >> /scratch16/abattle4/rebecca/final/eqtlgen/all_results_coloc.txt", sep=" "))\
\}\
\
system(paste("rm /scratch16/abattle4/rebecca/final/eqtlgen/geno_QTL/", snp, ,"*", sep=""))\
system(paste("rm /scratch16/abattle4/rebecca/final/eqtlgen/snplist_all", variant, sep="_"))\
writeLines(capture.output(sessionInfo()), "/scratch16/abattle4/rebecca/final/eqtlgen/coloc.session")\
\
q()\
n}