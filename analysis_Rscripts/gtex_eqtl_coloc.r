{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(data.table)\
library(coloc)\
\
sentinels<-fread("/scratch16/abattle4/rebecca/final/gwama/sentinels.txt")\
\
context<-fread("/scratch16/abattle4/rebecca/final/gtex/eqtl/sentinel_tissues.txt", header=F)\
colnames(context)<-c("gene", "chrpos", "tissue")\
\
i<-0\
for (snp in sentinels$rs_number)\{\
  i<-i+1\
  print(paste("Working on snp ", i, " out of 56"))\
p=sentinels[rs_number==snp,]$pos\
  \
#import snplist for meta-analysis\
meta<-fread(paste("/scratch16/abattle4/rebecca/final/gwama/2mb_snplist/", snp, sep=""))\
colnames(meta)<-c("chrpos", "chr", "pos", "ref", "other", "eaf", "beta", "se", "beta_95L", "beta_95U", "z", "pvalue", "-log10pvalue", "q_statistic", "q_pvalue", "i2", "n_studies", "n_samples", "effects")\
#this snplist file has a 2Mb window, reduce to 1Mb\
meta<-subset(meta, meta$pos > (p-500000) & meta$pos < (p+500000))\
\
c<-subset(context, context$chrpos==snp)\
\
for (t in unique(c$tissue))\{\
  \
  t_genes<-subset(c, c$tissue==t)\
  fwrite(unique(t_genes[,1]), paste("/scratch16/abattle4/rebecca/final/gtex/eqtl/genelist_", t, sep=""))\
  \
#gather eQTL data for all relevant genes for this tissue\
eqtls<-fread(cmd=paste("grep -w -F -f /scratch16/abattle4/rebecca/final/gtex/eqtl/genelist_", t, " /scratch16/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/", t, ".allpairs.txt", sep=""))\
eqtls[,c("chr", "pos", "a1", "a2", "build"):=tstrsplit(V2, "_", fixed=T)] \
eqtls$chr<-gsub("chr", "", eqtls$chr)\
eqtls$chrpos<-paste(eqtls$chr, eqtls$pos, sep="_")\
eqtls<-subset(eqtls, !is.na(eqtls$V9))\
    #Get rid of any rows that have an NA for the slope_se, this will mess up coloc\
  for (g in unique(t_genes$gene))\{\
    data1<-eqtls[V1==g,]\
\
    #remove chr_pos where multiple alleles are tested\
    dup<-which(duplicated(data1$chrpos))\
    dup<-data1[c(dup),]$chrpos\
    data1<-subset(data1, !(data1$chrpos %in% dup))\
\
    #dataset1 = meta-analysis\
    #dataset2 = gtex eQTL for gene g in tissue t\
    res<-coloc.abf(dataset1=list(beta=meta$beta, varbeta=(meta$se)^2, snp=meta$chrpos, sdY=0.93, type="quant"), dataset2=list(beta=data1$V8, varbeta=data1$V9^2, snp=data1$chrpos, sdY=0.99, type="quant"), p12=1e-6)\
    system(paste("echo", snp, g, t, res[[1]][[2]], res[[1]][[3]], res[[1]][[4]], res[[1]][[5]], res[[1]][[6]], ">> /scratch16/abattle4/rebecca/final/gtex/eqtl/all_results_coloc.txt", sep=" "))\
  \}\
\}\
system(paste("rm /scratch16/abattle4/rebecca/final/gtex/eqtl/genelist_*"))\
\}}