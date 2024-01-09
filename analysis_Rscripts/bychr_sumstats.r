{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(data.table)\
#for each chr read in the meta-analysis sumstats and merge with the ref file\
for (chr in c(1:22, "X"))\{\
print(chr)\
meta<-fread(paste("/scratch16/abattle4/rebecca/final/gwama/no_sb_bychr/no_sb_chr_pos_chr", chr, ".txt", sep=""))\
colnames(meta)<-c("chr_hg38pos", "chr", "pos", "ref", "other", "eaf", "beta", "se", "beta_95L", "beta_95U", "z", "pvalue", "-log10_pvalue", "q_statistic", "q_pvalue", "i2", "n_studies", "n_samples", "effects")\
#make a snpID column, note that "other" is actually the reference allele in the genome browser\
meta$snpID<-paste(meta$chr, meta$pos, meta$other, meta$ref, sep=":")\
#query the snpID list against the reference file we have to fetch matches\
fwrite(meta[,20], paste("/scratch16/abattle4/rebecca/temp/meta_rsid_chr_", chr, sep=""))\
rsid<-fread(cmd=paste("grep -w -F -f /scratch16/abattle4/rebecca/temp/meta_rsid_chr_", chr, " /data/abattle4/rebecca/meta/euro_topmed_hg19.txt", sep=""), header=F)\
#this is capturing a lot of rows that are the extra headers\
rsid<-subset(rsid, !(rsid$V20 == "rsid"))\
#reduce table to just rsid and snpID\
rsid<-rsid[,c(1,20)]\
colnames(rsid)<-c("snpID", "rsid")\
#combine meta-analysis sumstats with rsIDs\
data<-merge(meta, rsid, by="snpID")\
fwrite(data, paste("/scratch16/abattle4/rebecca/final/twas/bychr_sumstats/no_sb_chrpos_", chr, sep=""), col.names=T)\
\}\
\
system(paste("rm /scratch16/abattle4/rebecca/temp/*", sep=""))}