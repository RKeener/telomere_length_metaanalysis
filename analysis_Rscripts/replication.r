{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(data.table)\
\
#identify novel signals\
n<-fread("/Users/rebecca_keener/Documents/Telomere-GWAS-project/meta-analysis/final/novel_loci/sentinels_novelty.txt")\
n<-subset(n, novel=="novel")\
\
#which of these novel signal sentinels are present in the UKBB dataset?\
#I need the rsIDs\
rsid<-fread("/Users/rebecca_keener/Documents/Telomere-GWAS-project/meta-analysis/final/sentinels.csv")\
n<-n[,c(1,20)]\
rsid<-merge(rsid, n, by="rs_number")\
\
#look for these rsIDs in the UKBB sumstats\
fwrite(rsid[,2], "/Users/rebecca_keener/Documents/Telomere-GWAS-project/meta-analysis/final/replication_analysis/rsid_list.txt", col.names=F)\
data<-fread(cmd=paste("zgrep -w -F -f /Users/rebecca_keener/Documents/Telomere-GWAS-project/meta-analysis/final/replication_analysis/rsid_list.txt /Users/rebecca_keener/Desktop/datasets/metaanalysis/UKB_telomere_gwas_summarystats.tsv.gz", sep=""))\
colnames(data)<-c("rsid", "p_value", "chr", "bp", "effec_allele", "other_allele", "effect_allele_frequency", "beta", "se")\
\
#significance threshold for replication \
alpha=0.05/nrow(data)\
\
rep<-subset(data, data$p_value < alpha)\
\
#pull the beta and beta.se from the meta-analysis and all input studies, where available, for comparison\
snplist<-data.table(V1=c("10_97357340", "8_257819"))\
fwrite(snplist[,1], "/scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt")\
#read in data from meta-analysis\
meta<-fread(cmd=paste("grep -w -F -f /scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt /scratch16/abattle4/rebecca/final/gwama/no_sb.out", sep=""))\
meta$dataset<-"meta-analysis"\
colnames(meta)<-c("chrpos", "ref", "other", "eaf", "beta", "se", "beta_95L", "beta_95U", "z", "pvalue", "log10pvalue", "q_statistic", "q_pvalue", "i2", "n_studies", "n_samples", "effects", "dataset")\
meta<-meta[,c(1,5,6,10,18)]\
\
#read in data from Dorajoo et al.\
dor<-fread(cmd=paste("grep -w -F -f /scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt /scratch16/abattle4/rebecca/allsnps/dorajoo_allsnps.txt", sep=""))\
dor$dataset<-"Dorajoo et al"\
colnames(dor)<-c("rsid", "chr", "pos", "nea", "ea", "pvalue", "beta", "se", "p_het", "eaf", "chrpos", "n", "dataset")\
dor<-dor[,c(6,7,8,11,13)]\
\
#read in data from Li et al.\
li<-fread(cmd=paste("grep -w -F -f /scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt /scratch16/abattle4/rebecca/allsnps/li_allsnps.txt", sep=""))\
li$dataset<-"Li et al"\
colnames(li)<-c("rsid", "chr", "pos", "hg37pos", "ea", "nea", "eaf", "beta", "se", "pvalue", "n", "FDR", "chrpos", "dataset")\
li<-li[,c(8,9,10,13,14)]\
\
#read in data from Delgado et al.\
del<-fread(cmd=paste("zgrep -w -F -f /scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt /scratch16/abattle4/rebecca/allsnps/delgado_allsnps.txt.gz", sep=""))\
del$dataset<-"Delgado et al"\
colnames(del)<-c("snp", "chr", "pos", "chrhg19", "hg19pos", "EA", "NEA", "EAF", "beta", "se", "pvalue", "chrpos", "n", "dataset")\
del<-del[,c(9:12,14)]\
\
#read in TOPMed European\
euro<-fread(cmd=paste("zgrep -w -F -f /scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt /scratch16/abattle4/rebecca/allsnps/european_standardized_allsnps.txt", sep=""))\
euro$dataset<-"Taub et al European"\
euro<-euro[,c(11,23:26)]\
colnames(euro)<-c("pvalue", "beta", "se", "chrpos", "dataset")\
\
#read in TOPMed African\
african<-fread(cmd=paste("zgrep -w -F -f /scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt /scratch16/abattle4/rebecca/allsnps/african_standardized_allsnps.txt", sep=""))\
african$dataset<-"Taub et al African"\
african<-african[,c(11,23:26)]\
colnames(african)<-c("pvalue", "beta", "se", "chrpos", "dataset")\
\
#read in TOPMed Asian\
asian<-fread(cmd=paste("zgrep -w -F -f /scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt /scratch16/abattle4/rebecca/allsnps/asian_standardized_allsnps.txt", sep=""))\
asian$dataset<-"Taub et al Asian"\
asian<-asian[,c(11,23:26)]\
colnames(asian)<-c("pvalue", "beta", "se", "chrpos", "dataset")\
\
#read in TOPMed Hispanic/Latino\
his<-fread(cmd=paste("zgrep -w -F -f /scratch16/abattle4/rebecca/final/manuscript_figures/tmp/rep_snplist.txt /scratch16/abattle4/rebecca/allsnps/hispanic_standardized_allsnps.txt", sep=""))\
his$dataset<-"Taub et al Hispanic/Latino"\
his<-his[,c(11,23:26)]\
colnames(his)<-c("pvalue", "beta", "se", "chrpos", "dataset")\
\
#read in UKBB data\
ukbb1<-fread(cmd=paste("zgrep -w rs3008267 /data/abattle4/rebecca/meta/original_data/UKB_telomere_gwas_summarystats.tsv.gz", sep=""))\
ukbb1$chrpos<-"8_257819"\
ukbb2<-fread(cmd=paste("zgrep -w rs7923385 /data/abattle4/rebecca/meta/original_data/UKB_telomere_gwas_summarystats.tsv.gz", sep=""))\
ukbb2$chrpos<-"10_97357340"\
ukbb<-rbind(ukbb1, ukbb2)\
colnames(ukbb)<-c("rsid", "pvalue", "chr", "pos", "EA", "NEA", "EAF", "beta", "se", "chrpos")\
ukbb<-ukbb[,c(2,8:10)]\
ukbb$dataset<-"Codd et al"\
\
#combine results into one data frame\
data<-do.call("rbind", list(meta, dor, li, del, euro, african, asian, his, ukbb))\
\
fwrite(data, "/scratch16/abattle4/rebecca/final/replication_analysis/replicated_snps.csv", sep=",", col.names = T)}