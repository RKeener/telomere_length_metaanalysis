{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 args=commandArgs(trailingOnly=T)\
chr=args[1]\
position=args[2]\
variant=args[3]\
\
library(data.table)\
library(susieR)\
library(ggplot2)\
library(Rfast)\
\
print(variant)\
\
#read in TOPMed pooled GWAS summary statistics \
data<-fread(paste("/scratch16/abattle4/rebecca/full_042020/chr", chr, "_telomere_adjagesexstudyseqctrbatchPCs_minDP0.csv.gz", sep=""))\
\
#subset to locus\
coord<-fread("/scratch16/abattle4/rebecca/final/susie/susie_coordinates.csv")\
coord<-subset(coord, coord$chrpos==variant)\
coord$start<-as.numeric(coord$start)\
coord$end<-as.numeric(coord$end)\
data<-subset(data, data$pos >= coord$start & data$pos <= coord$end)\
data<-unique(data)\
\
#complete formatting of data\
qc<-fread("/scratch16/abattle4/rebecca/snpIDsToDrop_combined_unique_20200515.txt", header=F)\
data<-subset(data, !(data$snpID %in% qc$V1))\
data$beta<-data$Score/(data$Score.SE^2)\
data$se<-1/data$Score.SE\
data$z<-data$beta/data$se\
\
fwrite(data[,19], paste("/scratch16/abattle4/rebecca/final/susie/snplist_", variant, sep=""))\
\
#calculate R matrix\
#use 15K individuals from TOPMed for LD reference\
system(paste("/data/abattle4/rebecca/software/plink2 --pfile /scratch16/abattle4/rebecca/ld_082022/freeze.8.chr", chr, ".pass_and_fail.gtonly.minDP0_noSARP --extract /scratch16/abattle4/rebecca/final/susie/snplist_", variant, " --keep /scratch16/abattle4/rebecca/final/susie_15K/topmed_15K.ind --make-bed --silent --out /scratch16/abattle4/rebecca/final/susie/ld/", variant, sep=""))\
\
system(paste("/data/abattle4/rebecca/software/plink --bfile /scratch16/abattle4/rebecca/final/susie/ld/", variant, " --recode A --keep-allele-order --write-snplist --silent --out /scratch16/abattle4/rebecca/final/susie/ld/", variant, sep=""))\
\
r<-fread(paste("/scratch16/abattle4/rebecca/final/susie/ld/", variant, ".raw", sep=""))\
r<-r[,-c(1:6)]\
allzero<-which(r[,colSums(r)!=0])\
r<-r[,..allzero]\
allone<-which(r[,colSums(r)]!=nrow(r))\
r<-r[,..allone]\
alltwo<-which(r[,colSums(r)]!=2*nrow(r))\
r<-r[,..alltwo]\
\
check<-data.table(id=colnames(r),placeholder=rep(NA,length(colnames(r))))\
check[,c("snpID", "variant"):=tstrsplit(id, "_", fixed=T)]\
snplist<-check$snpID\
z<-subset(data$z, data$snpID %in% snplist)\
\
R<-data.matrix(r)\
rm(r)\
R<-cor(R, use="complete.obs")\
\
#n_samples is fixed in TOPMed data\
n=109122\
\
fitted<-susie_rss(z, R, L=10, n=n)\
top<-data.frame(snpID=snplist, pip=fitted$pip, zscore=z)\
fwrite(top, paste("/scratch16/abattle4/rebecca/final/susie/allsnps/", variant, "_new_allsnps.txt", sep=""), col.names=T, row.names=F, sep="\\t", quote = F)\
\
num<-length(fitted$sets$cs)\
print(paste(variant, " is predicted to have ", num, " credible sets.", sep=""))\
\
for (set in c(1:num))\{\
  i<-fitted$sets$cs[[set]]\
  for (n in c(1:length(i)))\{\
    j<-i[n]\
      if (set == 1 && n == 1)\{\
        df = data.frame(j, snplist[j], z[j], fitted$pip[[j]])\
        colnames(df)<-c("number", "snpID", "zscore", "pip")\
        df$set_num<-rep(set, nrow(df))\
        data2<-df\
      \}else\{\
        df = data.frame(j, snplist[j], z[j], fitted$pip[[j]])\
        colnames(df)<-c("number", "snpID", "zscore", "pip")\
        df$set_num<-rep(set, nrow(df))\
        data2<-rbind(df, data2)\
      \}\
  \}\
\}\
\
data2<-data.frame(data2)\
colnames(data2)<-c("number", "snpID", "zscore", "pip", "cs_num")\
fwrite(data2, paste("/scratch16/abattle4/rebecca/final/susie/credset/", variant, "_new_credset.txt", sep=""), col.names=T, row.names=F, sep="\\t", quote = F)\
\
#plot credset\
full_snp<-data[Score.pval==min(Score.pval),]$snpID\
\
system(paste("/data/abattle4/rebecca/software/plink --bfile /scratch16/abattle4/rebecca/final/susie/ld/", variant, " --r2 inter-chr --ld-snp ", full_snp, " --ld-window-r2 0 --silent --out /scratch16/abattle4/rebecca/final/susie/ld/", variant, sep=""))\
    #Change to Plink1 to calculate LD between each SNP with the lead SNP. Don't specify a window because we don't want to limit it.\
\
ld<-fread(paste("/scratch16/abattle4/rebecca/final/susie/ld/", variant, ".ld", sep=""))\
ld<-ld[,c(6,7)]\
colnames(ld)<-c("snpID", "r2")\
\
data<-merge(data, ld, by="snpID")\
\
cs<-subset(data, data$snpID %in% data2$snpID)\
credset<-merge(data2, cs, by="snpID")\
credset$cs_num<-as.character(credset$cs_num)\
\
png(filename=paste("/scratch16/abattle4/rebecca/final/susie/manplot/", variant, "_new_credset.png", sep=""))\
ggplot(data, aes(pos, -log10(Score.pval), color=r2)) +\
  geom_point() +\
  geom_point(data=credset, aes(pos, -log10(Score.pval)), inherit.aes = F, size=3, shape=23) +\
  scale_x_continuous(name="Position on chromosome") +\
  ylab("-log10(P-value)") +\
  ggtitle(paste("New SUSIE TOPMed GWAS + all TOPMed LD ref ", variant, sep="")) +\
  scale_color_gradientn(colours = rainbow(5)) \
dev.off()\
\
\
#to generate 99% credible set for KBTBD6/KBTBD7 locus\
if (variant==\'9313_41150640\'94)\{\
fitted<-susie_rss(z, R, L=10, n=n, coverage=0.99)\
fwrite(top, paste("/scratch16/abattle4/rebecca/final/susie/allsnps/", variant, "_new_allsnps_0.99credset.txt", sep=""), col.names=T, row.names=F, sep="\\t", quote = F)\
\
num<-length(fitted$sets$cs)\
\
for (set in c(1:num))\{\
  i<-fitted$sets$cs[[set]]\
  for (n in c(1:length(i)))\{\
    j<-i[n]\
      if (set == 1 && n == 1)\{\
        df = data.frame(j, snplist[j], z[j], fitted$pip[[j]])\
        colnames(df)<-c("number", "snpID", "zscore", "pip")\
        df$set_num<-rep(set, nrow(df))\
        data2<-df\
      \}else\{\
        df = data.frame(j, snplist[j], z[j], fitted$pip[[j]])\
        colnames(df)<-c("number", "snpID", "zscore", "pip")\
        df$set_num<-rep(set, nrow(df))\
        data2<-rbind(df, data2)\
      \}\
  \}\
\}\
\
data2<-data.frame(data2)\
colnames(data2)<-c("number", "snpID", "zscore", "pip", "cs_num")\
fwrite(data2, paste("/scratch16/abattle4/rebecca/final/susie/credset/", variant, "_new_credset_0.99credset.txt", sep=""), col.names=T, row.names=F, sep="\\t", quote = F)\
\}\
\
writeLines(capture.output(sessionInfo()), "/scratch16/abattle4/rebecca/final/susie/run_susie.session")\
\
q()\
n}