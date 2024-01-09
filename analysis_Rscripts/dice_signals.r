{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(data.table)\
\
sentinels<-fread("/scratch16/abattle4/rebecca/final/dice/sentinels.txt")\
\
celltypes<-c("B_CELL_NAIVE", "CD4_NAIVE", "CD4_STIM", "CD8_NAIVE", "M2", "MONOCYTES", "NK", "TFH", "TH17", "TH1", "TH2", "THSTAR", "TREG_MEM", "TREG_NAIVE")\
\
i<-0\
meta_snp<-vector()\
dice_cell<-vector()\
\
for (cell in celltypes)\{\
  dice<-fread(paste("/data/abattle4/rebecca/ref/dice/dice_hg19/", cell, "_hg19.txt", sep=""))\
  dice$V1<-gsub("chr", "", dice$V1)\
  dice[,c("Gene", "GeneSymbol", "pvalue", "beta", "statistic", "fdr"):=tstrsplit(V8, ";", fixed=T)]\
  dice$Gene<-gsub("Gene=","", dice$Gene)\
  dice$GeneSymbol<-gsub("GeneSymbol=", "", dice$GeneSymbol)\
  dice$pvalue<-gsub("Pvalue=", "", dice$pvalue)\
  dice$beta<-gsub("Beta=", "", dice$beta)\
  dice$statistic<-gsub("Statistic=", "", dice$statistic)\
  dice$fdr<-gsub("FDR=", "", dice$fdr)\
  \
  for (chrpos in sentinels$V1)\{\
    c<-sentinels[V1==chrpos,]$V19\
    p<-sentinels[V1==chrpos,]$V2\
    df <- subset(dice, dice$V1==c & dice$V2 > (p-1000000) & dice$V2 < (p+1000000))\
    if (min(df$pvalue) < 5e-3)\{\
      i<-i+1\
      meta_snp[[i]]<-chrpos\
      dice_cell[[i]]<-cell\
      fwrite(df, paste("/scratch16/abattle4/rebecca/final/dice/dice_snplist/", chrpos, "_", cell, sep=""), sep="\\t")\
    \}\
  \}\
\}\
\
df2<-data.frame(meta_snp, dice_cell)\
fwrite(df2, "/scratch16/abattle4/rebecca/final/dice/dice_signals_index.txt", col.names=T)\
q()\
n}