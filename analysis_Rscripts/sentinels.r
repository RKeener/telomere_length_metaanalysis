{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 AndaleMono;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(data.table)\
data<-fread("/scratch16/abattle4/rebecca/final/gwama/no_sb.out")\
colnames(data)[10]<-"pvalue"\
data<-data[data$pvalue < 5e-8,]\
\
data<-data[order(data$pvalue),]\
data[,c("chr", "pos"):=tstrsplit(rs_number, "_", fixed=T)]\
data$pos<-as.numeric(data$pos)\
\
for(i in data$rs_number)\{\
 loc=which(data$rs_number==i)\
 c=data$chr[loc]\
 p=as.numeric(data$pos[loc])\
 data1<-subset(data$rs_number, data$chr==c & data$pos>(p-1000000) & data$pos<(p+1000000))\
 data1<-data1[-1]\
 data<-subset(data, !(data$rs_number %in% data1))\
\}\
\
fwrite(data, "/scratch16/abattle4/rebecca/final/gwama/sentinels_unfiltered.txt", col.names=T, row.names=F, sep="/t")\
\
#return to this to adjust column order\
library(data.table)\
data<-fread("/scratch16/abattle4/rebecca/final/gwama/sentinels_unfiltered.txt")\
data2<-data[,c(1,18,19,2,3,10,4:9,11:17)]\
fwrite(data2, "/scratch16/abattle4/rebecca/final/gwama/sentinels_unfiltered.txt", col.names=T, row.names=F, sep=",")\
\
#return again to remove single SNP loci and treat \expnd0\expndtw0\kerning0
6_30796116, 6_29785571, and 6_28301047 as one peak.\kerning1\expnd0\expndtw0 \
sentinels<-fread("/scratch16/abattle4/rebecca/final/gwama/sentinels_unfiltered.txt")\
sentinels<-subset(sentinels, !(sentinels$rs_number %in% c("7_126320597", "3_123867808", "3_122360383", "3_119809416", "3_118785854", "3_117584223", "2_181649232", "6_29785571", "6_28301047")))\
fwrite(sentinels, "/scratch16/abattle4/rebecca/final/gwama/sentinels.txt")}