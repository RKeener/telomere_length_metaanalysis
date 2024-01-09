{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(data.table)\
\
sent<-fread("/scratch16/abattle4/rebecca/final/gwama/sentinels.txt")\
\
known<-fread("/scratch16/abattle4/rebecca/final/novel_loci/published_sentinels_2022.csv")\
\
i<-0\
novel<-vector()\
for (snp in sent$rs_number)\{\
  i<-i+1\
  s<-sent[rs_number==snp,]\
  start<-s$pos-500000\
  end<-s$pos+500000\
  c<-s$chr\
  d<-subset(known, known$Position > start & known$Position < end & known$chr == c)\
  if (nrow(d)==0)\{\
    novel[[i]]<-"novel"\
  \}else\{\
    novel[[i]]<-"previously observed"\
  \}\
\}\
\
sent$novel<-novel\
fwrite(sent, "/scratch16/abattle4/rebecca/final/novel_loci/sentinels_novelty.txt")}