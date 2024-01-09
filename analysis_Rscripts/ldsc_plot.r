{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(data.table)\
library(ggplot2)\
\
#use a supplementary table from Finucane 2018 to get cell types and categories\
fin<-fread("/Users/rebecca_keener/Documents/Telomere-GWAS-project/manuscript/generating_figures/local_analysis/finucane_2018_supp_table7.csv")\
yes<-subset(fin, fin$entex=="Yes")\
yes$Name<-paste(yes$tissue, "_ENTEX__", yes$mark, sep="")\
no<-subset(fin, fin$entex=="No")\
no$Name<-paste(no$tissue, no$mark, sep="__")\
fin<-rbind(yes, no)\
\
#read in LDSC results using chromatin marks\
meta<-fread("/Users/rebecca_keener/Documents/Telomere-GWAS-project/manuscript/generating_figures/local_analysis/gwama_euro_screen_Multi_tissue_chromatin.cell_type_results.txt")\
\
data<-merge(meta, fin, by="Name")\
data$Name<-factor(data$Name, levels=data[order(category),]$Name)\
\
ggplot(data, aes(x=category, y=-log10(Coefficient_P_value), color=category)) +\
  geom_point() +\
  ggtitle("Cell type enrichment for European meta-analysis by chromatin data") +\
  ylab("-log10(pvalue)") +\
  xlab("") +\
  geom_hline(yintercept=2.85, color="grey60", linetype="dashed") +\
  geom_text(data=subset(data, -log10(data$Coefficient_P_value) > 2.85), aes(x=category, y=-log10(Coefficient_P_value), color=category, label=Name)) +\
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none")\
ggsave("ldsc_chromatin.png", path="/Users/rebecca_keener/Documents/Telomere-GWAS-project/manuscript/panels/", type="png")\
\
#read in LDSC results using expression data\
meta<-fread("/Users/rebecca_keener/Documents/Telomere-GWAS-project/manuscript/generating_figures/local_analysis/gwama_euro_Multi_tissue_gene_expr.cell_type_results.txt")\
\
#use Finucane 2018 Supplemental Table 6 to get categories\
fin<-fread("/Users/rebecca_keener/Documents/Telomere-GWAS-project/manuscript/generating_figures/local_analysis/finucane_2018_supptab6_categories.csv")\
colnames(fin)[1]<-"Name"\
\
data<-merge(meta, fin, by="Name")\
data$Name<-factor(data$Name, levels=data[order(category),]$Name)\
\
ggplot(data, aes(x=category, y=-log10(Coefficient_P_value), color=category)) +\
  geom_point() +\
  ggtitle("Cell type enrichment for European meta-analysis by expression data") +\
  ylab("-log10(pvalue)") +\
  xlab("") +\
  geom_hline(yintercept=2.85, color="grey60", linetype="dashed") +\
  geom_text(data=subset(data, -log10(data$Coefficient_P_value) > 2.85), aes(x=category, y=-log10(Coefficient_P_value), color=category, label=Name)) +\
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none")\
ggsave("ldsc_expression.png", path="/Users/rebecca_keener/Documents/Telomere-GWAS-project/manuscript/panels/", type="png")}