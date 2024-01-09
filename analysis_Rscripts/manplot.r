{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 library(data.table)\
library(ggplot2)\
library(dplyr)\
\
data<-fread("/scratch16/abattle4/rebecca/final/gwama/no_sb.out")\
\
colnames(data)<-c("rs_number", "ref_allele", "other_allele", "eaf", "beta", "se", "beta_95L", "beta_95U", "z", "pvalue", "-log10pvalue", "q_statistic", "q_pvalue", "i2", "n_studies", "n_samples", "effects")\
data[,c("chr", "pos"):=tstrsplit(rs_number, "_", fixed=T)]\
\
#dplyr code requires dataset to be called data\
data<-subset(data, data$chr %in% c(1:22, "X"))\
data$chr<-sub("X", 23, data$chr)\
data$chr<-as.numeric(data$chr)\
data$pos<-as.numeric(data$pos)\
#dplyr code requires chromosome column be called CHR\
colnames(data)[18]<-"CHR"\
#dplyr code requires position column be called BP\
colnames(data)[19]<-"BP"\
#dplyr code requires p-value column be called P\
colnames(data)[10]<-"P"\
\
t <- data %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(chr_len)-chr_len) %>% dplyr::select(-chr_len) %>% left_join(data, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot)\
\
axisdf = t %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )\
\
ggplot(t, aes(x=BPcum, y=-log10(P))) + \
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +\
  scale_color_manual(values = rep(c("grey10", "grey60"), 22 )) +\
  scale_x_continuous(name="Position on Chromosome", label = axisdf$CHR, breaks= axisdf$center ) +\
  scale_y_continuous(name="-log10(P-value)") +\
  theme(legend.position="none", panel.grid.major.x = element_blank()) +\
  geom_hline(yintercept=5, color="blue") +\
  geom_hline(yintercept=9, color="red") +\
  ggtitle("TL Meta-analysis Manhattan Plot")\
ggsave("manplot_nolabels.png", device="png", path="/scratch16/abattle4/rebecca/final/manplot/", width=13, height=7, units="in")\
\
#add novel loci in color\
n<-fread("/scratch16/abattle4/rebecca/final/novel_loci/sentinels_novelty.txt")\
n<-subset(n, n$novel=="novel")\
\
i<-0\
for (snp in n$rs_number)\{\
  i<-i+1\
  s<-n[rs_number==snp,]\
  d<-subset(t, t$CHR==s$chr & t$BP > (s$pos-150000) & t$BP < (s$pos+150000))\
  if (i==1)\{\
    novel<-d\
  \}else\{\
    novel<-rbind(novel,d)\
  \}\
\}\
\
ggplot(t, aes(x=BPcum, y=-log10(P))) + \
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +\
  scale_color_manual(values = rep(c("grey10", "grey60"), 22 )) +\
  scale_x_continuous(name="Position on Chromosome", label = axisdf$CHR, breaks= axisdf$center ) +\
  scale_y_continuous(name="-log10(P-value)") +\
  theme(legend.position="none", panel.grid.major.x = element_blank(), panel.background = element_rect(fill="white")) +\
  geom_hline(yintercept=5, color="blue") +\
  geom_hline(yintercept=9, color="red") +\
  ggtitle("TL Meta-analysis Manhattan Plot") +\
  geom_point(data=novel, color="blue")\
ggsave("manplot_novel_blue.png", device="png", path="/scratch16/abattle4/rebecca/final/manplot/", width=13, height=7, units="in")}