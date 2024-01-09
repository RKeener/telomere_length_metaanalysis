{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 args=commandArgs(trailingOnly=T)\
variant=args[1]\
\
library(data.table)\
library(ggplot2)\
library(coloc)\
library(cowplot)\
\
\
#import snplist for meta-analysis\
data<-fread(paste("/scratch16/abattle4/rebecca/final/dice/snplist/" ,variant, sep=""), header=F)\
colnames(data)<-c("chr_hg38pos", "hg19pos", "ref_allele", "other_allele", "eaf", "beta", "se", "beta_95L", "beta_95U", "z", "pvalue", "-logpvalue", "q_statistic", "q_pvalue", "i2", "n_studies", "n_samples", "effects", "chr", "hg38pos")\
c<-unique(data$chr)\
c<-gsub("chr", "", c)\
data$chr_hg19pos<-paste(c, data$hg19pos, sep="_")\
\
#identify cell types of interest\
cell_list<-fread(cmd=paste("grep ", variant,  " /scratch16/abattle4/rebecca/final/dice/dice_signals_index.txt", sep=""), header=F)\
\
if (nrow(cell_list) == 0)\{\
  system(paste("echo", variant, ">> /scratch16/abattle4/rebecca/final/dice/variant_not_tested.txt", sep=" "))\
  print("No DICE cell types with signals")\
  PPH4=0\
\}else\{\
  #run coloc by cell type by gene\
  for(cell in cell_list$V2)\{\
\
  dice<-fread(paste("/scratch16/abattle4/rebecca/final/dice/dice_snplist/", variant, "_", cell, sep=""))\
  genes<-unique(dice$Gene)    \
\
    for(gene in genes)\{\
    print(paste(variant, gene, cell, sep=" "))\
    \
    eqtls<-fread(cmd=paste("grep -w ", gene, " /data/abattle4/rebecca/ref/dice/dice_hg19/", cell, "_hg19.txt", sep=""))\
    eqtls[,c("Gene", "GeneSymbol", "pvalue", "beta", "statistic", "fdr"):=tstrsplit(V8, ";", fixed=T)]\
    eqtls$Gene<-gsub("Gene=","", eqtls$Gene)\
    eqtls$GeneSymbol<-gsub("GeneSymbol=", "", eqtls$GeneSymbol)\
    eqtls$pvalue<-gsub("Pvalue=", "", eqtls$pvalue)\
    eqtls$beta<-gsub("Beta=", "", eqtls$beta)\
    eqtls$statistic<-gsub("Statistic=", "", eqtls$statistic)\
    eqtls$fdr<-gsub("FDR=", "", eqtls$fdr)\
    eqtls$pvalue<-as.numeric(eqtls$pvalue)\
    eqtls$beta<-as.numeric(eqtls$beta)\
    eqtls$statistic<-as.numeric(eqtls$statistic)\
    eqtls$se<-eqtls$beta/eqtls$statistic\
    eqtls$chr_hg19pos<-paste(c, eqtls$V2, sep="_")\
    #I've seen that for some SNPs the SE = 0, these will mess up coloc, must be removed.\
    eqtls<-subset(eqtls, !(eqtls$se==0))\
    \
    #do not evaluate genes with eQTLs in this cell type that have a max -log10(pvalue)>2\
    if (max(-log10(eqtls$pvalue)) < 2)\{\
      PPH4=0\
      print("max -log10(pvalue) < 2")\
    \}else\{\
      print("Proceed with coloc")\
      #I noticed that in some cases where coloc.abf gives me all NAs, the beta values for the DICE dataset are quite large, these ggplot() lines are to compare the beta distributions.\
      #not needed to run the analysis\
      #ggplot(eqtls, aes(beta)) + geom_density() + ggtitle(paste("Beta distribution in DICE for ", unique(eqtls$GeneSymbol), " near sentinel ", variant, sep=""))\
      #ggplot(data, aes(beta)) + geom_density() + ggtitle(paste("Beta distribution in meta-analysis near sentinel ", variant, sep=""))\
      \
      #in some cases there aren't any shared SNPs between the meta-analysis peak and eQTLGen\
      #these cases cannot be evaluated. Obvi\
      shared<-subset(data, data$chr_hg19pos %in% eqtls$chr_hg19pos)\
      \
      if (nrow(shared)>0)\{\
        print(paste(nrow(shared), " SNPs shared between meta-analysis and DICE", sep=""))\
        #in some cases DICE has multiple alleles tested at the same position. We need to resolve this\
        #remove duplicates in DICE\
        dup<-which(duplicated(eqtls$chr_hg19pos))\
        if (length(dup)>1)\{\
          eqtls<-eqtls[-dup,]\
        \}\
        \
        #run coloc\
        res<-coloc.abf(dataset1=list(beta=data$beta, varbeta=(data$se)^2, snp=data$chr_hg19pos, sdY=0.99, type="quant"), dataset2=list(beta=eqtls$beta, varbeta=(eqtls$se)^2, snp=eqtls$chr_hg19pos, sdY=0.99, type="quant"), p12=1e-6)\
        system(paste("echo", variant, gene, cell, nrow(shared), res[[1]][[2]], res[[1]][[3]], res[[1]][[4]], res[[1]][[5]], res[[1]][[6]], ">> all_results_coloc.txt", sep=" "))\
    \
        PPH4<-res[[1]][[6]]\
        \
        if (is.na(PPH4))\{\
          print("PPH4 was NA")\
          PPH4=0\
        \}\
      \}else\{\
        print("No shared SNPs between datasets")\
        PPH4=0\
      \}\
      print(paste("PPH4=", PPH4, sep=""))\
    \}\
    \
    #make Manhattan plots of the meta-analysis peak and DICE peak in cases where coloc PPH4 > 0.3\
    #this is extrememly conservative to ensure we get anything we might be interested in\
    #a strong result would have something closer to PPH4 ~ 0.7\
    if (PPH4 > 0.3)\{\
      print("Proceed with plots")\
      #need to fwrite a list of SNPs in TOPMed format (chr:pos:ref:alt)\
 data$snpID<-paste(c, data$hg38pos, data$other_allele, data$ref_allele, sep=":")\
      fwrite(data[,22], paste("/scratch16/abattle4/rebecca/final/dice/snplist_all", variant, sep="_"), quote=F, row.names=F, sep="\\t")\
      \
      hg38<-gsub(paste(c, "_", sep=""), "", variant)\
      full_snp<-data[data$hg38pos==hg38,]$snpID\
  \
system(paste("/data/abattle4/rebecca/software/plink2 --pfile /scratch16/abattle4/rebecca/ld_082022/freeze.8.chr", c, ".pass_and_fail.gtonly.minDP0_noSARP --extract /scratch16/abattle4/rebecca/final/dice/snplist_all_", variant, " --exclude /scratch16/abattle4/rebecca/snpIDsToDrop_combined_unique_20200515.txt --keep ~/scratch16-abattle4/rebecca/final/final_meta.inds --make-bed --silent --out /scratch16/abattle4/rebecca/final/dice/geno_QTL/", variant, "_", gene, sep=""))\
    #pass plink2 files to bed format\
\
system(paste("/data/abattle4/rebecca/software/plink --bfile /scratch16/abattle4/rebecca/final/dice/geno_QTL/", variant, "_", gene, " --r2 inter-chr --ld-snp ", full_snp, " --ld-window-r2 0 --silent --out /scratch16/abattle4/rebecca/final/dice/geno_QTL/", variant,"_", gene, sep=""))\
    #use plink1 to calculate LD\
\
ld<-fread(paste("/scratch16/abattle4/rebecca/final/dice/geno_QTL/", variant,"_",gene,".ld", sep=""))\
ld<-ld[,c(4:7)]\
if (c == "X")\{\
ld$CHR<-rep("X", nrow(ld))\
ld$chr_hg38pos<-paste(ld$CHR, ld$pos, sep="_")\
ld<-ld[,c(4,6)]\
\}else\{\
ld$chr_hg38pos<-paste(ld$CHR_B, ld$BP_B, sep="_")\
ld<-ld[,c(4:5)]\
\}     \
\
#Manhattan plot for meta-analysis signal\
data<-merge(data, ld, by="chr_hg38pos")\
excluded<-subset(data, !(data$chr_hg19pos %in% eqtls$chr_hg19pos))\
g1<-subset(data, data$chr_hg38pos==variant)\
\
p1<- ggplot(data, aes(x=as.numeric(hg19pos), y=-log10(as.numeric(pvalue)), color=R2)) +\
      geom_point() +\
      scale_y_continuous(name="-log10(p-value)") + \
      scale_x_continuous(name="Position on chromosome") +\
      geom_point(data=g1, aes(x=as.numeric(hg19pos), y=-log10(as.numeric(pvalue))), inherit.aes=F, fill="black",shape=23, size=3) +\
      geom_point(data=excluded, aes(x=as.numeric(hg19pos), y=-log10(as.numeric(pvalue))), inherit.aes=F, color="grey60", alpha=1/10) +\
      ggtitle(paste("Meta-analysis signal near ", variant, sep="")) +\
      scale_color_gradientn(colours = rainbow(5))\
\
#Manhattan plot for DICE eQTL, match x-axis\
#from ld keep chr_hg19pos, r2, chr_hg38pos\
ld<-data[,c(21:23)]\
eqtls<-merge(eqtls, ld, by="chr_hg19pos")\
g1<-subset(eqtls, eqtls$chr_hg38pos==variant)\
xmin<-min(data$hg19pos)\
xmax<-max(data$hg19pos)\
eqtls<-subset(eqtls, eqtls$V2 > xmin & eqtls$V2 < xmax)\
hgnc<-unique(eqtls$GeneSymbol)\
\
p2<- ggplot(eqtls, aes(x=as.numeric(V2), y=-log10(as.numeric(pvalue)), color=R2)) +\
  geom_point() +\
  scale_y_continuous(name="-log10(p-value)") + \
  scale_x_continuous(name="Position on chromosome") +\
  geom_point(data=g1, aes(x=as.numeric(V2), y=-log10(as.numeric(pvalue))), inherit.aes=F, fill="black",shape=23, size=3) + \
  ggtitle(paste(variant, cell, hgnc, sep=" ")) +\
  scale_color_gradientn(colours = rainbow(5))\
\
plot_grid(p1, p2, ncol=1)\
ggsave(paste(paste(variant, cell, gene, sep="_"), ".png", sep=""), device="png", width=6, height=4, units=c("in"), path="/scratch16/abattle4/rebecca/final/dice/man_plots/")\
\
\}\
\}\
\}\
\}\
  \
  system(paste("rm /scratch16/abattle4/rebecca/final/dice/snplist_all", variant, sep="_"))\
  system(paste("rm /scratch16/abattle4/rebecca/final/dice/geno_QTL/", variant, "*", sep=""))\
\
q()\
n}