{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 AndaleMono;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red255\green255\blue255;\red0\green0\blue0;
}
{\*\expandedcolortbl;;\cssrgb\c0\c1\c1;\cssrgb\c100000\c100000\c100000;\cssrgb\c0\c0\c0;
}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs24 \cf2 \cb3 \CocoaLigature0 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
A. Meta-analysis and lead SNP analysis\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\CocoaLigature1 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf2 \
#run meta-analysis\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf2 \CocoaLigature0 ~/work/marios/GWAMA/GWAMA -i /work-zfs/abattle4/rebecca/meta/GWAMA/strand/gwama.in -o /work-zfs/abattle4/rebecca/meta/GWAMA/strand/no_sam_braz/no_sb -qt --random\
\
#identify sentinel SNPs\
Rscript sentinels.r\
\
#identify novel loci\
Rscript novel.r\
\
#make Manhattan plot\
Rscript manplot.r\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
B. Replication analysis\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Rscript replication.r\
\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C. Colocalization analysis setup\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\
#create a version of the GWAMA output that has columns for chromosome and position\
input=no_sb.out\
cat $\{input\} | awk '\{print $0 "\\t" $1\}' >> temp.txt\
awk -F'_' -v OFS='\\t' '\{sub(/^[_ ]/,"");$18=$18\}1' temp.txt >> no_sb_chr_pos.txt\
#note that the chr column is called rs and the pos column is called number\
\
#generate snplist files for each signal to be used in colocalization analysis\
input=no_sb_chr_pos.txt\
\
mkdir 2mb_snplist\
\
while read p; do \
chr=$(echo $\{p\} | cut -d ',' -f2)\
pos=$(echo $\{p\} | cut -d ',' -f3)\
snp=$(echo $\{p\} | cut -d ',' -f1)\
cat $\{input\} | awk -F ' ' -v var1=$\{chr\} -v var2=$\{pos\} '\{OFS='\\t' ; if ($2==var1 && $3<var2+1000000 && $3>var2-1000000) print\}' | awk '!seen[$1]++' > ./2mb_snplist/$\{snp\}\
echo $\{snp\}\
done < sentinels.txt\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
D. GTEx eQTL colocalization analysis\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\
#identify relevant tissues for each signal\
rm ./tempfile\
\
awk -F ',' '\{print $1, $2"_"$3\}' /scratch16/abattle4/rebecca/final/gwama/sentinels.txt | awk '!seen[$1]++' > sentinel_chrpos\
\
while read p; do\
snp=$(echo $\{p\} | cut -d' ' -f2)\
for file in /scratch16/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt\
do\
tissue=$(echo $\{file\} | cut -d'/' -f8 | cut -d'.' -f1)\
grep 'chr'$\{snp\}'_' $\{file\} | awk -v var=$\{tissue\} '\{OFS="\\t"; print $0, var\}' >> tempfile\
done\
awk -v var=$\{snp\} '\{OFS="\\t"; print $2, var, $13\}' tempfile >> sentinel_tissues.txt\
rm tempfile\
done < sentinel_chrpos\
\
Rscript gtex_eqtl_coloc.r\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
E. GTEx sQTL colocalization analysis\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\
#identify relevant tissues for each signal\
rm ./tempfile\
\
awk -F ',' '\{print $1, $2"_"$3\}' /scratch16/abattle4/rebecca/final/gwama/sentinels.txt | awk '!seen[$1]++' > sentinel_chrpos\
\
while read p; do\
snp=$(echo $\{p\} | cut -d' ' -f2)\
for file in /scratch16/abattle4/lab_data/GTEx_v8/sqtl/GTEx_Analysis_v8_sQTL/*.v8.sqtl_signifpairs.txt\
do\
tissue=$(echo $\{file\} | cut -d'/' -f8 | cut -d'.' -f1)\
grep 'chr'$\{snp\}'_' $\{file\} | awk -v var=$\{tissue\} '\{OFS="\\t"; print $0, var\}' >> tempfile\
done\
awk -v var=$\{snp\} '\{OFS="\\t"; print $2, var, $13\}' tempfile >> sentinel_tissues.txt\
rm tempfile\
done < sentinel_chrpos\
\
Rscript gtex_sqtl_coloc.r\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
F. DICE eQTL colocalization analysis\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\
#LiftOver summary statistics to hg19\
\
#generate snplists for each meta-analysis signal (as in B)\
mkdir ./snplist\
\
input=no_sam_braz_hg19.txt\
\
while read p; do \
chr=$(echo $\{p\} | cut -d ',' -f19)\
pos=$(echo $\{p\} | cut -d ',' -f2)\
snp=$(echo $\{p\} | cut -d ',' -f1)\
cat $\{input\} | awk -F ' ' -v var1='chr'$\{chr\} -v var2=$\{pos\} '\{OFS='\\t' ; if ($19==var1 && $2<var2+500000 && $2>var2-500000) print\}' | awk '!seen[$1]++' > ./snplist/$\{snp\}\
echo $\{snp\}\
done < sentinels.txt\
\
#make snplist files for each of these regions in each DICE cell type\
Rscript dice_signals.r\
\
#run colocalization analysis\
while read p; do \
snp=$(echo $\{p\} | cut -d ',' -f1)\
Rscript dice_coloc.r $\{snp\}\
done < sentinels.txt\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
G. eQTLGen eQTL colocalization analysis\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\
#use sentinels files from DICE analysis\
\
#make eQTLGen snplists for each signal\
mkdir eqtlgen_snplist\
input=/scratch16/abattle4/lab_data/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz\
\
while read p; do \
chr=$(echo $\{p\} | cut -f19 -d ',')\
pos=$(echo $\{p\} | cut -f2 -d ',')\
snp=$(echo $\{p\} | cut -f1 -d ',')\
zcat $\{input\} | awk -F ' ' -v var1=$\{chr\} -v var2=$\{pos\} '\{OFS='\\t' ; if ($3==var1 && $4<var2+1000000 && $4>var2-1000000) print\}' | awk '!seen[$1]++' > ./eqtlgen_snplist/$\{snp\}\
done < sentinels.txt\
\
#run colocalization analysis\
mkdir geno_QTL\
\
while read p; do \
snp=$(echo $\{p\} | cut -f1 -d ',')\
pos=$(echo $\{p\} | cut -f2 -d ',')\
Rscript eqtlgen_byvariant_coloc.r $\{snp\} $\{pos\}\
done < sentinels.txt\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
H. TWAS analysis with FUSION\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#follow documentation on \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 \outl0\strokewidth0 \strokec4 http://gusevlab.org/projects/fusion/\cf2 \kerning1\expnd0\expndtw0 \CocoaLigature0 \outl0\strokewidth0 \
#download precomputed model\
mkdir YFS_WEIGHTS\
cd YFS_WEIGHTS\
wget https://data.broadinstitute.org/alkesgroup/FUSION/WGT/YFS.BLOOD.RNAARR.tar.bz2\
tar xjf YFS.BLOOD.RNAARR.tar.bz2\
\
#convert summary statistics to sumstats format\
ml r/4.0.2\
Rscript bychr_sumstats.r\
\
#catenate files\
cat bychr_sumstats/no_sb_chrpos_1 > no_sb_chrpos_rsid.txt\
#iterate for other chromosomes skipping the header\
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do\
cat bychr_sumstats/no_sb_chrpos_$\{chr\} | tail -n +2 >> no_sb_chrpos_rsid.txt\
done\
\
cd /scratch16/abattle4/rebecca/final/twas/\
\
ml anaconda\
conda activate ldsc\
python /data/abattle4/rebecca/ldsc/munge_sumstats.py --sumstats no_sb_chrpos_rsid_whitespace.txt --merge-alleles /data/abattle4/rebecca/ldsc_celltype/w_hm3.snplist --out tl_meta --a1-inc --snp rsid --a1 other --a2 ref --p pvalue --frq eaf --N-col n_samples --chunksize 500000 --ignore snpID,chr_hg38pos,chr,pos,beta,se,beta_95L,beta_95U,z,-log10_pvalue,q_statistic,q_pvalue,i2,n_studies,effects\
\
conda deactivate\
\
#compress input file\
gzip no_sb_chrpos_rsid_whitespace.txt\
\
#run FUSION\
ml r/4.0.2\
cd /scratch16/abattle4/rebecca/final/twas/\
mkdir yfs_output\
\
#YFS_Weights\
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do\
Rscript /data/abattle4/rebecca/software/fusion_twas-master/FUSION.assoc_test.R \\\
--sumstats tl_meta_noblank.sumstats \\\
--weights /data/abattle4/rebecca/software/fusion_twas-master/YFS_WEIGHTS/YFS.BLOOD.RNAARR.pos \\\
--weights_dir /data/abattle4/rebecca/software/fusion_twas-master/YFS_WEIGHTS/ \\\
--ref_ld_chr /data/abattle4/rebecca/software/fusion_twas-master/LDREF/1000G.EUR. \\\
--chr $\{chr\} \\\
--out yfs_output/tl_meta_yfs_$\{chr\}.dat\
done\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
I. SuSiE 95% credible set analysis\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
ml r/4.0.2\
ml gsl/2.7\
ml parallel\
\
mkdir ld\
mkdir allsnps\
mkdir credset\
mkdir manplot\
mkdir process_cmd\
mkdir tmp\
\
rm process_cmd/run_susie.txt\
\
while read p; do\
chr=$(echo $\{p\} | cut -d ',' -f1)\
pos=$(echo $\{p\} | cut -d ',' -f2)\
chrpos=$(echo $\{p\} | cut -d ',' -f3)\
Rscript run_susie.r $\{chr\} $\{pos\} $\{chrpos\} >> process_cmd/run_susie.txt\
done < susie_coordinates.csv\
\
parallel --tmpdir /scratch16/abattle4/rebecca/final/susie/tmp :::: process_cmd/run_susie.txt\
\
rm ld/*\
rm snplist_*\
\
\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
J. Stratified LDSC\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#use GWAMA to meta-analyze the European telomere length GWAS from TOPMed and Li et al.\
~/work/marios/GWAMA/GWAMA -i /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/gwama_li_top_euro/gwama.in -o /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/gwama_li_top_euro/gwama_euro -qt --random\
\
#regresssion\
ml anaconda\
conda activate ldsc\
\
#sumstats\
python /work-zfs/abattle4/rebecca/ldsc/munge_sumstats.py --sumstats /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/gwama_li_top_euro/gwama_euro.out --merge-alleles /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/w_hm3.snplist --out gwama_euro --a1-inc --snp rs_numer --a1 reference_allele --a2 other_allele --p p-value --frq eaf --N-col n_samples --chunksize 500000\
\
#partitioned heritability\
python /work-zfs/abattle4/rebecca/ldsc/ldsc.py --h2 gwama_euro.sumstats.gz --ref-ld-chr /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/baseline/baseline. --w-ld-chr /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/1000G_frq/1000G.mac5eur. --out gwama_euro\
conda deactivate\
\
#without changing the directory I kept getting and error for the index file even though I'm 100% the paths are correct. Must be something underlying the python scripts.\
cd /work-zfs/abattle4/rebecca/ldsc_celltype\
\
#cell type specificity\
cts_name=Multi_tissue_gene_expr\
/work-zfs/abattle4/rebecca/ldsc/ldsc.py --h2-cts /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/gwama_li_top_euro/gwama_euro.sumstats.gz --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. --out gwama_euro_$\{cts_name\} --ref-ld-chr-cts $cts_name.ldcts --w-ld-chr weights_hm3_no_hla/weights.\
\
cts_name=Multi_tissue_chromatin\
/work-zfs/abattle4/rebecca/ldsc/ldsc.py --h2-cts /work-zfs/abattle4/rebecca/ldsc_ph_walkthrough/gwama_li_top_euro/gwama_euro.sumstats.gz --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. --out gwama_euro_$\{cts_name\} --ref-ld-chr-cts $cts_name.ldcts --w-ld-chr weights_hm3_no_hla/weights.\
\
conda deactivate\
\
ldsc_plot.r\
\
\
\
}