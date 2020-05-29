Association of DG4R to CNA
================

Analysis of the association of DG4R to CNA
------------------------------------------

The analysis involved different steps.

1.  shuffle DG4r acorss genome
2.  overlap DG4r with individual copy number alteration in the respective PDTXs (i.e. check overlap DG4 vs different type of CNA in each model)
3.  collect the overlaps of the actual data
4.  collect the overlaps obtained using randomized data
5.  print data to file where one file capture the same type of CNA (ex. AMP, GAIN, NEUT, etc... ) across all individual PDRX
6.  plot

``` bash

#shuffle DG4
##############

path_analysis=~/20200219_PDTX_NewOld
cases=(DG4R_CG4R)
g=~/reference_genomes/robert_genome/hg19_robert_github.size.genome
n=(1 2 3 4 5 6 7 8 9 10); 

for case in ${cases[@]}
do
  path_DG4R=$path_analysis/$case
  cd $path_DG4R
  path_DG4R_shuffled=$path_DG4R/shuffled_DG4R
  mkdir $path_DG4R_shuffled
  for f in *_vs_all_up.bed; 
  do
    for i in "${!n[@]}"; 
    do 
      sbatch --mem 6G --wrap "bedtools shuffle -chromFirst -noOverlapping -seed ${n[i]} -i $f -g $g | sort -k1,1 -k2,2n > $path_DG4R_shuffled/${f%%.bed}.shuffle_${n[i]}.bed" 
    done; 
  done
done

# get the list of PTDX names
path_DG4R=$path_analysis/$case
cd $path_DG4R
PDTX=`ls *_vs_all_up.bed | sed 's/_vs_all_up.bed//g' | sed 's/peak_counts.norm_filtered.detable.//g'`
echo $PDTX
#139M_ 143_ 201_ 316_284 331_284 521_284 5_317 79_ 863_317 9_317 98_ AB551 AB555 AB577 AB580 AB636 AB790 PAR1006_3 PAR1022 STG139 STG195 STG282

# loop over PDTX LIST AND OVERLAP THE RELATIVE PDTX OVER THE CNA REGIONS OF INTEREST
for case in ${cases[@]}
do
  path_DG4R=$path_analysis/$case
  cd $path_DG4R
  path_DG4R_shuffled=$path_DG4R/shuffled_DG4R
  cd $path_DG4R_shuffled
  for pdtx in ${PDTX[@]}
  do
    shuffled_sets=`ls *$pdtx*shuffle*bed`
    echo $pdtx
    
    pdtx_spec_AMP=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*AMP.100kb.merged.bed`
    pdtx_spec_GAIN=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*GAIN.100kb.merged.bed`
    pdtx_spec_NEUT=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*NEUT.100kb.merged.bed`
    pdtx_spec_HETD=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*HETD.100kb.merged.bed`
    pdtx_spec_HOMD=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*HOMD.100kb.merged.bed`
      
    echo $pdtx_spec_AMP
    echo $pdtx_spec_GAIN
    echo $pdtx_spec_NEUT
    echo $pdtx_spec_HETD
    echo $pdtx_spec_HOMD
      
    echo "================="
    #for f in ${shuffled_sets[@]}
    #do
      #sbatch --mem 2 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_NEUT > ${f%%.bed}.int.NEUT.100kb.bed"
      #sbatch --mem 2 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_AMP > ${f%%.bed}.int.AMP.100kb.bed"
      #sbatch --mem 2 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_GAIN > ${f%%.bed}.int.GAIN.100kb.bed"
      #sbatch --mem 2 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_HETD > ${f%%.bed}.int.HETD.100kb.bed"
      #sbatch --mem 2 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_HOMD > ${f%%.bed}.int.HOMD.100kb.bed"
  
    #done
    actual_case=`ls $path_DG4R/peak_counts.norm_filtered.detable.*$pdtx*_vs_all_up.bed`
    ls $actual_case
    f_actual_case="$(basename -- $actual_case)"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_NEUT > ${f_actual_case%%.bed}.int.NEUT.100kb.bed"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_AMP > ${f_actual_case%%.bed}.int.AMP.100kb.bed"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_GAIN > ${f_actual_case%%.bed}.int.GAIN.100kb.bed"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_HETD > ${f_actual_case%%.bed}.int.HETD.100kb.bed"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_HOMD > ${f_actual_case%%.bed}.int.HOMD.100kb.bed"
  done
done


# for each pdtx and each CN Alteration build the fold enrichment (FR)
for case in ${cases[@]}
do
  path_DG4R=$path_analysis/$case
  cd $path_DG4R
  PDTX=`ls *_vs_all_up.bed | sed 's/_vs_all_up.bed//g' | sed 's/peak_counts.norm_filtered.detable.//g'`
  
  path_DG4R_shuffled=$path_DG4R/shuffled_DG4R
  cd $path_DG4R_shuffled
  CNA_type=(AMP GAIN NEUT HETD HOMD)
  for pdtx in ${PDTX[@]}
  do
    for cna in ${CNA_type[@]}
    do
      echo $pdtx 
      echo $cna
    #rob command
    #wc -l *${pdtx}*.int.$cna.100kb.bed  | head -n 11 | awk '{printf( "%s ", $1 ); } END {print "\n"}' | awk '{print $1/$2 "\t" $1/$3 "\t" $1/$4 "\t" $1/$5 "\t" $1/$6 "\t" $1/$7 "\t" $1/$8 "\t" $1/$9 "\t" $1/$10 "\t" $1/$11}' 
      wc -l *${pdtx}*int.$cna.100kb.bed |head -n 11|  awk '{printf( "%s ", $1 ); } END {print "\n"}' | head -n 1 | awk '{ if ($1==0 || $2==0 ) {print "0\t0\t0\t0\t0\t0\t0\t0\t\t0\t0"} else {print $1/$2 "\t" $1/$3 "\t" $1/$4 "\t" $1/$5 "\t" $1/$6 "\t" $1/$7 "\t" $1/$8 "\t" $1/$9 "\t" $1/$10 "\t" $1/$11}}' | tr '\t' '\n' > pasted_FR_${pdtx}.vs.all_int.${cna}.sites
    echo "=============="
    done
  done
done

CNA_sets=(AMP GAIN NEUT HETD HOMD)
for case in ${cases[@]}
do
  path_DG4R=$path_analysis/$case
  path_DG4R_shuffled=$path_DG4R/shuffled_DG4R
  cd $path_DG4R_shuffled
  for cna in ${CNA_sets[@]}
  do
    paste pasted*vs.all_*$cna.sites > PDTX_vs_all_$cna.sites
  done
done


# shuffle CG4r
# shuffle CRG4 10 tims
path_DG4R=~/20200219_PDTX_NewOld/DG4R_CG4R
path_DG4R_shuffled=~/20200219_PDTX_NewOld/DG4R_CG4R/shuffled_DG4R
g=~/reference_genomes/robert_genome/hg19_robert_github.size.genome
cd $path_DG4R
CG4=common.PDTX.hg19.peaks.sorted.bed  
n=(1 2 3 4 5 6 7 8 9 10); 
for i in "${!n[@]}"
do
  sbatch --mem 6G --wrap "bedtools shuffle -chromFirst -noOverlapping -seed ${n[i]} -i $CG4 -g $g | sort -k1,1 -k2,2n > $path_DG4R_shuffled/${CG4%%.bed}.shuffle_${n[i]}.bed" 
done

# overlap of CG4R with each individual tumor CNA regions
path_DG4R=~/20200219_PDTX_NewOld/DG4R_CG4R
cd $path_DG4R
PDTX=`ls *_vs_all_up.bed | sed 's/_vs_all_up.bed//g' | sed 's/peak_counts.norm_filtered.detable.//g'`
echo $PDTX
#005_317 009_317 139M 143 179 201 316_284 331_284 521_284 863_317 98 AB551 AB555 AB577 AB580 AB636 AB790 PAR1006 PAR1022 STG139 STG195 STG282

cd $path_DG4R_shuffled
common_G4_shuffled_sets=`ls common*shuffle*bed`
for pdtx in ${PDTX[@]}
do
  echo " ================ >>" $pdtx "<< ================"
  pdtx_spec_AMP=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*AMP.100kb.merged.bed`
  pdtx_spec_GAIN=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*GAIN.100kb.merged.bed`
  pdtx_spec_NEUT=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*NEUT.100kb.merged.bed`
  pdtx_spec_HETD=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*HETD.100kb.merged.bed`
  pdtx_spec_HOMD=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*HOMD.100kb.merged.bed`
    
  echo $pdtx_spec_AMP
  echo $pdtx_spec_GAIN
  echo $pdtx_spec_NEUT
  echo $pdtx_spec_HETD
  echo $pdtx_spec_HOMD
    
  #for f in ${common_G4_shuffled_sets[@]}
  #do
  #  
  #  echo "================="
  #  sbatch --mem 5 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_NEUT > ${f%%.bed}.int.${pdtx}.NEUT.100kb.bed"
  #  sbatch --mem 5 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_AMP > ${f%%.bed}.int.${pdtx}.AMP.100kb.bed"
  #  sbatch --mem 5 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_GAIN > ${f%%.bed}.int.${pdtx}.GAIN.100kb.bed"
  #  sbatch --mem 5 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_HETD > ${f%%.bed}.int.${pdtx}.HETD.100kb.bed"
  #  sbatch --mem 5 --wrap "bedtools intersect -wa -a $f -b $pdtx_spec_HOMD > ${f%%.bed}.int.${pdtx}.HOMD.100kb.bed"
#
#  done
  actual_case=$path_DG4R/common.PDTX.hg19.peaks.bed
  f_actual_case="$(basename -- $actual_case)"
  sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_NEUT > ${f_actual_case%%.bed}.int.${pdtx}.NEUT.100kb.bed"
  sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_AMP > ${f_actual_case%%.bed}.int.${pdtx}.AMP.100kb.bed"
  sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_GAIN > ${f_actual_case%%.bed}.int.${pdtx}.GAIN.100kb.bed"
  sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_HETD > ${f_actual_case%%.bed}.int.${pdtx}.HETD.100kb.bed"
  sbatch --mem 2G --wrap "bedtools intersect -wa -a $actual_case -b $pdtx_spec_HOMD > ${f_actual_case%%.bed}.int.${pdtx}.HOMD.100kb.bed"
done


path_DG4R=~/20200219_PDTX_NewOld/DG4R_CG4R
path_DG4R_shuffled=~/20200219_PDTX_NewOld/DG4R_CG4R/shuffled_DG4R
cd $path_DG4R
PDTX=`ls *_vs_all_up.bed | sed 's/_vs_all_up.bed//g' | sed 's/peak_counts.norm_filtered.detable.//g'`
cd $path_DG4R_shuffled
CNA_type=(AMP GAIN NEUT HETD HOMD)
for pdtx in ${PDTX[@]}
do
  for cna in ${CNA_type[@]}
  do
    echo $pdtx 
    echo $cna
    #rob command
    #wc -l common*${pdtx}*.int.$cna.100kb.bed  | head -n 11 | awk '{printf( "%s ", $1 ); } END {print "\n"}' | awk '{print $1/$2 "\t" $1/$3 "\t" $1/$4 "\t" $1/$5 "\t" $1/$6 "\t" $1/$7 "\t" $1/$8 "\t" $1/$9 "\t" $1/$10 "\t" $1/$11}' 
    wc -l common*int.*${pdtx}*.$cna.100kb.bed |head -n 11|  awk '{printf( "%s ", $1 ); } END {print "\n"}' | head -n 1 | awk '{ if ($1==0 || $2==0 ) {print "0\t0\t0\t0\t0\t0\t0\t0\t0\t0"} else {print $1/$2 "\t" $1/$3 "\t" $1/$4 "\t" $1/$5 "\t" $1/$6 "\t" $1/$7 "\t" $1/$8 "\t" $1/$9 "\t" $1/$10 "\t" $1/$11}}' | tr '\t' '\n'  >  common_pasted_FR_${pdtx}.vs.all_int.${cna}.sites
    echo "=============="
  done
done
```

Use the data to produce plot. Check before that all \*sites files have 0 when value is missing

``` r
#' @title association ΔG4r to CNA regions
#' @description import process and visualzie fold-enrichments of ΔG4r at CNA regions for each PDTX
#' @param frequencies of observed (actual case) overlaps of ΔG4r at CNA regions and frequencies of obaserveved (random cases) overlaps of ΔG4r at CNA regions
#' @return boxplot and summary txt file
#' @author angela simeone  \email{angela.simeone@cruk.cam.ac.uk}


library(reshape)
library(ggplot2)
library(forcats)
library(gridExtra)

setwd('~/PDTX_dir/PDTX_dm6_inputSubtraction_/DG4R/shuffled_DG4R')
AMP <- read.table('PDTX_vs_all_AMP.sites',stringsAsFactors = F)
GAIN <- read.table('PDTX_vs_all_GAIN.sites',stringsAsFactors = F)
NEUT <- read.table('PDTX_vs_all_NEUT.sites',stringsAsFactors = F)
HETD <- read.table('PDTX_vs_all_HETD.sites',stringsAsFactors = F)
HOMD <- read.table('PDTX_vs_all_HOMD.sites',stringsAsFactors = F)
colnames(AMP) <- colnames(GAIN) <- colnames(NEUT) <- colnames(HETD) <- colnames(HOMD) <-  c("139M_","143_","201_","316_284","331_284","521_284","5_317","179_","863_317","9_317","98_","AB551","AB555","AB577","AB580","AB636","AB790","PAR1006_3","PAR1022","STG139","STG195","STG282")

AMP_df <- melt(as.data.frame(AMP))
dg4r_AMP <- ggplot(data=AMP_df, aes(variable, value)) + geom_boxplot() + ylim(0,7)+ ggtitle('DG4r AMP')+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

GAIN_df <- melt(as.data.frame(GAIN))
dg4r_GAIN <-ggplot(data=GAIN_df, aes(variable, value)) + geom_boxplot() + ylim(0,7)+ ggtitle('DG4r GAIN')+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

NEUT_df <- melt(as.data.frame(NEUT))
dg4r_NEUT <- ggplot(data=NEUT_df, aes(variable, value)) + geom_boxplot() + ylim(0,7)+ ggtitle('DG4r NEUT')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

HETD_df <- melt(as.data.frame(HETD))
dg4r_HETD <- ggplot(data=HETD_df, aes(variable, value)) + geom_boxplot() + ylim(0,7)+ ggtitle('DG4r HETD')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

HOMD_df <- melt(as.data.frame(HOMD))
dg4r_HOMD <- ggplot(data=HOMD_df, aes(variable, value)) + geom_boxplot() + ylim(0,7) + ggtitle('DG4r HOMD')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_AMP <- read.table('CG4R_PDTX_vs_all_AMP.sites',stringsAsFactors = F)
CG4_GAIN <- read.table('CG4R_PDTX_vs_all_GAIN.sites',stringsAsFactors = F)
CG4_NEUT <- read.table('CG4R_PDTX_vs_all_NEUT.sites',stringsAsFactors = F)
CG4_HETD <- read.table('CG4R_PDTX_vs_all_HETD.sites',stringsAsFactors = F)
CG4_HOMD <- read.table('CG4R_PDTX_vs_all_HOMD.sites',stringsAsFactors = F)
colnames(CG4_AMP) <- colnames(CG4_GAIN) <- colnames(CG4_NEUT) <- colnames(CG4_HETD) <- colnames(CG4_HOMD) <-  c('005_317','009_317','139M','143','179','201','316_284','331_284','521_284','863_317','98','AB551','AB555','AB577','AB580','AB636','AB790','PAR1006','PAR1022','STG139','STG195','STG282')

CG4_AMP_df <- melt(as.data.frame(CG4_AMP))
cg4r_AMP <- ggplot(data=CG4_AMP_df, aes(variable, value)) + geom_boxplot() + ylim(0,7)+ ggtitle('CG4r AMP')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_GAIN_df <- melt(as.data.frame(CG4_GAIN))
cg4r_GAIN <-ggplot(data=CG4_GAIN_df, aes(variable, value)) + geom_boxplot() + ylim(0,7)+ ggtitle('CG4r GAIN')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_NEUT_df <- melt(as.data.frame(CG4_NEUT))
cg4r_NEUT <-ggplot(data=CG4_NEUT_df, aes(variable, value)) + geom_boxplot() + ylim(0,7)+ ggtitle('CG4r NEUT')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_HETD_df <- melt(as.data.frame(CG4_HETD))
cg4r_HETD <-ggplot(data=CG4_HETD_df, aes(variable, value)) + geom_boxplot() + ylim(0,7)+ ggtitle('CG4r HETD')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_HOMD_df <- melt(as.data.frame(CG4_HOMD))
cg4r_HOMD <-ggplot(data=CG4_HOMD_df, aes(variable, value)) + geom_boxplot() + ylim(0,7) + ggtitle('CG4r HOMD')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)


#GG <- grid.arrange(dg4r_HOMD, dg4r_HETD,dg4r_NEUT,dg4r_GAIN,dg4r_AMP,cg4r_HOMD, cg4r_HETD,cg4r_NEUT,cg4r_GAIN,cg4r_AMP,nrow=2)
GG <- grid.arrange(dg4r_HOMD, dg4r_HETD,dg4r_NEUT,dg4r_GAIN,dg4r_AMP,nrow =1)
ggsave(filename = 'PDTX_relationthip_to_CNA.pdf',plot=GG, dpi = 150,width=12,height=8)

CG4_AMP <- read.table('CG4R_PDTX_vs_all_AMP.sites',stringsAsFactors = F)
CG4_GAIN <- read.table('CG4R_PDTX_vs_all_GAIN.sites',stringsAsFactors = F)
CG4_NEUT <- read.table('CG4R_PDTX_vs_all_NEUT.sites',stringsAsFactors = F)
CG4_HETD <- read.table('CG4R_PDTX_vs_all_HETD.sites',stringsAsFactors = F)
CG4_HOMD <- read.table('CG4R_PDTX_vs_all_HOMD.sites',stringsAsFactors = F)
colnames(CG4_AMP) <- colnames(CG4_GAIN) <- colnames(CG4_NEUT) <- colnames(CG4_HETD) <- colnames(CG4_HOMD) <-  c('005_317','009_317','139M','143','179','201','316_284','331_284','521_284','863_317','98','AB551','AB555','AB577','AB580','AB636','AB790','PAR1006','PAR1022','STG139','STG195','STG282')

CG4_AMP_df <- melt(as.data.frame(CG4_AMP))
cg4r_AMP <- ggplot(data=CG4_AMP_df, aes(variable, value)) + geom_boxplot() + ylim(0,0.5)+ ggtitle('CG4r AMP')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_GAIN_df <- melt(as.data.frame(CG4_GAIN))
cg4r_GAIN <-ggplot(data=CG4_GAIN_df, aes(variable, value)) + geom_boxplot() + ylim(0,0.5)+ ggtitle('CG4r GAIN')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_NEUT_df <- melt(as.data.frame(CG4_NEUT))
cg4r_NEUT <-ggplot(data=CG4_NEUT_df, aes(variable, value)) + geom_boxplot() + ylim(0,0.5)+ ggtitle('CG4r NEUT')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_HETD_df <- melt(as.data.frame(CG4_HETD))
cg4r_HETD <-ggplot(data=CG4_HETD_df, aes(variable, value)) + geom_boxplot() + ylim(0,0.5)+ ggtitle('CG4r HETD')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)

CG4_HOMD_df <- melt(as.data.frame(CG4_HOMD))
cg4r_HOMD <-ggplot(data=CG4_HOMD_df, aes(variable, value)) + geom_boxplot() + ylim(0,0.5) + ggtitle('CG4r HOMD')+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ geom_hline(yintercept=20, linetype="dashed", color = "red", size=1)


#GG <- grid.arrange(dg4r_HOMD, dg4r_HETD,dg4r_NEUT,dg4r_GAIN,dg4r_AMP,cg4r_HOMD, cg4r_HETD,cg4r_NEUT,cg4r_GAIN,cg4r_AMP,nrow=2)

GG_CG4r <- grid.arrange(cg4r_HOMD, cg4r_HETD, cg4r_NEUT,cg4r_GAIN,cg4r_AMP,nrow=1)
ggsave(filename = 'DG4_CG4r_relationship_to_CNA.pdf',plot=GG, dpi = 150,width=12,height=8)
ggsave(filename = 'CG4r_relationship_to_CNA.pdf',plot=GG_CG4r, dpi = 150,width=12,height=8)


## for figure 4
inds_rob <- c(grep("521",colnames(AMP)),grep("5_317",colnames(AMP)),grep("9_317",colnames(AMP)),grep("139M_",colnames(AMP)),grep("143_",colnames(AMP)),grep("201_",colnames(AMP)),grep("316_284",colnames(AMP)),grep("331_284",colnames(AMP)),grep("98_",colnames(AMP)))
library(robustbase)

AMP_median <- colMedians(data.matrix(AMP[,inds_rob]))
GAIN_median <- colMedians(data.matrix(GAIN[,inds_rob]))
NEUT_median <- colMedians(data.matrix(NEUT[,inds_rob]))
HETD_median <- colMedians(data.matrix(HETD[,inds_rob]))

write.table(t(AMP_median), file="selected_AMP_median.txt",quote = F, col.names = NA,sep="\t")
write.table(t(GAIN_median), file="selected_GAIN_median.txt",quote = F, col.names = NA,sep="\t")
write.table(t(NEUT_median), file="selected_NEUT_median.txt",quote = F, col.names = NA,sep="\t")
write.table(t(HETD_median), file="selected_HETD_median.txt",quote = F, col.names = NA,sep="\t")
```
