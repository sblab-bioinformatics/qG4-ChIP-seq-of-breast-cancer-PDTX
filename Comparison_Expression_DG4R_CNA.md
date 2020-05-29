Association of G4 signal to CNA and expression
================
angela simeone

Here we want to assess the association G4 signal to CNA and expression in each tumor samples.

To do so, we had to prepare 3 differnet types of data: PART1. Extract the G4 signal for promoters overlapping DG4 regions PART2. Extract information about pronomters and memebership to each of the CNA regions PART3. extract information about promoter and expression groups

PART1: extract G4 signal
------------------------

``` bash

promoter=~/PDTX/Annotation_robert/hg19.gene_name.promoters.bed
out_folder=~/PDTX/PDTX_dm6_inputSubtraction_/DG4_expression_CNA

for file in *lcpm.bed
do
  intersectBed -a $promoter -b $file -wa -wb > $out_folder/Promoter_BG4signal_${file%%.bed}.txt
done
```

Read the ouput in R and combine (by median) the BG4 signal (it is a logCPM in various peaks... )

``` r
setwd('~/PDTX/PDTX_dm6_inputSubtraction_/DG4_expression_CNA')
library(tidyverse)
list_files <- list.files(pattern='up.lcpm.txt')

for (i in (1:length(list_files))) {
  
 temp_file <- read.delim(list_files[i],sep = "\t",stringsAsFactors = F,header = F)
 temp_file$promoter <- paste0(temp_file$V1,"_",temp_file$V2,"_",temp_file$V3,"_",temp_file$V4)
 
 combined_data <- temp_file %>% dplyr::group_by(promoter) %>% summarise(median_G4=median(V8,na.rm = T))
 combined_data_to_print <- cbind(str_split_fixed(combined_data$promoter,"_",4),combined_data$median_G4)
 
 file_out <- gsub('up.lcpm.txt','up.lcpm.median.txt',list_files[i])
 write.table(combined_data_to_print,file=file_out,quote=F, sep = "\t", col.names = F, row.names = F)
}
```

PART2: extract promoters overlapping each of the CNA regions
------------------------------------------------------------

``` bash
path_analysis=~/PDTX_dir
cases=(DG4R_CG4R)
g=/scratcha/sblab/simeon01/reference_genomes/robert_genome/hg19_robert_github.size.genome
promoters=~/PDTX_dir/hg19.gene_name.promoters.bed

#path_DG4R=$path_analysis/$case
#cd $path_DG4R


for case in ${cases[@]}
do
  path_DG4R=$path_analysis/$case
  cd $path_DG4R
  PDTX=`ls *_vs_all_up.bed | sed 's/_vs_all_up.bed//g' | sed 's/peak_counts.norm_filtered.detable.//g'`
  echo $PDTX
  path_DG4R_shuffled=$path_DG4R/shuffled_DG4R
  cd $path_DG4R_shuffled
  for pdtx in ${PDTX[@]}
  do
      pdtx_spec_AMP=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*AMP.100kb.merged.bed`
      pdtx_spec_GAIN=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*GAIN.100kb.merged.bed`
      pdtx_spec_NEUT=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*NEUT.100kb.merged.bed`
      pdtx_spec_HETD=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*HETD.100kb.merged.bed`
      pdtx_spec_HOMD=`ls ~/PDTX_DG4R_CG4R/20191216_PDTX_bam_5M_CNA/*$pdtx*HOMD.100kb.merged.bed`
      
    echo $pdtx
    # for each pdtx levels, read .tab file, extract bed format including the
    actual_case=`ls $path_DG4R/*$pdtx*_vs_all_up.bed`
    f_actual_case="$(basename -- $actual_case)"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $promoters -b $actual_case| sortBed -i - | bedtools intersect -a - -b $pdtx_spec_NEUT -wa  > ./promoter_analysis/${f_actual_case%%.bed}.promoters.NEUT.100kb.bed"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $promoters -b $actual_case| sortBed -i - | bedtools intersect -a - -b $pdtx_spec_AMP -wa > ./promoter_analysis/${f_actual_case%%.bed}.promoters.AMP.100kb.bed"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $promoters -b $actual_case| sortBed -i - | bedtools intersect -a - -b $pdtx_spec_GAIN -wa > ./promoter_analysis/${f_actual_case%%.bed}.promoters.GAIN.100kb.bed"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $promoters -b $actual_case| sortBed -i - | bedtools intersect -a - -b $pdtx_spec_HETD -wa > ./promoter_analysis/${f_actual_case%%.bed}.promoters.HETD.100kb.bed"
    sbatch --mem 2G --wrap "bedtools intersect -wa -a $promoters -b $actual_case| sortBed -i - | bedtools intersect -a - -b $pdtx_spec_HOMD -wa > ./promoter_analysis/${f_actual_case%%.bed}.promoters.HOMD.100kb.bed"
    echo " *** pdtx finished *** "
    echo " *** =============== *** "
  done
done

CNA=(NEUT AMP GAIN HETD HOMD)
# when it is finished for each case combine the promoters - CNA by PDTX into a single file
for case in ${cases[@]}
do
  path_DG4R=$path_analysis/$case
  cd $path_DG4R
  PDTX=`ls *_vs_all_up.bed | sed 's/_vs_all_up.bed//g' | sed 's/peak_counts.norm_filtered.detable.//g'`
  echo $PDTX
  path_DG4R_shuffled=$path_DG4R/shuffled_DG4R/promoter_analysis
  cd $path_DG4R_shuffled
  for pdtx in ${PDTX[@]}
  do
    echo $pdtx
    # for each pdtx levels, read .tab file, extract bed format including the
    actual_case=`ls ./*$pdtx*.bed`
    out_summary_cna_file_pdtx=${pdtx}.promoters.allCNA.100kb.txt
    touch $out_summary_cna_file_pdtx
    for cna in ${CNA[@]}
    do
      temp_file=`ls *$pdtx*.promoters.${cna}.100kb.bed`
      echo $temp_file
      wc -l  $temp_file
      if [ -s "$temp_file" ]
      then
        awk -v cna_case="$cna" '{print $0"\t"cna_case}' $temp_file>> $out_summary_cna_file_pdtx
      fi
    done
  done
done
```

PART3: expression
-----------------

For this, I can re-use what develped before

Refer to the R script that combines together the stuff computed in PART1,PART2,PART3

``` bash
promoter=~/PDTX/Annotation_robert/hg19.gene_name.promoters.bed
cd ~/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction
mkdir promoter_lfc
for file in *vs_all_up.lcpm.bed
do
  intersectBed -a $promoter -b $file -wa -wb > ./promoter_lfc/${file%%.lcpm.bed}.lcpm.promoter.bed
done
```

Refer to the R script [**Analysis\_promoters\_G4levels\_expression\_CNA**](./Analysis_promoter_G4levels_expression_CNA.R) that combines all data previously prepared.
