Characterization of DG4R and CG4R
================

This document illustrates the characterization of DG4r previously identified.

Pairwise comparison of PDTXs DG4r by computing Jaccard index
------------------------------------------------------------

``` bash

## == cross compare DG4R
mkdir intervene_DG4R
mkdir intervene_DG4R_jaccard
#echo sbatch --mem 4G --wrap "intervene venn -i *up.bed --filenames -o $cell_line_path$case/DG4R/intervene_DG4R"
sbatch --mem 4G --wrap "intervene pairwise  -i *up.bed --filenames --compute frac --htype pie -o $path_DG4R_CG4R/intervene_DG4R"
sbatch --mem 4G --wrap "intervene pairwise  --compute=jaccard --genome hg19 --fontsize 3  -i *up.bed --filenames --htype pie -o $path_DG4R_CG4R/intervene_DG4R_jaccard"
```

characterize overlap between PDTX DG4r and CG4r and breast cancer drivers ( -- Intervene -- )
---------------------------------------------------------------------------------------------

``` bash
pdtx_path=/scratchb/sblab/simeon01/20200219_PDTX_NewOld
breast_cancer_driver=/scratcha/sblab/simeon01/reference_genomes/robert_genome/robert_annotations_paper/regions_driver_breastcancer.hg19.sort.bed
cases=(DG4R_CG4R)
for case in ${cases[@]}
do
  cd $pdtx_path/$case/DG4R
  out_folder=$pdtx_path/$case/DG4R/intervene_DG4r_CG4r_breast_cancer_drivers/
  common_CG4=$pdtx_path/$case/common.PDTX.hg19.peaks.sorted.bed
  sbatch --mem 6G --wrap "intervene pairwise --compute=fisher --genome hg19 -i $breast_cancer_driver $common_CG4 $pdtx_path/$case/*_vs_all_up.bed -o $out_folder"
done
```

characterize overlap between PDTX DG4r and CG4r and breast cancer drivers ( -- gat -- )
---------------------------------------------------------------------------------------

Fold enrichments of PDTX DG4r at 45 breast cancer drivers are estimated using using GAT (Genomic Association Tester).

<https://readthedocs.org/projects/gat/downloads/pdf/latest/>

``` bash

#prepare annotation file for GAT 

awk '{print $0"\tBREAST_DRIV"}' /scratcha/sblab/simeon01/reference_genomes/robert_genome/robert_annotations_paper/regions_driver_breastcancer.hg19.sort.bed > /scratcha/sblab/simeon01/reference_genomes/robert_genome/robert_annotations_paper/annotation_regions_driver_breastcancer.hg19.sort.bed

workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
anno=/scratcha/sblab/simeon01/reference_genomes/robert_genome/robert_annotations_paper/annotation_regions_driver_breastcancer.hg19.sort.bed
for case in ${cases[@]}
do
  cd $pdtx_path/$case
  out_folder=$pdtx_path/$case/DG4R/intervene_DG4r_CG4r_breast_cancer_drivers/
  for file in *bed
  do
    echo $file
    cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen} -n 10000 --annotations=$anno --segments=${file} -S $out_folder/${file%%.bed}_on_BREAST_CANC_genes.dat -t 4"
    echo ${file%%.bed}_${anno%%.bed}.dat
    echo $cmd_gat
    sbatch --mem 3G --wrap "$cmd_gat"
    echo " "
  done
done

# collect all fold enrichments from .dat file and assemple a single file with all fold enrichments for each PDTX. The plot will report the fold enrichment of the PDTX peaks over the regions annotated as cacer drivers
for case in ${cases[@]}
do
  out_folder=$pdtx_path/$case/DG4R/intervene_DG4r_CG4r_breast_cancer_drivers
  cd $out_folder
  out_folder=$pdtx_path/$case/DG4R/intervene_DG4r_CG4r_breast_cancer_drivers
  touch $out_folder/fold_change_over_random_for_BREAST_CANC_genes.txt
  for file in *dat
  do
    m=${file##peak_counts.norm_filtered.detable.}
    f=${m%%_vs_all_up_on_BREAST_CANC_genes.dat}
    fc=`cat  $file | cut -f 8,10 | tail -n 1`
    echo $f$'\t'$fc >> $out_folder/fold_change_over_random_for_BREAST_CANC_genes.txt
  done
done
```

In R load the fold changes and plot a barplot

``` r
library(ggplot2)
setwd("/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction/DG4R/intervene_DG4r_CG4r_breast_cancer_drivers")
fc_table <- read.table('fold_change_over_random_for_BREAST_CANC_genes.txt', stringsAsFactors = F)
colnames(fc_table) <- c('PDTX','fc_over_random','log10_pvalue')
fc_table$log10_pvalue <- -log10(fc_table$log10_pvalue)

pdf('fold_enrichment_over_random_45_cancer_drivers.pdf')
ggplot(fc_table,aes(x=reorder(PDTX, fc_over_random),y=fc_over_random,fill=log10_pvalue)) + geom_bar(stat="identity") + geom_hline(yintercept=1, linetype="dashed", color = "red", size=1)+ xlab('PDTX') + ylab('fold enrichment G4 \n at breast cancer drivers over random') + ggtitle('45 common breast cancer drivers')+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()
dev.off()
```

G4 motif analysis
-----------------

After identifications of DG4R, generate CG4R. analyze the G4 motifs prevalence and enrichments at DG4R regions. Steps:

-   randomize DG4R
-   extract fasta from actual case and random cases
-   use R function that reads fasta and compute prevalence and enrichments
-   generate plots with prevalence and enrichments (barplots and pie charts).

``` bash

# CG4R generation
cd /Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction
subtractBed -A -a ../hg19_old_new_q005.all_peaks.25M.over99nt.sorted.bed -b *_vs_all_up.bed | sort -k1,1 -k2,2 > common.PDTX.hg19.peaks.bed

# go to the folder with the bed file of PDTX DG4R
Rscript /Users/simeon01/Documents/PDTX/Breast_cell_lines/breast_cells_drosophila_qBG4/barplot_N_peaks.R /Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction

# after moving the files on the cluster
### generate folder with DG4R and produce permuations of the same file
consensus_hg19=/scratchb/sblab/simeon01/20191113_winnie_katie_oct19_nov19_PDTX/fastq/trimmed/aligned/merged_lines/consensus_hg19/hg19_old_new_q005.all_peaks.25M.over99nt.sorted.bed

path_DG4R_CG4R=/scratchb/sblab/simeon01/20200219_PDTX_NewOld/DG4R_CG4R

g=/scratcha/sblab/simeon01/reference_genomes/robert_genome/hg19.size.genome
hg19_white=/scratcha/sblab/simeon01/reference_genomes/robert_genome/hg19.whitelist.bed
fasta=/scratcha/sblab/simeon01/reference_genomes/hg19_gen/hsa.hg19.fa

seed_random=(1 2 3 4 5)
cd $path_DG4R_CG4R
cp $consensus_hg19 $path_DG4R_CG4R

## == shuffle common CG4R
for n in ${seed_random[@]}
  do
    # permute common
  bedtools shuffle -seed ${n} -i $path_DG4R_CG4R/common.PDTX.hg19.peaks.bed -g $g -incl $hg19_white > $path_DG4R_CG4R/common.PDTX.hg19.peaks_shuffle${n}.bed
done

## == shuffle common DG4R
for n in ${seed_random[@]}
do
  # permute consensus
  consensus_base=`basename $consensus_hg19`
  bedtools shuffle -seed ${n} -i $consensus_hg19 -g $g -incl $hg19_white > $path_DG4R_CG4R/${consensus_base%%.bed}_shuffle${n}.bed
done
```

``` bash
g=/scratcha/sblab/simeon01/reference_genomes/robert_genome/hg19.size.genome
fasta=/scratcha/sblab/simeon01/reference_genomes/hg19_gen/hsa.hg19.fa
cd /scratchb/sblab/simeon01/20200219_PDTX_NewOld/DG4R_CG4R

outout_fasta_dir=/scratchb/sblab/simeon01/20200219_PDTX_NewOld/DG4R_CG4R/hg19_old_new_q005_fasta
mkdir $outout_fasta_dir
all_peaks=hg19_old_new_q005.all_peaks.25M.over99nt.sorted.bed 

sbatch --mem 4G --wrap "bedtools getfasta -fi $fasta -bed ${all_peaks} -fo $outout_fasta_dir/${all_peaks%%.bed}.fa"

for file in hg19_old_new_q005.all_peaks.25M.over99nt.sorted_shuffle*.bed
do
  echo $file
  sbatch --mem 4G --wrap "bedtools getfasta -fi $fasta -bed ${file} -fo $outout_fasta_dir/${file%%.bed}.fa"
done
```

``` bash
# == Run locally the R script that perform the analysis for G4 specific structures
cd /Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction/DG4R/hg19_old_new_q005_fasta

Rscript ~/sequence_hits_analysis_ChIP_C.R hg19_old_new_q005.all_peaks.25M.over99nt.sorted.fa hg19_old_new_q005.all_peaks.25M.over99nt.sorted_shuffle1.fa hg19_old_new_q005.all_peaks.25M.over99nt.sorted_shuffle2.fa hg19_old_new_q005.all_peaks.25M.over99nt.sorted_shuffle3.fa
```

Generate plots after sequence motifs analysis (G4 motifs, loops etc...)

``` r
library(gridExtra)
library(tidyverse)
piePlot <- function(count, categories) {
    dat <- data.frame(count = count, category = categories)
    dat$fraction <- dat$count / sum(dat$count)
    dat$ymax <- cumsum(dat$fraction)
    dat$ymin <- c(0, head(dat$ymax, n = -1))
    dat$label <- factor(paste(dat$category, dat$count), levels = paste(dat$category, dat$count))
    plot <-
        ggplot(dat, aes(
            fill = label, # fill by label not category
            ymax = ymax,
            ymin = ymin,
            xmin = 0,
            xmax = 1
        )) +
        geom_rect() +
        coord_polar(theta = "y") +
        theme(legend.position="top") + theme_void() # no need for labels anymore
    plot
}

generate_structure_plots <- function(dataframe_peaks,outuput_file,outuput_file_pie,outuput_file_fe, label_to_print,labels_legend){
  
  #pie chart
  pie_all_peaks <- ggplot(dataframe_peaks, aes(x = "", y = actual_count, fill = motif_name)) +
  geom_bar(width = 1, stat = "identity", 
           color = "white") +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values=c("skyblue4","tomato3","darkseagreen3","mediumpurple1","lightskyblue","coral","azure3"),
                    label=labels_legend)+
  #geom_text(aes(label = actual_count), position = position_stack(vjust = 0.5)) +
  ggtitle(label_to_print) + 
  theme(legend.position ="rigth") +
  theme_void()
  
  bar_fold_enrichments <- ggplot(dataframe_peaks, 
       aes(x = motif_name, y = fold_change_over_random, fill = motif_name)) +
  geom_bar(stat = "identity", color = "white",
           width = 0.7) + 
  scale_fill_manual(values=c("skyblue4","tomato3","darkseagreen3","mediumpurple1","lightskyblue","coral","azure3"),
                    label=label_to_print)+
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.8) +
  ggtitle(paste0("ESC - fold_enrich")) + ylab('fold_enrichm') + 
  theme(legend.position = "none",
        text = element_text(size=12),
        #labels=labels_legend,
        axis.text.x = element_text(angle = 90))
  
  Final_plot <- grid.arrange(pie_all_peaks,bar_fold_enrichments)
  ggsave(file=outuput_file,Final_plot)
  ggsave(file=outuput_file_pie,pie_all_peaks)
  ggsave(file=outuput_file_fe,bar_fold_enrichments)
  rm(Final_plot)
}
setwd('/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction/DG4R/hg19_old_new_q005_fasta')
assign('hg19_all_peaks_motifs', get(load('hg19_old_new_q005.all_peaks.25M.over99nt.sorted.results_motifs.Rdata')))

hg19_all_peaks_motifs <- hg19_all_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))

labels_legend_hg19_all_peaks_motifs <- paste0(hg19_all_peaks_motifs$actual_count,' (',hg19_all_peaks_motifs$percentage,'%)')


generate_structure_plots(hg19_all_peaks_motifs,
                         'hg19_old_new_q005.all_peaks.25M.over99nt_all_peaks_motifs.pdf',
                         'hg19_old_new_q005.all_peaks.25M.over99nt_all_peaks_motifs_pie_counts.pdf',
                         'hg19_old_new_q005.all_peaks.25M.over99nt_all_peaks_motifs_fold_enrich.pdf',
                         'hg19_old_new_q005.all_peaks.25M.over99nt_all_peaks',labels_legend_hg19_all_peaks_motifs)
write.table(results_motifs,file='hg19_old_new_q005.all_peaks.25M.over99nt.sorted.results_motifs.csv',sep = ",",row.names =F, col.names = T)
```

Prepare data to check association with expression
-------------------------------------------------

The following step is to annotate promoters by PDTXs' DG4R.

The R scripts [**Association DG4r to expression groups**](./Compare_expressionValues_DG4R_CG4R.R) perform the analyisis on the association of DG4r to highly expressed genes by PDTX.

``` bash
# run this locally

promoter=/Users/simeon01/Documents/PDTX/Annotation_robert/hg19.gene_name.promoters.bed
path_analysis=/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation
cases=(PDTX_only_dm6_inputSubtraction)

#1.  first annotate peaks across 
cd $path_analysis
for case in ${cases[@]}
  do
  cd $path_analysis/$case
  sortBed -i common.PDTX.hg19.peaks.bed > common.PDTX.hg19.peaks.sorted.bed
  ls *up.bed > list_pdtx.txt
  a=`cat list_pdtx.txt`
  bedtools annotate -i $promoter -files *up.bed -counts -names $a > hg19.gene_name.promoters.DG4r.bed
  bedtools annotate -i $promoter -files common.PDTX.hg19.peaks.sorted.bed -counts -names common > hg19.gene_name.promoters.CG4r.bed
done
```

original: Analysis\_PDTX\_NO\_cells\_new\_old\_PDTX.\*
