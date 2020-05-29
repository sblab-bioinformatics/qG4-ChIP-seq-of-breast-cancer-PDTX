#' @title compare ΔG4r and expression levels
#' @description compare ΔG4r and the expression levels of genes having them in promoters
#' @param promoters list, expression levels, list of promoters overlapping at least one ΔG4r
#' @return barplot and txt files with frequency of genes with ΔG4r stratified by expression level (high-med-low)
#' @author angela simeone  \email{angela.simeone@cruk.cam.ac.uk}


library(tidyverse)
library(Hmisc)
library(data.table)

# a=`cat list_pdtx.txt`
# bedtools annotate -i hg19.gene_name.promoters.bed -files /Users/simeon01/Documents/PDTX/REVISIONS/Giovanni_script/clean_25M_union_norm5_dm6_only_analysis/*up.bed -counts -names $a > hg19.gene_name.promoters.DG4r_up.bed

# cd /Users/simeon01/Documents/PDTX/DG4R_regions
# ls *_specific_G4.bed > list_DG4r_ROBERT.txt
# b=`cat list_DG4r_ROBERT.txt`
# cd /Users/simeon01/Documents/PDTX/Annotation_robert
# bedtools annotate  -i hg19.gene_name.promoters.bed -files /Users/simeon01/Documents/PDTX/DG4R_regions/*_specific_G4.bed -counts -names $b > hg19.gene_name.promoters.DG4r_ROBERT_original.bed

## ==================== use this block to look at **NEW + OLD** PDTX ======================
# NEW + OLD PDTX

path_peaks <- '/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction/'
setwd(path_peaks)

list_peaks <-list.files(path=path_peaks,pattern = 'up.bed')
tot_N_peaks_all <- c()
for (i in (1:length(list_peaks))){
  con <- file(paste0(path_peaks,list_peaks[i])) 
  tot_N_peaks_all <- c(tot_N_peaks_all, length(readLines(con)))
  close(con)
}
names(tot_N_peaks_all) <- gsub('peak_counts.norm_filtered.detable.','',list_peaks)
names(tot_N_peaks_all) <- gsub('_vs_all_up.bed','',names(tot_N_peaks_all))

DG4_annotation_on_promoters <- read.delim(paste0(path_peaks,'hg19.gene_name.promoters.DG4r.bed'),stringsAsFactors = F, header = T)

# fix all values to 1 ( we want to count only one time the event of DG4 overlapping promoters)
# temp_DG4_annotation_on_promoters <- DG4_annotation_on_promoters[5:dim(DG4_annotation_on_promoters)[2]]
# temp_DG4_annotation_on_promoters[temp_DG4_annotation_on_promoters>1] <- 1
# DG4_annotation_on_promoters[5:dim(DG4_annotation_on_promoters)[2]] <- temp_DG4_annotation_on_promoters
# rm(temp_DG4_annotation_on_promoters)
# keep original names
orig_col_names_DG4 <- colnames(DG4_annotation_on_promoters)[5:dim(DG4_annotation_on_promoters)[2]]
CG4_annotation_on_promoters <- read.delim(paste0(path_peaks,'hg19.gene_name.promoters.CG4r.bed'),stringsAsFactors = F, header = T)
#colnames(DG4_annotation_on_promoters) <- gsub('PDTX_only_dm6_inputSubtraction_','',colnames(DG4_annotation_on_promoters))
colnames(DG4_annotation_on_promoters) <- gsub('PDTX_NormalCellsOnly_dm6_inputSubtraction_2stepNormalization_','',colnames(DG4_annotation_on_promoters))
colnames(DG4_annotation_on_promoters) <- gsub('__vs_all_up.bed','',colnames(DG4_annotation_on_promoters))
colnames(DG4_annotation_on_promoters) <- gsub('_vs_all_up.bed','',colnames(DG4_annotation_on_promoters))


## =====================================================================================


# ## ==================== use this block to look at **OLD** PDTX ============================
# #OLD PDTX
# path_peaks <- '/Users/simeon01/Documents/PDTX/DG4R_regions/'
# DG4_annotation_on_promoters <- read.delim('/Users/simeon01/Documents/PDTX/Annotation_robert/hg19.gene_name.promoters.DG4r_ROBERT_original.bed',stringsAsFactors = F, header = T)
# # keep original names
# orig_col_names_DG4 <- colnames(DG4_annotation_on_promoters)[5:dim(DG4_annotation_on_promoters)[2]]
# colnames(DG4_annotation_on_promoters) <- gsub('_specific_G4.bed','',colnames(DG4_annotation_on_promoters))
# ## =====================================================================================


# promoters infos +-1kb 
promoters <- read.table('/Users/simeon01/Documents/PDTX/Annotation_robert/hg19.gene_name.promoters.bed',stringsAsFactors = F)


# load expression data
ET <- read.delim('/Users/simeon01/Documents/PDTX/Expression_Data_OSCAR/NEW_UPDATED_OSCAR/ExpModelsData_all_plus_PARsamples.txt', stringsAsFactors = F, sep = "\t",header = T)

# ======  use this renaming if we go for new + old PDTX ===========
# ET and DG4 do not have the exact same name. Fix this.
colnames(ET)[which(colnames(ET)=='STG331.T1.P02..36366.T1.FF1')] <- '331_284'
colnames(ET)[which(colnames(ET)=='HCI009.M1.P03..35978.T5.FF1')] <- '9_317'
colnames(ET)[which(colnames(ET)=='VHIO179.M1.P04..35942.T1.FF1')] <- '179'
colnames(ET)[which(colnames(ET)=='VHIO098.M1.P02..39866.T1.FF1')] <- '98'
colnames(ET)[which(colnames(ET)=='STG316.T1.PPrimary..T1.FF1')] <- '316_284'
colnames(ET)[which(colnames(ET)=='HCI005.M1.P03..39710.T8.FF1')] <- '5_317'
colnames(ET)[which(colnames(ET)=='AB521.M1.P04..22113.T1.FF1')] <- '521_284'
colnames(ET)[which(colnames(ET)=='STG139M.M1.P07..91802.T1.FF1')] <- '139M_'
#colnames(ET)[which(colnames(ET)=='AB863.T2.P01..19920.T1.FF1')] <- '863_317'
colnames(ET)[which(colnames(ET)=='STG139.T1.X02.00.000000.T')] <- 'STG139'
colnames(ET)[which(colnames(ET)=='HCI005.M1.P03..39710.T8.FF1')] <- '5_317'
colnames(ET)[which(colnames(ET)=='STG201.T1.X00.00.000000.T')] <- '201'

colnames(ET)[which(colnames(ET)=='AB636.M1.P02..49124.T1.FF1')] <- 'AB636'
colnames(ET)[which(colnames(ET)=='STG282.M1.P04..12621.T1.FF1')] <- 'STG282'
colnames(ET)[which(colnames(ET)=='AB577.M1.P01..25994.T1.FF1')] <- 'AB577'
colnames(ET)[which(colnames(ET)=='AB551.M3.P02..34058.T2.FF1')] <- 'AB551'
colnames(ET)[which(colnames(ET)=='STG195.M1.P04..45218.T1.FF1')] <- 'STG195'
colnames(ET)[which(colnames(ET)=='AB863.T2.PPrimary..T1.FF1')] <- '863_317'
colnames(ET)[which(colnames(ET)=='STG143.T1.P03..29103.T1.FF1')] <- '143'
colnames(ET)[which(colnames(ET)=='AB580.T1.P01..23191.T1.FF1')] <- 'AB580'
colnames(ET)[which(colnames(ET)=='AB555.T1.x00.00.000000.T')] <- 'AB555'
colnames(ET)[which(colnames(ET)=='STG335.T1.P02..12259.T12.FF1')] <- 'STG335'


colnames(DG4_annotation_on_promoters) <- gsub('peak_counts.norm_filtered.detable.','',colnames(DG4_annotation_on_promoters))
colnames(DG4_annotation_on_promoters)[which(colnames(DG4_annotation_on_promoters)=='PAR1006_3')] <- 'PAR1006'
colnames(DG4_annotation_on_promoters)[which(colnames(DG4_annotation_on_promoters)=='179')] <- '179'
colnames(DG4_annotation_on_promoters)[which(colnames(DG4_annotation_on_promoters)=='179')] <- '179'
colnames(DG4_annotation_on_promoters)[which(colnames(DG4_annotation_on_promoters)=='139')] <- 'STG139'
colnames(DG4_annotation_on_promoters)[which(colnames(DG4_annotation_on_promoters)=='139M')] <- '139M_'

#colnames(ET)<- c( "X","Geneid","AB521","AB863.T2.PPrimary..T1.FF1","HCI005.M1.P03..39710.T8.FF1","HCI009.M1.P03..35978.T5.FF1","STG139M.M1.P07..91802.T1.FF1", "STG143.T1.P03..29103.T1.FF1","STG201.T1.X00.00.000000.T","STG316.T1.PPrimary..T1.FF1","VHIO098.M1.P02..39866.T1.FF1","VHIO179.M1.P04..35942.T1.FF1","AB790.T1.P01..45197.T1.FF1","AB863.T2.P01..19920.T1.FF1","AB577.M1.P01..25994.T1.FF1","AB580.T1.P01..23191.T1.FF1","AB582.T1.PPrimary..T1.FF1","AB636.M1.P02..49124.T1.FF1","AB551.M3.P02..34058.T2.FF1","AB555.T1.x00.00.000000.T","STG139.T1.X02.00.000000.T","STG282.M1.P04..12621.T1.FF1","STG195.M1.P04..45218.T1.FF1","STG335.T1.P02..12259.T12.FF1","STG331.T1.P02..36366.T1.FF1","PAR1006" ,"PAR1040","PAR1022" )

new_col_names_EXPRESSION <- c()
new_reordered_EXPRESSION <- c()
flag=1
for (i in (5:length(colnames(DG4_annotation_on_promoters)))) {
  
  curr_pdtx <- colnames(DG4_annotation_on_promoters)[i]
  print("--- ")
  print(i)
  print(curr_pdtx)
  #idx <- which(colnames(ET) == curr_pdtx)
  idx <- grep(curr_pdtx,colnames(ET))
  print(colnames(ET)[idx])
  print("=====")
  if (length(idx)==1) {
    print(' == > name found <== ')
    new_col_names_EXPRESSION <- c(new_col_names_EXPRESSION,curr_pdtx)
    new_reordered_EXPRESSION <- cbind(new_reordered_EXPRESSION,ET[,idx])
    flag = flag+1
  }else {print('not found ****************')}
  
}
# ======  END use this renaming if we go for new + old PDTX  END ===========
# ==========================================================================
colnames(new_reordered_EXPRESSION) <- new_col_names_EXPRESSION
ET <- cbind(ET[,2],new_reordered_EXPRESSION)
colnames(ET)[1] <- 'Gene'

# left join
colnames(DG4_annotation_on_promoters)[4] <- 'Gene'
ET_df <-as.data.frame(ET)
DG4_annotation_and_expression <- left_join(DG4_annotation_on_promoters,ET_df, by='Gene')

# for each cancer type create a dynamic structures that contains the fraction of genes with paak in
# 3 differeent category of exprssion( low-med-high)

#to identify PDTX model name check colnames

PDTX_list <- colnames(DG4_annotation_on_promoters)[5:length(colnames(DG4_annotation_on_promoters))]
print(PDTX_list)

#PDTX_list <- gsub('_dg4','',PDTX_list)

summary_freq <- data.frame(expression_subsetting = c(1,2,3))
flag=1;

for (i in 1:length(PDTX_list)) {
  curr_PDTX <- PDTX_list[i]
  print(curr_PDTX)
  subset_DG4_annotation_and_expression <- DG4_annotation_and_expression %>% dplyr::select(contains(curr_PDTX))
  dim(subset_DG4_annotation_and_expression)
  if (dim(subset_DG4_annotation_and_expression)[2] <2) 
  {
    next 
  }
  
  tot_N_peaks <- tot_N_peaks_all[grep(curr_PDTX,names(tot_N_peaks_all))]#tot_N_peaks_all[i]
  subset_DG4_annotation_and_expression <-  drop_na(subset_DG4_annotation_and_expression)
  dim(subset_DG4_annotation_and_expression)
  
  subset_DG4_annotation_and_expression[which(subset_DG4_annotation_and_expression[,1]>=1),1] <- 1
  
  colnames(subset_DG4_annotation_and_expression) <- c('G4_presence','expression')
  
  subset_DG4_annotation_and_expression$expression_subsetting <- as.numeric(cut2(as.numeric(subset_DG4_annotation_and_expression$expression), g=3))
  
  summary_cases <- subset_DG4_annotation_and_expression %>% group_by(expression_subsetting) %>% tally()
  subset_DG4_annotation_and_expression %>% group_by(expression_subsetting) %>% tally()
  
  # here you can choose what fraction to use
  ## fraction = N promoters with quadruplex over those annotated
  #A <-  subset_DG4_annotation_and_expression %>% group_by(expression_subsetting) %>% summarise(Frequency = sum(G4_presence)/sum(subset_DG4_annotation_and_expression[,1]))
  ## fraction = N promoters with quadruplex over total number of peaks
  A <-  subset_DG4_annotation_and_expression %>% dplyr::group_by(expression_subsetting) %>% dplyr::summarise(Frequency = sum(G4_presence)/tot_N_peaks)
  write.table(subset_DG4_annotation_and_expression,file=paste0('subset_DG4_annotation_and_expression_',PDTX_list[i],'.csv'),quote = F,sep= ",", col.names = NA)
  summary_freq <- left_join(summary_freq,A ,by='expression_subsetting')
  colnames(summary_freq)[flag+1] <- curr_PDTX
  
  flag <- flag+1
  print(flag)
}


df_summary_freq <- as.data.frame(t(summary_freq))
df_summary_freq <- df_summary_freq[-1,]
colnames(df_summary_freq) <- c('Low','Med','High')
df_summary_freq$pdtx <- rownames(df_summary_freq)
df_summary_freq$pdtx <- factor(df_summary_freq$pdtx, levels = df_summary_freq$pdtx[order(df_summary_freq$High)])
library(reshape)
df_summary_freq_melt <- melt(df_summary_freq, id.vars='pdtx')

save(df_summary_freq_melt,file='df_summary_freq_melt.Rdata')
save(df_summary_freq,file='df_summary_freq.Rdata')
write.table(df_summary_freq,file='df_summary_freq.csv', quote = F, sep= ",",col.names = NA)
write.table(df_summary_freq_melt,file='df_summary_freq_melt.csv', quote = F, sep= ",",col.names = NA)

ggplot(df_summary_freq_melt, aes(x=pdtx,y=value, group = variable)) + geom_bar(aes(fill = variable),stat="identity",position="dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('fraction of peaks in promoters') + ylab('PDTX') +
  ggtitle('fraction of DG4R positive promoters stratified by gene expression levels') +labs(x = "", fill = "Expression\nlevel group")
ggsave('Fraction_DG4R_positive_promoters_stratified_by_expression.pdf')

## ========================================
## check expression of selected helicases
## ========================================


# load Helicases info
helicases <- read.delim('helicase_list_angela.csv', sep=",",stringsAsFactors = F)
colnames(helicases)[1] <- "Gene"
# merge helicases with Expression data ET


Helicases_expression <- left_join(helicases,ET_df, by='Gene')
Helicases_expression$Gene_genefunction <- paste0(Helicases_expression$Gene," ++ ",Helicases_expression$G4_Function)
Helicases_expression <- select(Helicases_expression,-c(GENE_ID,uniprot,protein,G4_Function,Gene_genefunction))
Helicases_expression <- as_tibble(Helicases_expression)
Helicases_expression <- mutate_if(Helicases_expression, is.factor, ~ as.numeric(as.character(.x)))
Helicases_expression <- Helicases_expression %>% dplyr::group_by(Gene) %>% summarise_all(funs(mean(., na.rm=TRUE)))

Helicases_expression_subset <- data.matrix(Helicases_expression[,2:23])
rownames(Helicases_expression_subset) <- Helicases_expression$Gene
# for each helicase plot expression vs N peaks
vect1=tot_N_peaks_all
for (i in (1:length(rownames(Helicases_expression_subset)))){
  helicase_name <- rownames(Helicases_expression_subset)[i]
  pdf(paste0('./Helicases/Helicases_expression_vs_N_DG4r',helicase_name,'.pdf'))
  vect2=Helicases_expression_subset[i,]
  cor_exp_ndg4r <- round(cor(as.numeric(vect1),as.numeric(vect2)),2)
  par(mfrow=c(1,1))
  plot(vect1, vect2, main=paste0(helicase_name,'\nr:',cor_exp_ndg4r), xlab='N DG4r',ylab='gene expr.')
  dev.off()
  #readline(prompt="Press [enter] to continue")
}
  


for (i in (2:23)){
  #plot(1:22,Helicases_expression[i,6:(6+21)], main=Helicases_expression$Gene[i])
  pdf(paste0('./Helicases/Helicases_expression_',colnames(ET_df)[i],'.pdf'))
  par(mfrow=c(1,2))
  boxplot(as.numeric(levels(ET_df[,i])),horizontal = F, outline = T,names = colnames(ET_df)[i],
          show.names=TRUE, main=colnames(ET_df)[i],
          xlab='gene expr. (RPKM)'); 
  expression_hel <- Helicases_expression_subset[,i-1]
  index_order <- order(expression_hel,decreasing = F)
  abline(h=expression_hel, col = "red")
  barplot(expression_hel[index_order] ,names.arg =Helicases_expression$Gene[index_order], las=2, horiz = T,
          main=colnames(ET_df)[i],xlim=c(0,max(expression_hel)))
  #readline(prompt="Press [enter] to continue")
  dev.off()
}

A <- t(Helicases_expression[,c(1,6:dim(Helicases_expression)[2])])
A <- A[-1,]
Helicases_expression_sel_df <- melt.array(A)




colnames(Helicases_expression_sel_df)[2] <- "id_n"
mapping_gene_index <- data.frame(Gene=Helicases_expression$Gene,id_n =1:length(Helicases_expression$Gene))
Helicases_expression_sel_fields2_df <- left_join(Helicases_expression_sel_df,mapping_gene_index,by="id_n")


Helicases_expression_sel_fields2_df %>% ggplot(aes(x=X1,y=value)) + geom_line()
