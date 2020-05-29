#' @title compare G4 signal expression and presence of  ΔG4r in promoters
#' @description compare G4 signal and expression levels of genes with ΔG4r in promoters
#' @param promoters list, expression levels, list of promoters overlapping at least one ΔG4r and G4 signal (from differential analysis)
#' @return heatmap and txt files with frequency of genes with ΔG4r stratified by expression level (high-med-low) and with color codifying G4 levels
#' @author angela simeone  \email{angela.simeone@cruk.cam.ac.uk}



library(tidyverse)
library(Hmisc)
library(data.table)

## ## ## ## ## ## ## ## ## ## ## ##  ## 
##  load promoter information        ## 
## and relative G4 levels as log(cpm)## 
## ## ## ## ## ## ## ## ## ## ## ##  ## 

# promoter files- read files with peaks overlapping promoters having also promoters coordinates!!!! 

path_peaks <- '/PDTX_dm6_inputSubtraction_/promoter_lfc/' # this folder contains bed files that contain the normalized G4 signal at  ΔG4r produced at the differential analysis step
setwd(path_peaks)
Promoters_annotations <- read.delim("hg19.gene_name.promoters.bed",stringsAsFactors = F, header = F)
Promoters_annotations$promoter_id <- paste(Promoters_annotations$V4,Promoters_annotations$V1,Promoters_annotations$V2,Promoters_annotations$V3,sep="_")
colnames(Promoters_annotations)[4] <- 'Gene'
list_peaks <-list.files(path=path_peaks,pattern = 'lcpm.promoter.bed')

for (i in (1:length(list_peaks))){
  temp <- read.table(list_peaks[1],sep = "\t",header = F,stringsAsFactors = F)
  curr_pdtx_temp <- gsub('peak_counts.norm_filtered.detable.','',list_peaks[i])
  curr_pdtx <- gsub('_vs_all_up.lcpm.promoter.bed','',curr_pdtx_temp)
  curr_pdtx_label <- paste0("G4_cpm.",curr_pdtx)
  
  summarized_promoter_G4_temp<- data.frame(promoter_id=paste(temp$V4,temp$V1,temp$V2,temp$V3,sep="_"),
                                           G4_cpm=exp(temp$V8))
  
  #summarise by median of G4 lcpm
  summarized_promoter_G4 <- summarized_promoter_G4_temp %>% dplyr::group_by(promoter_id) %>% dplyr::summarise(G4_median_cmp=median(G4_cpm))
  
  #combine with promoter matrix
  Promoters_annotations <- left_join(Promoters_annotations,summarized_promoter_G4,by='promoter_id')
  colnames(Promoters_annotations)[grep('G4_median_cmp', colnames(Promoters_annotations))] <- curr_pdtx_label
}
Promoters_annotations_noNA <-  Promoters_annotations %>% drop_na()

## ## ## ## ## ## ## ## ## ##  ## 
#      load expression data    ##
## ## ## ## ## ## ## ## ## ##  ## 

ET <- read.delim('ExpModelsData_all_plus_PARsamples.txt', stringsAsFactors = F, sep = "\t",header = T)

# ======  use this renaming if we go for new + old PDTX ===========
# ET and DG4 do not have the exact same name. Fix this.
colnames(ET)[which(colnames(ET)=='STG331.T1.P02..36366.T1.FF1')] <- '331_284'
colnames(ET)[which(colnames(ET)=='HCI009.M1.P03..35978.T5.FF1')] <- '9_317'
colnames(ET)[which(colnames(ET)=='VHIO179.M1.P04..35942.T1.FF1')] <- '179_'
colnames(ET)[which(colnames(ET)=='VHIO098.M1.P02..39866.T1.FF1')] <- '98_'
colnames(ET)[which(colnames(ET)=='STG316.T1.PPrimary..T1.FF1')] <- '316_284'
colnames(ET)[which(colnames(ET)=='HCI005.M1.P03..39710.T8.FF1')] <- '5_317'
colnames(ET)[which(colnames(ET)=='AB521.M1.P04..22113.T1.FF1')] <- '521_284'
colnames(ET)[which(colnames(ET)=='STG139M.M1.P07..91802.T1.FF1')] <- '139M_'
#colnames(ET)[which(colnames(ET)=='AB863.T2.P01..19920.T1.FF1')] <- '863_317'
colnames(ET)[which(colnames(ET)=='STG139.T1.X02.00.000000.T')] <- 'STG139'
colnames(ET)[which(colnames(ET)=='HCI005.M1.P03..39710.T8.FF1')] <- '5_317'
colnames(ET)[which(colnames(ET)=='STG201.T1.X00.00.000000.T')] <- '201_'
colnames(ET)[which(colnames(ET)=='AB636.M1.P02..49124.T1.FF1')] <- 'AB636'
colnames(ET)[which(colnames(ET)=='STG282.M1.P04..12621.T1.FF1')] <- 'STG282'
colnames(ET)[which(colnames(ET)=='AB577.M1.P01..25994.T1.FF1')] <- 'AB577'
colnames(ET)[which(colnames(ET)=='AB551.M3.P02..34058.T2.FF1')] <- 'AB551'
colnames(ET)[which(colnames(ET)=='STG195.M1.P04..45218.T1.FF1')] <- 'STG195'
colnames(ET)[which(colnames(ET)=='AB863.T2.PPrimary..T1.FF1')] <- '863_317'
colnames(ET)[which(colnames(ET)=='STG143.T1.P03..29103.T1.FF1')] <- '143_'
colnames(ET)[which(colnames(ET)=='AB580.T1.P01..23191.T1.FF1')] <- 'AB580'
colnames(ET)[which(colnames(ET)=='AB555.T1.x00.00.000000.T')] <- 'AB555'
colnames(ET)[which(colnames(ET)=='STG335.T1.P02..12259.T12.FF1')] <- 'STG335'
colnames(ET)[which(colnames(ET)=='PAR1006')] <- 'PAR1006_3'

new_col_names_EXPRESSION <- c()
new_reordered_EXPRESSION <- c()
flag=1
for (i in (1:22)) {
  
  curr_pdtx <- gsub('G4_cpm.','',grep('G4_cpm.',colnames(Promoters_annotations), value = T))[i]
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

ET_df <-as.data.frame(ET)
# loop over the expression matrix and create the 3 groups
ET_df_updated <- data.frame(Gene=ET_df$Gene)

for (i in (2:23)){
  curr_PDTX <- colnames(ET_df)[i]
  subset_ET <- ET_df %>% dplyr::select("Gene",contains(curr_PDTX))
  colnames(subset_ET)[2] <- "expression"
  subset_ET_df <- data.frame(Gene=subset_ET$Gene,expression=as.numeric(as.vector(subset_ET$expression)))
  #subset_ET <- as.tibble(subset_ET)
  subset_ET_df2 <- subset_ET_df %>% dplyr::group_by(Gene) %>% summarise(medianExpression= median(expression, na.rm = T))
  
  subset_ET_df2$expression_subsetting <- as.numeric(cut2(as.numeric(subset_ET_df2$medianExpression), g=3))
  colnames(subset_ET_df2)[2] <- paste0('medianExpression.',curr_PDTX)
  colnames(subset_ET_df2)[3] <- paste0('expressionGroup.',curr_PDTX)
  
  ET_df_updated <- dplyr::left_join(ET_df_updated,subset_ET_df2, by="Gene")
}

G4_levels_and_expression <- left_join(Promoters_annotations_noNA,ET_df_updated, by='Gene')

## ## ## ## ## ## ## ## ## ##  ## 
## Promoter at CNA import data ##
## ## ## ## ## ## ## ## ## ##  ## 
#to identify PDTX model name check colnames
setwd('~/PDTX_dm6_inputSubtraction_/Promoter_based_CNA_expression_G4level/promoter_CNA')

cna_promoter_annotations <- list.files(path = ".",pattern=".100kb.sorted.txt")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ##
## loop over PDTX that interactively          ##
## produce the overlap between                ##
## the 3 parts!!                              ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ##

PDTX_summary_pdtx <- data.frame(exp_cna_group=c("1_NEUT","2_NEUT","3_NEUT",
                                                "1_AMP","2_AMP","3_AMP",
                                                "1_GAIN","2_GAIN","3_GAIN",
                                                "1_HETD","2_HETD","3_HETD",
                                                "1_HOMD","2_HOMD","3_HOMD"))
for (i in (1:22)) {
 
  temp_cna_promoter_annotations <- read.table(cna_promoter_annotations[i],stringsAsFactors = F)
  temp_cna_promoter_annotations$promoter_id <- paste(temp_cna_promoter_annotations$V4,
                                                     temp_cna_promoter_annotations$V1,
                                                     temp_cna_promoter_annotations$V2,
                                                     temp_cna_promoter_annotations$V3,sep="_")
  
  colnames(temp_cna_promoter_annotations)[4]<- "Gene"
  colnames(temp_cna_promoter_annotations)[5]<- "cna_type"
  head(temp_cna_promoter_annotations)
  
  curr_pdtx <- gsub('.promoters.allCNA.100kb.sorted.txt','',cna_promoter_annotations[i])
  temp_selected_G4_levels_and_expression <- G4_levels_and_expression %>% select(c("Gene","promoter_id",contains(curr_pdtx)))
  
  temp_CNA_G4_levels_and_expression <- inner_join(temp_selected_G4_levels_and_expression,temp_cna_promoter_annotations,by = "Gene")
  colnames(temp_CNA_G4_levels_and_expression)[5] <- 'expressionGroup'
  colnames(temp_CNA_G4_levels_and_expression)[3] <- 'G4_cpm'
  head(temp_CNA_G4_levels_and_expression)
  #plot
  summary_pdtx <- temp_CNA_G4_levels_and_expression %>% group_by(expressionGroup,cna_type) %>% summarise(integralG4=median(G4_cpm)) %>% drop_na() 
  summary_pdtx$pdtx <- paste(curr_pdtx,summary_pdtx$expressionGroup,summary_pdtx$cna_type,sep = "_")
  summary_pdtx$exp_cna_group <- paste(summary_pdtx$expressionGroup,summary_pdtx$cna_type,sep = "_")
  summary_pdtx %>% ggplot(aes(x= pdtx,y=1, fill = integralG4)) +geom_tile() + coord_flip()
  write.table(summary_pdtx,file=paste0(curr_pdtx,'summary_G4median_by_CNA_expressiongroups.txt'),quote = F, sep = "\t",col.names = NA)
  summary_pdtx %>% ggplot(aes(x= pdtx,y=1, fill = integralG4)) +geom_tile() + ylab("medianG4") + coord_flip() 
  ggsave(paste0(curr_pdtx,"summary_G4median_by_CNA_expressiongroups.pdf"),width = 6, height = 17)
  label=paste0("integralG4.",curr_pdtx)
  # other approach is to combine into a single matrix
  PDTX_summary_pdtx <- left_join(PDTX_summary_pdtx,summary_pdtx,by="exp_cna_group") %>% select(exp_cna_group,starts_with("integralG4")) %>% rename(!!label:=integralG4)

}

write.table(PDTX_summary_pdtx,file='All_PTDX_summary_G4median_by_CNA_expressiongroups.txt',quote = F, sep = "\t",col.names = NA)

PDTX_summary_pdtx <- PDTX_summary_pdtx %>% arrange(exp_cna_group)
PDTX_summary_pdtx_M <- data.matrix(PDTX_summary_pdtx[,2:23])
rownames(PDTX_summary_pdtx_M) <- PDTX_summary_pdtx$exp_cna_group
colnames(PDTX_summary_pdtx_M) <- colnames(PDTX_summary_pdtx_M)[2:23]
f1 = colorRamp2(seq(0, max(PDTX_summary_pdtx_M, na.rm = T), length = 3), c("yellow", "orange", "red"))
pdf('All_PTDX_summary_G4median_by_CNA_expressiongroups.pdf')
ht1 = Heatmap(PDTX_summary_pdtx_M, 
              col=f1,
              na_col="lightyellow",
              cluster_rows=F,cluster_columns=F,
              row_names_gp = gpar(fontsize = 7),
              column_names_gp=gpar(fontsize = 6),
              row_split=c(rep('C_lowExp',5),rep('B_medExp',5),rep('A_highExp',5)),
              name = "Promoter G4-signal \n (median cpm)") 
draw(ht1)
dev.off()
