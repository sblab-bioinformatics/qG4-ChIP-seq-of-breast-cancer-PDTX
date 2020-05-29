#' @title TF fold-enrichments at ΔG4r
#' @description explore TF fold-enrichments at ΔG4r
#' @param files downloaded from ChIP_atlas and pre-processed 
#' @return FR matrix and correlations matrices as heatmap and as txt
#' @author angela simeone  \email{angela.simeone@cruk.cam.ac.uk}

setwd('~/ChIP_atlas')
rm(list = ls())

# == process the txt filed keeping onyl those with logPvalue < -3, remove ing
# for file in *up.txt
# do
# awk -F "\t" '{if ($9 <- 3) print $1"\t"$3"\t"$9"\t"$11}' $file |grep -v inf | sed 's/SRX150480_E2F4/SRX150480_E2F4.SRC/g' > ${file%%.txt}.tab2
# done

## == create the list of unique TF symbols to use for the building up the matrix with all pairwise correlations
# rm temp_list_tf_robstyle.txt
# touch temp_list_tf_robstyle.txt
# for file in *up.txt; do  awk  -F "\t"  '{print $3}' $file >> temp_list_tf_robstyle.txt ; done
# sort temp_list_tf_robstyle.txt | uniq  > list_tf_robstyle.txt

## == import libraries
library(tidyverse)
library(corrplot)
library("Hmisc")
files=list.files(pattern='*tab2')

## == reads list of unique TF list (gene names)
temp1 <- read.table('list_tf_robstyle.txt',sep = "\t", header = F, stringsAsFactors = F)
head(temp1)

## == initialize the TF fold enrichments across various PDTX
MM_all <- temp1
colnames(MM_all)<- c('tf')
for (i in 1:length(files))  {
  print(i)
  print(files[i])
  # read file and select only tf and fold enrichment
  temp_read_full <- read.table(files[i],sep = "\t", header = F, stringsAsFactors = F)
  temp_read <- temp_read_full[,c(2,4)]
  colnames(temp_read)[1:2] <- c('tf','fe')
  #temp_read$id_tf <- paste(temp_read_full[,1],temp_read_full[,2],sep="_")
  print(dim(temp_read))
  head(temp_read)
  print(dim(temp_read))
  
  # group by tf and summarize fe with the median value
  temp_read_median <- temp_read %>% dplyr::group_by(tf) %>% dplyr::summarise(meandian_fe=median(fe, na.rm = T))
  
  # join the big summary table with the current summarized PDTX fe
  MM_all <- left_join(MM_all,temp_read_median,by="tf")
}

# drop NA
#MM_all_no_NA <- MM_all %>%  drop_na()
# update column names
colnames(MM_all)[2:dim(MM_all)[2]] <- gsub('_vs_all_up.tab2','',files)

rownames(MM_all) <- MM_all$tf
dim(MM_all)



## == remove mC and 8-Hydroxydeoxyguanosine
# MM_all <- MM_all[-c(grep('mC',MM_all[,1]),grep('8-Hydroxydeoxyguanosine',MM_all[,1]),grep('Epitope tags',MM_all[,1])),]
# dim(MM_all)

## == sostitute NA with zeros (0) --> those zeros results from the fact that some of the TF are absent from the rest
MM_all_noNA <- MM_all
MM_all_noNA[is.na(MM_all)] <- 0
#MM_all_numeric_temp <- MM_all %>%  mutate_all(~replace(., is.na(.), 0))
MM_all_numeric_temp <- data.matrix(MM_all_noNA[,2:dim(MM_all_noNA)[2]])
rownames(MM_all_numeric_temp) <- MM_all[,1]
dim(MM_all_numeric_temp)

## == loop to keep only the rows where minimum fe is ==20
MM_all_new <- c()
MM_all_new_rownames <- c()
for (i in c(1:dim(MM_all_numeric_temp)[1])){
  tot_zeros <- length(which(MM_all_numeric_temp[i,]==0))
  print(paste0("total zeros ",tot_zeros))
  min_row <- mean(as.numeric(MM_all_numeric_temp[i,]))
  print(i)
  print(min_row)
  if (tot_zeros < 22){
    print("in th loop")
     MM_all_new <- rbind(MM_all_new,MM_all_numeric_temp[i,])
     MM_all_new_rownames <- c(MM_all_new_rownames,rownames(MM_all_numeric_temp)[i])
   }
}
dim(MM_all_new)
rownames(MM_all_new) <- MM_all_new_rownames

write.table(MM_all_new, file='Fold_enrichments_DG4R_at_TFbindingRegions.processed.csv',quote = F, sep = ",", col.names = NA)
# MM_all_numeric_temp <- MM_all_new %>%  mutate_all(~replace(., is.na(.), 0))
# MM_all_numeric <- as.matrix(MM_all_numeric_temp[,2:dim(MM_all_numeric_temp)[2]])
# rownames(MM_all_numeric) <- rownames(MM_all_new)

## == read rob data

MM_rob <- read.table('matrix.RHH.processed_AS.TF.DG4R.enrichment_no_empty_lines.txt', sep = "\t",stringsAsFactors = F, header = T)


#write.table(MM_all_new,file = 'matrix_TF_fold_enrichments_cross_PDTX.csv',sep = ",",quote = F,col.names = NA)
## == compute spearman correlation between TF
#C <- cor(t(MM_all_numeric_temp[,2:23]))
C <- cor(t(data.matrix(MM_all_new)),method = "spearman")
C_rob <- rcorr(as.matrix(t(MM_all_new)), type = "spearman")
dim(C)
dim(C_rob$r)

## == plot and save spearman correlation between TF
#C_drop_NA <- as.data.frame(C) %>% drop_na()
#C_drop_NA <- as.data.frame(C) %>% purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=20)
col1 <- colorRampPalette(c( "red", "yellow", "blue"))

pdf('TB_chipAtlas_FE_corrplot_median.pdf')
CorrTF <- corrplot(C, type="lower", order = "hclust", hclust.method = "complete",addrect = 3, col = col1(1000),tl.cex=0.3, tl.col = "black", tl.pos = "ld", addgrid.col=NA)
#CorrTF <- corrplot(C, order = "hclust", hclust.method = "complete",addrect = 4, col = col1(1000),tl.cex=0.3, tl.col = "black", addgrid.col=NA)
dev.off()

pdf('TB_chipAtlas_FE_corrplot_median_ROB.pdf')
#CorrTF_ROB <- corrplot(C_rob$r, method = "circle",order = "hclust", hclust.method="ward.D2",addrect =2 ,col = col1(100),tl.cex=0.3)
CorrTF_ROB <- corrplot(C_rob$r, type="lower", order = "hclust", hclust.method = "complete", col = col1(1000), tl.col = "black", tl.pos = "ld",
         p.mat = C_rob$P, sig.level = 0.9, insig = "blank", addgrid.col=NA,tl.cex=0.3)
dev.off()

## == import annotations PDTX
tumor_annotations <- read.csv('PDTX_annotations_classification.txt',sep="\t",header = T, stringsAsFactors = F)

## == rename some of the pdtx
tumor_annotations$Tumour[grep('316',tumor_annotations$Tumour)] <- '316_284'
tumor_annotations$Tumour[grep('863',tumor_annotations$Tumour)] <- '863_317'
tumor_annotations$Tumour[grep('HCI009',tumor_annotations$Tumour)] <- '9_317'
tumor_annotations$Tumour[grep('HCI005',tumor_annotations$Tumour)] <- '5_317'
tumor_annotations$Tumour[grep('STG331',tumor_annotations$Tumour)] <- '331_284'
tumor_annotations$Tumour[grep('STG139M',tumor_annotations$Tumour)] <- '139M_'
tumor_annotations$Tumour[grep('AB521',tumor_annotations$Tumour)] <- '521_284'
tumor_annotations$Tumour[grep('PAR1006',tumor_annotations$Tumour)] <- "PAR1006_3"
tumor_annotations$Tumour[grep('STG201',tumor_annotations$Tumour)] <- "201_"
tumor_annotations$Tumour[grep('STG143',tumor_annotations$Tumour)] <- "143_"
tumor_annotations$Tumour[grep('VHIO179',tumor_annotations$Tumour)] <- "179_"
tumor_annotations$Tumour[grep('VHIO098',tumor_annotations$Tumour)] <- "98_"

## == compute spearman correlation between PDTX
C_pdtx <- cor(data.matrix(MM_all_new),method = 'spearman')^2

## == rename correlation between PDTX 
new_PDTX_names_annotations <- colnames(C_pdtx)
#loop over col and rownames to include tumor annotations
for (i in 1:length(colnames(C_pdtx))){
  print(colnames(C_pdtx)[i])
  A <- grep(colnames(C_pdtx)[i],tumor_annotations$Tumour)
  print(A)
  new_PDTX_names_annotations[i] <- paste0(colnames(C_pdtx)[i]," == ",tumor_annotations$annotations[A])
  if (is_empty(A)){
    print(paste0(colnames(C_pdtx)[i],"not found"))
}
}
colnames(C_pdtx) <- rownames(C_pdtx) <- new_PDTX_names_annotations
dim(C_pdtx)

## == plot and save spearman correlation between TF
col2 <- colorRampPalette(c("white", "orange", "yellow", "blue"))
pdf('TB_chipAtlas_corrplot_median_PDTX.pdf')
corrplot(C_pdtx, method = "circle",order = "hclust", hclust.method="ward.D",col = col2(50),tl.cex=0.8,addrect = 5,cl.lim = c(0,1))
dev.off()


## == export correlations to csv files == 
write.table(C, file='Corr_Spearm_TF_Fold_enrichments_DG4R_at_TFbindingRegions.processed.csv',quote = F, sep = ",", col.names = NA)
write.table(C, file='Corr_Spearm_PDTX_Fold_enrichments_DG4R_at_TFbindingRegions.processed.csv',quote = F, sep = ",", col.names = NA)

# ## == % for this section, I manually selected the TF that are in each of the rectangles 
# ##    % because I cannot export the outcome of corrplot after stating addrect
# 
# ## ==  individual Fold Enrichments programs (2 programs)
# Program1 <- data.frame(tf=c("TP53","FOSL2","TP63","CEBPB","FOS","FOSL1","TLE3","PGR","EP300","FOXA1","AHR","ARID1A",'GATA3',"ESRRA","ESR1","TFAP2A","NR2F2","NR3C1","KMT2C","YAP1"))
# Program2 <- data.frame(tf=setdiff(rownames(MM_all_new),Program1$tf))
# 
# ## ==  Extract fold-enrichment set1 and plot
# Program1_dt <- left_join(Program1,MM_all,by='tf')
# Program1_M <- data.matrix(Program1_dt[,-1])
# rownames(Program1_M) <- Program1_dt$tf
# colnames(Program1_M) <- new_PDTX_names_annotations
# Program1_M[is.na(Program1_M)] <- 0
# col3=colorRampPalette(c("blue","yellow","red"))
# pdf('set1.pdf')
# corrplot(Program1_M,is.corr=F,col = col3(50),method = "color", tl.cex=0.3)
# dev.off()
# 
# ## ==  Extract fold-enrichment set2 and plot
# Program2_dt <- left_join(Program2,MM_all,by='tf')
# Program2_M <- data.matrix(Program2_dt[,-1])
# colnames(Program2_M) <- new_PDTX_names_annotations
# rownames(Program2_M) <- Program2_dt$tf
# Program2_M[is.na(Program2_M)] <- 0
# col3=colorRampPalette(c("blue","yellow","red"))
# pdf('set2.pdf')
# corrplot(Program2_M,is.corr=F,col = col3(50),method = "color", tl.cex=0.3)
# dev.off()
# 
# ## ==  concatenate set1 and set2
# PM <- rbind(Program1_M,Program2_M)
# 
# # to extract names in the 2 cluster use dimnames(CorrTF)??
# 
# ## ==  visualize set1 and set2 with complexHeatmap function Heatmap
# library(ComplexHeatmap)
# library(circlize)
# 
# f1 = colorRamp2(seq(min(PM), max(PM), length = 3), c("cornflowerblue", "yellow", "red"))
# ht1 = Heatmap(Program1_M,col = f1, cluster_rows=F,cluster_columns=F,row_names_gp = gpar(fontsize = 5),column_names_gp=gpar(fontsize = 5),name = "set1")
# ht2 = Heatmap(Program2_M,col = f1, cluster_rows=F,cluster_columns=F,row_names_gp = gpar(fontsize = 5),column_names_gp=gpar(fontsize = 5),name = "set2")
# 
# ht1 
# ht2
# 
# pdf('complex_heatmap_set1.pdf')
# draw(ht1)
# dev.off()
# 
# pdf('complex_heatmap_set2.pdf')
# draw(ht2)
# dev.off()


## attempts to boxplot?!
# group each of the 2groups by classes 
# boxplot(as.vector(Program1_M),as.vector(Program2_M))

               
# boxplot(as.vector(Program1_M),as.vector(Program2_M), outline = F)
# boxplot(as.vector(Program1_M[,-which(grepl("IC10",colnames(Program1_M)))]),as.vector(Program1_M[,which(grepl("IC10",colnames(Program1_M)))]), outline = T)
# boxplot(as.vector(Program2_M[,-which(grepl("IC10",colnames(Program2_M)))]),as.vector(Program2_M[,which(grepl("IC10",colnames(Program2_M)))]), outline = T)
# boxplot(as.vector(Program1_M[,which(grepl("ER-",colnames(Program1_M)))]),as.vector(Program1_M[,which(grepl("ER+",colnames(Program1_M)))]), outline = F)
# boxplot(as.vector(Program2_M[,which(grepl("ER-",colnames(Program2_M)))]),as.vector(Program2_M[,which(grepl("ER+",colnames(Program2_M)))]), outline = F)
# 

# import ER
ET <- read.delim('/Users/simeon01/Documents/PDTX/Expression_Data_OSCAR/NEW_UPDATED_OSCAR/ExpModelsData_all_plus_PARsamples.txt', stringsAsFactors = F, sep = "\t",header = T)

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
head(ET)
ET2 <- ET[,-1]
colnames(ET2) <- c("Geneid",paste0('Expr_',colnames(ET2)[2:dim(ET2)[2]]))
ET_df <- data.frame(ET2)

ET_df_nodup <- ET_df %>% dplyr::group_by(Geneid) %>% summarise_all(funs(mean(., na.rm=TRUE)))
colnames(ET_df_nodup) <- colnames(ET_df)

MM_all_new_temp <- MM_all_new
colnames(MM_all_new_temp) <- paste0('FE_',colnames(MM_all_new_temp))
MM_all_new_df <- data.frame(MM_all_new_temp)
MM_all_new_df$Geneid <- rownames(MM_all_new_temp)
  
MM_FE_EXP <- left_join(MM_all_new_df,ET_df_nodup, by="Geneid")

dim(MM_FE_EXP)

# split again in 2 matrices
MM_FE <- data.matrix(MM_FE_EXP[,1:22])
MM_Exp <- data.matrix(MM_FE_EXP[,24:49])
vect_sorted <- c()
for (i in 1:22) {
  
  temp<- colnames(MM_all_new)[i]
  
  if(length(grep(temp,colnames(MM_Exp)))>0)
  {
    vect_sorted <- c(vect_sorted,grep(temp,colnames(MM_Exp)))
    print(grep(temp,colnames(MM_Exp), value = T))
    } else {print(paste0('not_found: ',temp))}
}

length(vect_sorted)

MM_Exp_updated <- MM_Exp[,vect_sorted]
C_fe_expr_temp <- c()
for (i in 1:134) {
  vect1=MM_FE[i,]
  vect2=MM_Exp_updated[i,]
  vect1[is.na(vect1)] <- 0
  vect2[is.na(vect2)] <- 0
  cor_fe_expr <- round(cor(as.numeric(vect1),as.numeric(vect2),method = "spearman"),2)
  C_fe_expr_temp <- c(C_fe_expr_temp,cor_fe_expr)
  # pdf(paste0('scatter_plot_FE_expr',MM_FE_EXP$Geneid[i],'.pdf'))
  # plot(vect1, vect2, main=paste0('Symb:',MM_FE_EXP$Geneid[i],'\n',cor_fe_expr), xlab='FE',ylab='log2(Expr)')
  # dev.off()
}
names(C_fe_expr_temp) <- MM_FE_EXP$Geneid
C_fe_expr <- data.frame(Cspeaman=C_fe_expr_temp,Geneid=MM_FE_EXP$Geneid)

head(cbind(MM_FE,MM_Exp_updated,C_fe_expr))

# write.table(cbind(MM_FE,MM_Exp_updated,C_fe_expr),file='TF_FE_and_Expression_and_Cspearman.csv',quote = F,sep=",",col.names = NA)
# pdf('barplot_Cspearman_FE_expr.pdf')
# barplot(C_fe_expr$Cspeaman, ylab="corelation",xlab='PDTXs')
# dev.off()

# import Robert annotations
Rob_anno <- read.delim('/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/Figure.3ROB_input/PDTX_annotations.txt' ,header = F, stringsAsFactors = F,sep = "\t")

R_anno_df <- as.data.frame(t(Rob_anno[,-1]))
head(R_anno_df)
colnames(R_anno_df) <- Rob_anno[,1]

#transpose the matrices with FE and expression in order to have TF on the columns and PDTX id on the rows and then join annotations

#matrix with FE
MM_Exp_updated_t <- t(data.matrix(MM_Exp_updated))
colnames(MM_Exp_updated_t) <- rownames(C_fe_expr)
MM_Exp_updated_t_df <- as.data.frame(MM_Exp_updated_t)
MM_Exp_updated_t_df$PDTX <- rownames(MM_Exp_updated_t_df)
MM_Exp_updated_t_df_anno <- left_join(MM_Exp_updated_t_df,R_anno_df,by="PDTX")
MM_Exp_updated_t_df_anno <- MM_Exp_updated_t_df_anno[ -c(18,19),]

#matrix with expression
MM_FE_t <- t(data.matrix(MM_FE))
colnames(MM_FE_t) <- rownames(C_fe_expr)
MM_FE_t_df <- as.data.frame(MM_FE_t)
MM_FE_t_df$PDTX <- rownames(MM_FE_t)
MM_FE_t_df_anno <- left_join(MM_FE_t_df,R_anno_df,by="PDTX")


# now import the various clusters (there are 7)
list_cluster_files <- list.files(path = "/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/Figure.3ROB_input/", pattern = "CLUSTER_.*.listTF.txt")
setwd("/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/Figure.3ROB_input/")

cluster_TF_bound <- lapply(list_cluster_files, read.delim, sep=',')

## for each cluster do the following plots
## ER+vs ER- for FE
## IC8/1 vs IC10/9 for FE
##
## ER+vs ER- for expr
## IC8/1 vs IC10/9 for expr

library(gridExtra)
library(grid)


for (i in (1:7)){
  file_tag <- gsub('_sorted.matrix_clust_TF_spearman_wardd2_AS_listTF.txt','',list_cluster_files[i])
  tmp <- unlist(levels(cluster_TF_bound[[i]]$PDTX))
  print(file_tag)
  print(length(tmp))
  #FE -- ER
  tmp1_FE <- MM_FE_t_df_anno %>% select( one_of(c(tmp,'ER-status'))) %>% mutate(ER= ifelse(`ER-status`=='neg','b_ER-','a_ER+')) %>% select(-c(`ER-status`))
  tmp2_FE <- melt(tmp1_FE, id="ER")
  write.table(tmp2_FE,file=paste0(file_tag,"_FE_ER.csv"),quote = F, sep= ",", col.names = NA)
  p_FE <- ggplot(tmp2_FE,aes(y=value, x=ER)) + geom_boxplot()
  p1 <- p_FE + stat_compare_means(method = "wilcox.test") + ggtitle('FE')

  #FE -- IC
  tmp1_FE_IC <-  MM_FE_t_df_anno %>% select( one_of(c(tmp,'IC'))) %>%
    mutate(IC_type= ifelse(IC=="IC8" |  IC=="IC1",'a_IC8_1',
                           ifelse(IC=="IC10" | IC=="IC9",'b_IC10_9','NA'))) %>% select(-c("IC"))
  tmp2_FE_IC <- melt(tmp1_FE_IC, id="IC_type")
  write.table(tmp2_FE_IC,file=paste0(file_tag,"_FE_IC.csv"),quote = F, sep= ",", col.names = NA)
  p_FE_IC <- tmp2_FE_IC %>% filter (IC_type!="NA") %>% ggplot(aes(y=value, x=IC_type)) + geom_boxplot()
  p2 <- p_FE_IC + stat_compare_means(method = "wilcox.test") + ggtitle('FE')


  #EXP -- ER
  tmp1_Exp <- MM_Exp_updated_t_df_anno %>% select( one_of(c(tmp,'ER-status'))) %>% mutate(ER= ifelse(`ER-status`=='neg','b_ER-','a_ER+')) %>% select(-c(`ER-status`))
  tmp2_Exp <- melt(tmp1_Exp, id="ER")
  write.table(tmp2_Exp,file=paste0(file_tag,"_expression_ER.csv"),quote = F, sep= ",", col.names = NA)
  p_Exp <- ggplot(tmp2_Exp,aes(y=value, x=ER)) + geom_boxplot()
  p3 <- p_Exp + stat_compare_means(method = "wilcox.test") + ggtitle('Expr')

  #EXP -- IC
  tmp1_Exp_IC <-  MM_Exp_updated_t_df_anno %>% select( one_of(c(tmp,'IC'))) %>%
    mutate(IC_type= ifelse(IC=="IC8" | IC=="IC1",'a_IC8_1',
                           ifelse(IC=="IC10" | IC=="IC9",'b_IC10_9','NA'))) %>% select(-c("IC"))
  tmp2_Exp_IC <- melt(tmp1_Exp_IC, id="IC_type")
  write.table(tmp2_Exp_IC,file=paste0(file_tag,"_expression_IC.csv"),quote = F, sep= ",", col.names = NA)
  p_Exp_IC <- tmp2_Exp_IC %>% filter (IC_type!="NA") %>% ggplot(aes(y=value, x=IC_type)) + geom_boxplot()
  p4 <-p_Exp_IC + stat_compare_means(method = "wilcox.test") + ggtitle('Expr')

  p_final <- grid.arrange(p1, p2, p3, p4, ncol=2)
  ggsave(paste0(file_tag,"_boxplot.jpg"), p_final)
}
  
                    