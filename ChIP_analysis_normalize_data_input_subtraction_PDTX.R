#' @title drosophila normalization after input subtraction
#' @description Generate drosophila normalized human signal after input subtraction and drosophila normalization
#' @param human peak counts, sample file, normalization factor file
#' @return sample files and normalzation factors for drosophila normalization step
#' @author angela simeone  \email{angela.simeone@cruk.cam.ac.uk}

base_dir <- '/PDTX_dir'

setwd(base_dir)

output_dir_name <- 'PDTX_dm6_inputSubtraction_/'

library(gplots)
library(edgeR)
library(dendextend)
library(plyr)

# ====================
# == INPUT FILES =====
# ====================
# tab-separated file file containing peak file and read counts in each peak, for each samples: see peak_example.bedgraph
peak_file_counts_name <- 'temp_concat_hg19_libraries.txt'

# tab-separated file containing sample information see sample_info_example.txt
sample_file_name <- '/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/sample_file_PDTX_only.txt'

# tab-separated file containing the normalization factor for each sample (only one row): norm_example.txt
norm_file <- '/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_old_new_dm6_norm5_peak_recovery.txt'
#norm_file <- '/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_old_new_no_dm6.txt'
#norm_file <- 'PDTX_old_new_dm6_total_recovery.txt'
#norm_file <- 'PDTX_old_new_no_dm6.txt'



# ==============================
# ===== CONTROL PARAMETERS =====
# ==============================
log_file_name <- paste(gsub('\\/','',output_dir_name),'.log',sep='')
do_subtract_input <- 0
do_evalute_after_subtract_input <- 0
do_cnv <- 0
prefix <- output_dir_name
save_plots <- 1
save_tables <- 1
save_norm_summary_files <- 1
system(paste('mkdir ',prefix,sep=''))

# ====== Load input files ======
peaks <- read.table(peak_file_counts_name,stringsAsFactors = F,header=F,sep='\t')
sample_file <- read.table(sample_file_name,stringsAsFactors = F,header=F,sep='\t')
norm_factor <- as.numeric(read.table(norm_file,stringsAsFactors = F,header=F))

samples <- as.character(sample_file[1,])
group <- as.character(sample_file[2,])
libSize <- as.numeric(sample_file[3,])
y <- peaks[,-c(1,2,3)]
coords <- peaks[,c(1,2,3)]
names(y) <-samples
row.names(y) <- paste(peaks[,1],peaks[,2],peaks[,3],sep='_')
y_orig <- y
row_inds <- rep(-1,0)

sink(log_file_name)
print('subtract input log:')
if (do_subtract_input && dim(sample_file)[1] >= 5)
{
  tmp1 <- peaks[,c(4:(4+dim(sample_file)[2]-1))]
  colnames(tmp1) <- samples
  # do rpm for all
  for (i in 1:dim(tmp1)[2]) {
    tmp1[,i] <- tmp1[,i] / (libSize[i] / 1000000)
  }
  # do respective input subtraction for all (in the RPM domain)
  tmp2 <- tmp1
  for (i in 1:dim(tmp2)[2]) {
    col_ind <- grep(sample_file[5,i],colnames(tmp1))
    #tmp2[,i] <- tmp1[,i] - tmp1[,col_ind]
    tmp2[,i] <- pmax(0,tmp1[,i] - tmp1[,col_ind])
  }  
  # do reverse rpm for all and go back to original counts
  for (i in 1:dim(tmp2)[2]) {
    y[,i] <- tmp2[,i] * (libSize[i] / 1000000)
  }
} else {
  print("Input subtraction not performed!")
}

# filter out selected samples and inputs
inds <- rep(-1,0)
if (dim(sample_file)[1] >= 4)
  inds <- c(inds,which(sample_file[4,] == 1))
inds <- c(inds,grep('input',samples))

print(paste('files removed: ',paste(samples[inds],collapse=', '),sep=''))
sink()


if (length(row_inds) > 0)
{
  y <- y[-row_inds,]
  y_orig <- y_orig[-row_inds,]
  coords <- coords[-row_inds,]
  if (do_cnv && dim(sample_file)[1] >= 4)
    cnv_orig <- cnv_orig[-row_inds,]
}
if (length(inds) > 0)
{
  y <- y[,-inds]
  y_orig <- y_orig[,-inds]
}

samples <- samples[-inds]
group <- group[-inds]
norm_factor <- norm_factor[-inds]
libSize <- libSize[-inds]

rep_names <- samples
for (i in 1:length(group))
{
  suf <- 1+(i+3)%%4
  rep_names[i] <- paste(samples[i],".",suf,sep="")
}
colnames(y) <- rep_names

# y2 = raw data with cpm (will be used for distance matrix analysis)
y2 <- y_orig
for (i in 1:dim(y2)[2])
{
  y2[,i] <- 1000000 * y_orig[,i] / libSize[i]
}
dist_euc_no_cpm <- as.matrix(dist(t(y2)))

# y3 = raw data with cpm and normalization (can be cnv+dm6 or dm6 only) (will be used for distance matrix analysis)
y3 <- y
for (i in 1:dim(y3)[2])
{
  y3[,i] <- 1000000 * y[,i] / libSize[i]
  y3[,i] <- y3[,i] * norm_factor[i]
}
dist_euc_tot_cpm <- as.matrix(dist(t(y3)))

# y1 = raw data no cpm with normalization (can be cnv+dm6 or dm6 only) (will be used for dba analysis)
y1 <- y
for (i in 1:dim(y1)[2])
{
  y1[,i] <- y[,i] * norm_factor[i]
}

# y4 = row data - obtained after input subtraction - with cpm 
y4 <- y
for (i in 1:dim(y2)[2])
{
  y4[,i] <- 1000000 * y[,i] / libSize[i]
}

dist_euc_cpm_input_subtract <- as.matrix(dist(t(y4)))


if (save_plots) {
  # hierarchical clustering
  #pdf(paste(prefix,'hclust.norm.pdf',sep=''),12,12)
  par(mfrow=c(2,1),cex=0.7,font=1)
  hc <- hclust(dist(t(y2)))
  dend <- as.dendrogram(hc)
  pal <- c(1:length(unique(group)))
  
  pal_seq_mapping <- data.frame(sample_type=unique(group), color=pal)
  group_df <- data.frame(sample_type=group)
  colors_mapped <- dplyr::left_join(group_df,pal_seq_mapping,by="sample_type")
  # hierarchical clustering
  pdf(paste(prefix,'hclust.norm.pdf',sep=''),12,12)
  
  colorCodes <- colors_mapped$color
  labels_colors(dend) <- colorCodes[order.dendrogram(dend)]
  
  plot(dend,main='no normalization, cpm')
  hc <- hclust(dist(t(y3)))
  dend <- as.dendrogram(hc)
  labels_colors(dend) <- colorCodes[order.dendrogram(dend)]
  plot(dend,main='dm6 normalization, cpm')
  dev.off()
  
  # PCA similarity analysis  
  pdf(paste(prefix,'pca.norm.pdf',sep=''),12,6)
  par(mfrow=c(1,2))
  pcaResult<-prcomp(t(y2))
  plot(pcaResult$x,
       main= 'Principal components, no normalization, cpm',
       xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
       ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
       type= 'n'
  )
  #text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= gsub('\\..*','',rownames(pcaResult$x)), cex= 0.7, col= c(rep('firebrick1',8),rep('green',4), rep('orange',8),rep('darkgreen',8),rep('blue',8), rep('darkgrey',4),rep('black',4), rep('navy',8), rep('purple',4)))
  text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= gsub('\\..*','',rownames(pcaResult$x)), cex= 0.7,col = colorCodes)
  pcaResult<-prcomp(t(y3))
  plot(pcaResult$x,
       main= 'Principal components, dm6 normalization, cpm',
       xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
       ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
       type= 'n'
  )
  text(x= pcaResult$x[,1], y= pcaResult$x[,2],col = colorCodes,labels= gsub('\\..*','',rownames(pcaResult$x)), cex= 0.7)
  dev.off()
  
  # MDS similarity analysis
  pdf(paste(prefix,'mds_libsizecpm.pdf',sep=''),12,6)
  par(mfrow=c(1,2))
  d <- DGEList(counts = y, group = group, lib.size = libSize)
  par(las= 1, cex= 0.8)
  my_labels <- group
  plotMDS(d, xaxp= c(-2, 2, 8),  yaxp= c(-2, 2, 8),
          labels= my_labels,# col= c(rep('firebrick1',8),rep('green',4), rep('orange',8),rep('darkgreen',8),rep('blue',4), rep('darkgrey',4),rep('black',4), rep('purple',8)),
          col = colorCodes,
          main= 'MDS plot, un-normalized libsize-cpm', cex= 1)
  grid()
  d <- DGEList(counts = y1, group = group, lib.size = libSize)
  par(las= 1, cex= 0.8)
  my_labels <- group
  plotMDS(d, xaxp= c(-2, 2, 8),  yaxp= c(-2, 2, 8),
          labels= my_labels,# col= c(rep('firebrick1',8),rep('green',4), rep('orange',8),rep('darkgreen',8),rep('blue',4), rep('darkgrey',4),rep('black',4), rep('purple',8)),
          col = colorCodes,
          main= 'MDS plot, dm6-normalized libsize-cpm', cex= 1)
  grid()
  dev.off()
}

# see difference between cpm normalization and full normalization distance matrices
# ** add-on: seee difference between cpm normalization with input subtration and 
# cpm normalization with input subtraction and dm6 correction **
if (do_evalute_after_subtract_input ==0) {
  amax <- max(dist_euc_tot_cpm)
  # crucial: divide matrix by maximal distance
  tmp1 <- (dist_euc_tot_cpm / amax)
  amax <- max(dist_euc_no_cpm)
  tmp2 <- (dist_euc_no_cpm / amax)}
if (do_evalute_after_subtract_input==1) {
  amax <- max(dist_euc_tot_cpm)
  # crucial: divide matrix by maximal distance
  tmp1 <- (dist_euc_tot_cpm / amax)
  amax <- max(dist_euc_cpm_input_subtract)
  tmp2 <- (dist_euc_cpm_input_subtract / amax)
}
# if (do_evalute_after_subtract_input==1) {
#   amax <- max(dist_euc_cpm_input_subtract)
#   # crucial: divide matrix by maximal distance
#   tmp1 <- (dist_euc_cpm_input_subtract / amax)
#   amax <- max(dist_euc_no_cpm)
#   tmp2 <- (dist_euc_no_cpm / amax)
# }


tmp12 <- tmp1 / tmp2
# values along the diagonal will have 0 distance, so the ratio is not a number and we put it to 1 (i.e. unchanged)
tmp12[which(is.na(tmp12))] <- 1
# crucial: scaling factor for color coding
s <- max(max(log2(tmp12)),abs(min(log2(tmp12))))

if (save_plots) {
  # plot normalized distance matrix as heatmap
  pdf(paste(prefix,'heatmap.cpm_vs_norm','.pdf',sep=''),6,6)
  par(mfrow=c(1,1))
  heatmap.2(log2(tmp12),Rowv = F, Colv = F, dendrogram = 'none', trace = 'none',cexRow = 0.17,cexCol = 0.17,breaks=seq(-s,s,length.out = 101),col=redgreen(100))
  dev.off()
}

# calculate improvement factors
# first create a mapping that codify with integer the type of group

mapping_group <- data.frame(symb=unique(group),unique_r=as.numeric(factor(unique(group))))
tmp12_labels <- data.frame(symb=group)
tmp12_labels_updated <- tmp12_labels %>% dplyr::left_join(mapping_group,by="symb")

unique_r = unique(group)
#table_r = rbind(label=unique_r, count=sapply(unique_r,function(x)sum(group==x)))
table_r = table(tmp12_labels_updated$symb)

inds_counts <- as.numeric(table_r[2,])
count <- 1;
belong <- tmp12_labels_updated$unique_r

# for (i in 1:length(inds_counts))
# {
#   aind <- count:(count+inds_counts[i]-1)
#   inds_list[[i]] <- aind
#   count <- count + length(aind)
#   belong <- c(belong,rep(i,length(aind)))
# }
inds_counts <- unique(tmp12_labels_updated$unique_r)
inds_list <- vector("list", length(inds_counts))
for (i in 1:length(inds_counts))
{
  aind <- which(tmp12_labels_updated$unique_r==inds_counts[i])
  inds_list[[i]] <- aind
}

sum_diag <- 0
sum_out_diag <- 0
improvement_cpm_biol_all <- rep(NA,length(inds_list))

for (l in 1:length(inds_list))
{
  a_mat <- log2(tmp12)[inds_list[[l]],inds_list[[l]]]
  a_mat <- a_mat[upper.tri(a_mat,diag = F)]
  improvement_cpm_biol_all[l] <- sum(log2(tmp12)[1,-inds_list[[l]]]) / (length(belong)-length(inds_list[[l]])) - sum(a_mat) / length(a_mat)
  sum_diag <- sum_diag + sum(a_mat) / length(a_mat)
  sum_out_diag <- sum_out_diag + sum(log2(tmp12)[1,-inds_list[[l]]]) / (length(belong)-length(inds_list[[l]]))
}

improvement_cpm_biol <- sum_out_diag - sum_diag # save number
if (save_norm_summary_files)
  write.table(improvement_cpm_biol_all,paste(prefix,'improvement_cpm_biol_all.txt',sep=''),col.names=F,row.names=F,quote=F,sep='\t')

improvement_cpm_biol_pairs <- rep(NA,length(inds_list))
for (i in 1:length(inds_list))
{
  tmp <- rep(NA,length(inds_list)-1)
  my_inds <- inds_list[[i]]
  a_mat <- log2(tmp12)[my_inds,my_inds]
  a_mat <- a_mat[upper.tri(a_mat,diag = F)]
  sum_diag <- sum(a_mat) / length(a_mat)
  all_other <- inds_list[-belong[i]]
  for (k in seq(1,length(all_other))){
    a_mat <- log2(tmp12)[my_inds,all_other[[k]]]
    sum_out_diag <- sum(a_mat) / (length(a_mat))
    tmp[k] <- sum_out_diag - sum_diag
  }
  improvement_cpm_biol_pairs[i] <- mean(tmp)
}
if (save_norm_summary_files)
  write.table(improvement_cpm_biol_pairs,paste(prefix,'improvement_cpm_pairs.biol.txt',sep=''),col.names=F,row.names=F,quote=F,sep='\t')
improvement_cpm_biol_pairs_avg <- mean(improvement_cpm_biol_pairs)  # save number


unique_r = unique(samples)
table_r = rbind(label=unique_r, count=sapply(unique_r,function(x)sum(samples==x)))
inds_counts <- as.numeric(table_r[2,])
count <- 1;
inds_list <- vector("list", length(inds_counts))
belong <- rep(-1,0)
for (i in 1:length(inds_counts))
{
  aind <- count:(count+inds_counts[i]-1)
  inds_list[[i]] <- aind
  count <- count + length(aind)
  belong <- c(belong,rep(i,length(aind)))
}

sum_diag <- 0
sum_out_diag <- 0
improvement_cpm_tech_all <- rep(NA,length(inds_list))
for (l in 1:length(inds_list))
{
  a_mat <- log2(tmp12)[inds_list[[l]],inds_list[[l]]]
  a_mat <- a_mat[upper.tri(a_mat,diag = F)]
  improvement_cpm_tech_all[l] <- sum(log2(tmp12)[1,-inds_list[[l]]]) / (length(belong)-length(inds_list[[l]])) - sum(a_mat) / length(a_mat)
  sum_diag <- sum_diag + sum(a_mat) / length(a_mat)
  sum_out_diag <- sum_out_diag + sum(log2(tmp12)[1,-inds_list[[l]]]) / (length(belong)-length(inds_list[[l]]))
}
improvement_cpm_tech <- sum_out_diag - sum_diag  # save number
if (save_norm_summary_files)
  write.table(improvement_cpm_tech_all,paste(prefix,'improvement_cpm_tech_all.txt',sep=''),col.names=F,row.names=F,quote=F,sep='\t')

improvement_cpm_tech_pairs <- rep(NA,length(inds_list))
for (i in 1:length(inds_list))
{
  tmp <- rep(NA,length(inds_list)-1)
  my_inds <- inds_list[[i]]
  a_mat <- log2(tmp12)[my_inds,my_inds]
  a_mat <- a_mat[upper.tri(a_mat,diag = F)]
  sum_diag <- sum(a_mat) / length(a_mat)
  all_other <- inds_list[-belong[i]]
  for (k in seq(1,length(all_other))){
    a_mat <- log2(tmp12)[my_inds,all_other[[k]]]
    sum_out_diag <- sum(a_mat) / (length(a_mat))
    tmp[k] <- sum_out_diag - sum_diag
  }
  improvement_cpm_tech_pairs[i] <- mean(tmp)
}
if (save_norm_summary_files)
  write.table(improvement_cpm_tech_pairs,paste(prefix,'improvement_cpm_pairs.tech.txt',sep=''),col.names=F,row.names=F,quote=F,sep='\t')
improvement_cpm_tech_pairs_avg <- mean(improvement_cpm_tech_pairs)  # save number

sink(paste(prefix,'improvement_factor_all.average','.txt',sep=''),append = T)
print(paste('improvement factor biological replicates = ', improvement_cpm_biol,sep=''))
print(paste('improvement factor biological replicates by pairs = ', improvement_cpm_biol_pairs_avg,sep=''))
print(paste('improvement factor technical replicates = ', improvement_cpm_tech,sep=''))
print(paste('improvement factor technical replicates by pairs = ', improvement_cpm_tech_pairs_avg,sep=''))
sink()

# barplot with improvement factors
improvement_factors_df <- data.frame(IF=c('IF_biol','IF_biol_pairs_avg','IF_tech','IF_tech_pairs_avg'),
                                     value=c(improvement_cpm_biol,improvement_cpm_biol_pairs_avg,improvement_cpm_tech,improvement_cpm_tech_pairs_avg))

improvement_factors_individ_df <- list("improvement_cpm_tech_all" = improvement_cpm_tech_all,"improvement_cpm_biol_all" = improvement_cpm_biol_all)
save(improvement_factors_individ_df,file=paste0(prefix,'list_all_improvement.RData'))
save(improvement_factors_df,file=paste0(prefix,'list_all_improvement_factors_df.RData'))

library(ggplot2)
print(print(improvement_factors_df))
g <- ggplot(improvement_factors_df, aes(x=IF,y=value)) +
  geom_bar(stat="identity",width = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylim(0,7) 
ggsave(paste0(prefix,'Improvement_factors.pdf'))

if (save_tables) {
  colnames(y_orig) <- colnames(y)
  #write.table(format(cbind(coords,cnv_orig),digits=3),paste(prefix,'cnv_counts.original_filtered.tab',sep=''),col.names=c('chr','start','end',colnames(cnv_orig)),row.names=F,sep='\t',quote=F)
  write.table(format(cbind(coords,y_orig),digits=3),paste(prefix,'peak_counts.original_filtered.tab',sep=''),col.names=c('chr','start','end',colnames(y_orig)),row.names=F,sep='\t',quote=F)
  write.table(format(cbind(coords,y2),digits=3),paste(prefix,'peak_counts.cpm_filtered.tab',sep=''),col.names=c('chr','start','end',colnames(y2)),row.names=F,sep='\t',quote=F)
  write.table(format(cbind(coords,y3),digits=3),paste(prefix,'peak_counts.cpm_norm_filtered.tab',sep=''),col.names=c('chr','start','end',colnames(y3)),row.names=F,sep='\t',quote=F)
  write.table(format(cbind(coords,y1),digits=3),paste(prefix,'peak_counts.norm_filtered.tab',sep=''),col.names=c('chr','start','end',colnames(y1)),row.names=F,sep='\t',quote=F)
  filtered_sample <- rbind(samples,group,libSize)
  write.table(filtered_sample,file=paste0(prefix,'peak_counts.cpm_filtered_sample_file.txt'),row.names=F,col.names = F,sep='\t',quote=F)
  
}
print(improvement_factors_df)