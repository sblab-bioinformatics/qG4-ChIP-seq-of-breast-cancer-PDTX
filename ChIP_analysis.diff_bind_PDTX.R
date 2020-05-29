#' @title differential binding analysis
#' @description perform iterative differential binding analysis
#' @param processed drosophila and input-subtracted counts
#' @return differential binding sites for each PDTX
#' @author angela simeone  \email{angela.simeone@cruk.cam.ac.uk}


# CAREFUL: output_dir_name and base_dir DEFINED IN PREVIOUS SCRIPT ChIP_analysis.normalize_data.R
base_dir <- '/PDTX_dir'

setwd(base_dir)

output_dir_name <- 'PDTX_dm6_inputSubtraction_/'

adir <- paste(base_dir,output_dir_name,sep='')
setwd(adir)
log_file_name <- paste(gsub('\\/$','',adir),'.log',sep='')

library(gplots)
library(edgeR)
library(dendextend)
library(RColorBrewer)

# ====================
# == INPUT FILES =====
# ====================

# tab-separated file file containing peak file and read counts in each peak, for each samples: see peak_example.bedgraph
processed_peak_file_counts_name <- paste(adir,'peak_counts.norm_filtered.tab',sep='')

# tab-separated file containing sample information see sample_info_example.txt
sample_file_name <- "/PDTX_dir/sample_file_PDTX_only.txt"

# tab-separated file containing the normalization factor for each sample: norm_example.txt
output_prefix <- paste(tools::file_path_sans_ext(processed_peak_file_counts_name),'.detable',sep='')


# ==============================
# ===== CONTROL PARAMETERS =====
# ==============================
save_plots <- 0
save_tables <- 1 
th <- 0.01
FC <- 0.6
de_type <- 0 # 0: split; 1: merged; 2:batch
#de_type <- 2# 0: split; 1: merged; 2:batch

# ====== Load input files ======
peaks <- read.table(processed_peak_file_counts_name,stringsAsFactors = F,header=T,sep='\t')
colnames(peaks) <- gsub('^X','',colnames(peaks))
peaks <- peaks[-grep('chr.*._.*',peaks[,1]),] # to exclude those in chr_

sample_file <- read.table(sample_file_name,stringsAsFactors = F,header=F,sep='\t')

samples <- as.character(sample_file[1,])
group <- as.character(sample_file[2,])
batch <- as.character(sample_file[6,])
libSize <- as.numeric(sample_file[3,])
y <- peaks[,-c(1,2,3)]
coords <- peaks[,c(1,2,3)]
row.names(y) <- paste(peaks[,1],peaks[,2],peaks[,3],sep='_')

inds <- rep(-1,0)
if (dim(sample_file)[1] >= 4)
  inds <- c(inds,which(sample_file[4,] == 1))
inds <- c(inds,grep('input',samples))
sink(log_file_name,append = T)
print('Input samples and file removed:')
samples[sort(inds)]
sink()

samples <- samples[-inds]
group <- group[-inds]
libSize <- libSize[-inds]
batches <- batch[-inds]

# ======== edgeR differential binding analysis ===========
if (de_type == 0) {
  # Combine all the experimental factors (tumour and date) into one combined factor (3.3.1 of edgeR user guide)
  de<- DGEList(counts=y, group=samples,lib.size = libSize)
  count_per_million=cpm(y,lib.size = libSize)
  agroup <- as.character(samples)
  design <- model.matrix(~0+agroup)
  conds_b <- as.character(group[which(!duplicated(group))])
  #conds <- as.character(samples[which(!duplicated(samples))])
  conds <- as.character(colnames(design))
}
if (de_type == 1) {
  # Combine all tumours in one factors (flatten date information)
  de<- DGEList(counts=y, group=group,lib.size = libSize)
  agroup <- as.character(group)
  design <- model.matrix(~0+agroup)
  conds_b <- as.character(group[which(!duplicated(group))])
  conds <- conds_b
}
if (de_type == 2) {
  # Model explicitely all the experimental factors (tumour and date) into a linear model (3.4.2 of edgeR user guide)
  de<- DGEList(counts=y, group=group,lib.size = libSize)
  de<- calcNormFactors(de, method= 'none')
  Batch <- factor(batches)
  agroup <- as.character(samples)
  design <- model.matrix(~0+agroup+Batch)
  
  # this case seems to generate 
  # Error in glmFit.default(y, design = design, dispersion = dispersion, offset = offset,  : 
  #                           Design matrix not of full rank.  The following coefficients not estimable:
  #                           Batchwk - 
  #Tumour <- factor(group)
  #design <- model.matrix(~0+Tumour+Batch)
  
  
  conds_b <- as.character(group[which(!duplicated(group))])
  conds <- conds_b
}

de<- estimateGLMCommonDisp(de, design)
de<- estimateGLMTagwiseDisp(de, design)
fit<- glmFit(de, design)

#colnames(design)[1:length(conds)] <- conds ##here the sorting is not necessarely right!!
#colnames(design)[1:length(conds)] <- gsub('agroup','',colnames(design)[1:length(conds)])

n <- length(unique(group))
detables <- vector("list", n) 
names(detables) <- conds_b

peaks_up_down <- matrix(0, dim(y)[1],n)
colnames(peaks_up_down) <- conds_b
row.names(peaks_up_down) <- row.names(y)

sink(log_file_name,append = T)
for (i in 1:length(detables))
{
  print(conds_b[i])
  mycond <- conds_b[i]
  # inds_yes <- grep(mycond,conds) ############## i have deleted this
  # inds_no <- setdiff(seq(1:length(conds)),c(inds_yes))############## i have deleted this
  inds_yes <- grep(mycond,colnames(design))
  inds_no <- setdiff(seq(1:dim(design)[2]),c(inds_yes))
  
  print(paste('Current condition being tested for dba: ',mycond,sep=''))
  print(paste('Indexes belonging: ',paste(inds_yes,collapse=' '),sep=''))
  conds[inds_yes]
  print(paste('Indexes NOT belonging: ',paste(inds_no,collapse=' '),sep=''))
  conds[inds_no]
  contrasts <- rep(0,dim(design)[2])
  contrasts[inds_yes] <- 1/length(inds_yes)
  contrasts[inds_no] <- -1/length(inds_no)
  print(contrasts)
  # contrasts[inds_yes] <- 1
  # contrasts[inds_no] <- -1
  
  print(paste('Contrast vector: ',paste(format(as.numeric(contrasts),digits=3),collapse=' '),sep=''))
  lrt <- glmLRT(fit,contrast=contrasts)
  detable<- topTags(lrt, n= nrow(y))$table
  length(which(detable$logFC>0.6 & detable$FDR<= 0.05))
  detable$ID <- rownames(detable)
  up <- detable$ID[which(detable$logFC >= FC & detable$FDR <= th)]
  down <- detable$ID[which(detable$logFC <= -FC & detable$FDR <= th)]
  peaks_up_down[which(row.names(peaks_up_down) %in% up),i] <- 1
  peaks_up_down[which(row.names(peaks_up_down) %in% down),i] <- -1
  detables[[i]] <- detable
  
}
sink()

if (save_tables)
{
  for (i in 1:length(detables))
    write.table(detables[[i]],paste(output_prefix,'.',names(detables[i]),'_vs_all.tab',sep=''),col.names=colnames(detables[[i]]),row.names=F,sep='\t')
  
  peaks_up_down_df <- as.data.frame(peaks_up_down,stringsAsFactors = F)
  peaks_up_down_df$ID <- row.names(peaks_up_down_df)
  write.table(peaks_up_down_df,paste(output_prefix,'.peaks_table_one_vs_all.tab',sep=''),col.names=colnames(peaks_up_down_df),row.names=F,sep='\t')
  save(detables,peaks_up_down,file=paste(output_prefix,'.detables_one_vs_all.RData',sep=''))
  
  # save also bed file for peak summary
  peaks_up_down_df <- as.data.frame(cbind(peaks[,1:3],peaks_up_down),stringsAsFactors=F)
  write.table(peaks_up_down_df,paste(output_prefix,'.peaks_table_one_vs_all.bed',sep=''),col.names=colnames(peaks_up_down_df),row.names=F,sep='\t')
}

pdf(paste(output_prefix,'.diff_peaks_n_barplot_nonorm.pdf',sep=''),12,6)
par(mfrow=c(1,3))
bars <- rep(0,n)
for (i in 1:n)
  bars[i] <- length(which(peaks_up_down[,i] == 1))
barplot(bars,names=conds_b,las=2,ylab='# sepcific up-regulated peaks')
bars2 <- rep(0,n)
for (i in 1:n)
  bars2[i] <- length(which(peaks_up_down[,i] == -1))
barplot(bars2,names=conds_b,las=2,ylab='# sepcific downs-regulated peaks')
barplot(rbind(bars,bars2),names=conds_b,las=2)
legend(1,10000,c('up','down'),text.col = grey.colors(2),fill=grey.colors(2))
dev.off()

all <- (rowSums(abs(peaks_up_down)))
table(all)
inv <- row.names(peaks_up_down)[which(all == 0)]
inv_means <- rowMeans(y[which(row.names(y) %in% inv),])
noinv_means <- rowMeans(y[which(!row.names(y) %in% inv),])

pdf(paste(output_prefix,'.invariant_boxplot_nonorm.pdf',sep=''),4,8)
par(mfrow=c(2,1))
boxplot(log2(inv_means),log2(noinv_means),outline=F,names=c('invariant','variant'),ylab='all tumors mean signal (log2)')
boxplot((inv_means),(noinv_means),outline=F,names=c('invariant','variant'),ylab='all tumors mean signal')
dev.off()


# =====================================
# ===== in bash do the following: =====
# =====================================

## in bash do the following:
#
#for f in `ls *_vs_all.tab| grep -v table_one`; do cat $f | awk -F "\t" '{ if(($1 >= 0.6) && ($5 < 0.05)) { print } }' | awk -F "\t" '{print $6}'|  sed  's/_/     /g' | sed 's/"//g' |awk '{print $1"\t"$2"\t"$3}' |  sortBed -i - > ${f%%.tab}_up.bed; done
#for f in `ls *_vs_all.tab| grep -v table_one`; do cat $f | awk -F "\t" '{ if($1 >= 0.6 && $5 < 0.05) print}' | awk -F "\t" '{print $6}'|  sed  's/_/     /g' | sed 's/"//g' |awk '{print $1"\t"$2"\t"$3}' |  sortBed -i - > ${f%%.tab}_up.bed; done


# for f in `ls *_vs_all.tab| grep -v table_one`; 
# do 
# cat $f | awk -F "\t" '{ if(($1 >= 1.5) && ($5 < 0.05)) { print } }' | awk -F "\t" '{print $6}'|  sed  's/_/     /g' | sed 's/"//g' |awk '{print $1"\t"$2"\t"$3}' |  sortBed -i - > ./fc_2/${f%%.tab}.1.5_up.bed; 
# done



## in bash generate bed files with also logcpm on the 4th columnfor 
#for f in `ls *_vs_all.tab| grep -v table_one`; do cat $f | awk -F "\t" '{ if(($1 >= 1) && ($5 < 0.05)) { print } }' | awk -F "\t" '{print $6}'|  sed  's/_/     /g' | sed 's/"//g' |awk '{print $1"\t"$2"\t"$3}' |  sortBed -i - > ${f%%.tab}_up_lfc1.bed; done




# 
# # =================================
# # ===== PAIR BY PAIR ANALYSIS =====
# # =================================
# # 1) for each pair (identified with naming automatically)
# # 2) do diff expr
# # 3) if genes is up write peak name in the list relative to the given condition (if it's down write in the other list)
# detables_pairs <- vector("list", (n**2 - n) / 2)
# detables_pairs_names <- rep('', (n**2 - n) / 2)
# all_up <- vector("list", n)
# all_peak_tables <- vector("list", n)
# all_up_unique <- vector("list", n)
# peaks_up_down_unique <- matrix(0, dim(y[1]),n)
# peaks_up_down_count <- matrix(0, dim(y[1]),n)
# row.names(peaks_up_down_unique) <- row.names(y)
# row.names(peaks_up_down_count) <- row.names(y)
# 
# sink(log_file_name,append = T)
# count <- 1
# for (i in 1:(length(conds_b)-1)) {
#   for (k in (i+1):length(conds_b)) {
#     cond1 <- conds_b[i]; cond2 <- conds_b[k]
#     detables_pairs_names[count] <- paste(cond1,'vs',cond2,sep='_')
#     print (paste(i,"versus",k))
#     inds_yes <- grep(cond1,conds); inds_no <- grep(cond2,conds)
#     contrasts <- rep(0,dim(design)[2])
#     contrasts[inds_yes] <- 1/length(inds_yes)
#     contrasts[inds_no] <- -1/length(inds_no)
# 
#     print(paste('Currecnt condition being tested for dba: ',cond1,' vs ',cond2,sep=''))
#     print(paste('Indexes belonging to condition_1: ',paste(inds_yes,collapse=' '),sep=''))
#     print(paste('Indexes belonging to condition_2: ',paste(inds_no,collapse=' '),sep=''))
#     print(paste('Contrast vector: ',paste(format(as.numeric(contrasts),digits=3),collapse=' '),sep=''))
# 
#     lrt <- glmLRT(fit,contrast=contrasts)
#     detable<- topTags(lrt, n= nrow(y))$table
#     detable$ID <- rownames(detable)
#     up <- detable$ID[which(detable$logFC >= FC & detable$FDR <= th)]
#     down <- detable$ID[which(detable$logFC <= -FC & detable$FDR <= th)]
#     all_up[[i]] <- c(all_up[[i]],up)
#     all_up[[k]] <- c(all_up[[k]],down)
#     detables_pairs[[count]] <- detable
#     count <- count + 1
#   }
# }
# names(detables_pairs) <- detables_pairs_names
# sink()
# 
# # 4) after all pairs have been done, take unique of peak name per tumor and fill similar table to peaks_up_down_n
# #     1 will mean that peak has figured has up-regulated in that condition when compared to at least another condition
# #     (alternative count how many times in pair comparisons and write the number here instead on 1)
# tmp1 <- as.data.frame(peaks_up_down_count,stringsAsFactors = F)
# colnames(tmp1) <- conds_b
# tmp1$id <- row.names(peaks_up_down_count)
# for (i in 1:length(conds_b)) {
#   all_up_unique[[i]] <- unique(all_up[[i]])
#   all_peak_tables[[i]] <- table(all_up[[i]])
#   peaks_up_down_unique[which(row.names(peaks_up_down_unique) %in% all_up_unique[[i]]),i] <- 1
#   tmp2 <- as.data.frame(all_peak_tables[[i]],stringsAsFactors = F)
#   colnames(tmp2) <- c('Var1',paste(conds_b[i],'_count',sep=''))
#   tmp1 <- merge(tmp1,tmp2,by.x='id',by.y='Var1',all=T)
#   tmp1[which(!is.na(tmp1[,n+2])),(i+1)] <- tmp1[which(!is.na(tmp1[,n+2])),n+2]
#   tmp1 <- tmp1[-(n+2)]
# }
# tmp1$id <- factor(tmp1$id, levels=row.names(peaks_up_down_unique))
# output <- tmp1[order(tmp1$id),]
# peaks_up_down_count <- output[,2:(n+1)]
# row.names(peaks_up_down_count) <- output[,1]
# 
# # save to file: peaks_up_down_count, peaks_up_down_unique;
# 
# if (save_tables)
# {
#   for (i in 1:length(detables_pairs))
#     write.table(detables_pairs[[i]],paste(output_prefix,'.',names(detables_pairs[i]),'.tab',sep=''),col.names=colnames(detables_pairs[[i]]),row.names=F,sep='\t')
#   peaks_up_down_count_df <- as.data.frame(peaks_up_down_count,stringsAsFactors = F)
#   peaks_up_down_count_df$ID <- row.names(peaks_up_down_count_df)
#   write.table(peaks_up_down_count_df,paste(output_prefix,'.peaks_table_pairs.tab',sep=''),col.names=colnames(peaks_up_down_count_df),row.names=F,sep='\t')
#   save(detables_pairs,peaks_up_down_count,file=paste(output_prefix,'.detables_one_vs_one.RData',sep=''))
#   colnames(peaks_up_down_unique) <- colnames(peaks_up_down_count)
#   peaks_up_down_unique_df <- as.data.frame(peaks_up_down_unique,stringsAsFactors = F)
#   peaks_up_down_unique_df$ID <- row.names(peaks_up_down_unique_df)
#   write.table(peaks_up_down_unique_df,paste(output_prefix,'.peaks_table_pairs_unique.tab',sep=''),col.names=colnames(peaks_up_down_count_df),row.names=F,sep='\t')
#   # save also bed file for peak summary
#   peaks_up_down_count_df <- as.data.frame(cbind(peaks[,1:3],peaks_up_down_count),stringsAsFactors=F)
#   write.table(peaks_up_down_count_df,paste(output_prefix,'.peaks_table_one_vs_one_count.bed',sep=''),col.names=colnames(peaks_up_down_count_df),row.names=F,sep='\t')
#   peaks_up_down_unique_df <- as.data.frame(cbind(peaks[,1:3],peaks_up_down_unique),stringsAsFactors=F)
#   write.table(peaks_up_down_unique_df,paste(output_prefix,'.peaks_table_one_vs_one_unique.bed',sep=''),col.names=colnames(peaks_up_down_unique_df),row.names=F,sep='\t')
# }
# # system('gzip *.tab')
# # system('gzip *.bed')