#' @title Generate sample and generate normalization factors
#' @description Generate sample file for running drosophila normalization
#' @param total number of reads (in drosophila ROI) file
#' @return sample files and normalzation factors for drosophila normalization step
#' @author angela simeone  \email{angela.simeone@gmail.com}


## =========================================================
## ===  GENERATE SAMPLE FILE AND PRINT TO FILE === 
## =========================================================

rm(list = ls())

initial_file <- read.table('tot_num_reads_hg19_libraries.txt',stringsAsFactors = F)

initial_file <- as.data.frame(initial_file)
colnames(initial_file) <- c('lib','tot_reads')

initial_file$sample_name <- gsub('_R.*','',initial_file$lib)
initial_file$sample_name <- gsub('_ChIP.*','',initial_file$sample_name)

initial_file$sample_name <- gsub('_input.*.','_input',initial_file$sample_name)

# # normal cases
# normal_cases <-"MCF10A|ATCC|ZenBioA|MLE-F_B"
# initial_file$group <- initial_file$sample_name

# for each PDTX greate a group label
PDTX_input_temp <- gsub('.all.hg19','',initial_file$lib)
PDTX_input=c('139M_181_input','139M_284_input','143_284_input','143_317_input','AB636_input','AB577_input','179_181_input','179_284_input','PAR1022_input','201_181_input','201_284_input','AB580_input','316_284_input','331_284_input','521_284_input','5_317_input','863_317_input','9_317_input','98_181_input','98_284_input','STG139_input','AB555_input','AB551_input','STG282_input','AB790_input','STG195_input','PAR1006_3_input','STG197_3_input')
PDTX_pd <- gsub("_input","",PDTX_input)

PDTX_all <- c(PDTX_pd,PDTX_input)

initial_file$sample <-initial_file$sample_name

for (i in (1:length(PDTX_all))) {
  curr_label <- PDTX_all[i]
  initial_file$sample[grep(curr_label,initial_file$sample_name)] <- curr_label
}


for (i in (1:length(PDTX_pd))) {
  curr_label <- PDTX_pd[i]
  initial_file$group[grep(curr_label,initial_file$sample_name)] <- curr_label
}

initial_file$input <- initial_file$group

for (i in (1:length(PDTX_pd))) {
  curr_label <- PDTX_all[i]
  curr_label_input <- PDTX_input[i]
  initial_file$input[grep(curr_label,initial_file$sample_name)] <- curr_label_input
}


#flag: 0 - include ; 1 - exclude
initial_file$flag <- rep(0,length(initial_file$group))
initial_file$flag[grep("STG197",initial_file$group)] <- 1

# some groups represent the same group
initial_file$group[grep("139M_181|139M_284",initial_file$group)] <- "139M_"
initial_file$group[grep("143_284|143_317",initial_file$group)] <- "143_"
initial_file$group[grep("201_181|201_284",initial_file$group)] <- "201_"
initial_file$group[grep("98_181|98_284",initial_file$group)] <- "98_"
initial_file$group[grep("179_181|179_284",initial_file$group)] <- "179_"

# create a flag for cells and tissues
initial_file$tissueCells_flag <- rep(0,length(initial_file$group))
#initial_file$tissueCells_flag[grep('Normal|MCF7|MDAMB23|MLE',initial_file$group)] <- 1

# reorder columns
initial_file_ordered <-cbind(initial_file$sample,initial_file$group,initial_file$tot_reads,initial_file$flag,initial_file$input)

write.table(initial_file_ordered,'tot_num_reads_hg19_libraries_intermediate1.txt',quote=  F,sep = "\t",col.names = F,row.names = F)

write.table(t(initial_file_ordered),'sample_file_PDTX_only.txt',quote=  F,sep = "\t",col.names = F,row.names = F)


## =========================================================
## ===  GENERATE NORMALIZATION FACTORS AND PRINT TO FILE === 
## =========================================================

## norm_factors from drosophila - paired to sample file
dm6_stats <- read.table('dm6_stats_with_fnames.txt',stringsAsFactors = F)
colnames(dm6_stats) <- c('lib','tot','peak')
dm6_stats$fRIP <- dm6_stats$peak/dm6_stats$tot

#paste samples file and dm6stats
sample_matrix_and_dm6 <- cbind(initial_file,dm6_stats)



#exclude inputs and also experiments with flag==1
inds_to_exclude <- union(grep('input',sample_matrix_and_dm6$sample_name),which(sample_matrix_and_dm6$flag==1))

sample_matrix_and_dm6_max_peak <- max(sample_matrix_and_dm6$fRIP[-inds_to_exclude])
sample_matrix_and_dm6_max_tot <- max(sample_matrix_and_dm6$tot_reads[-inds_to_exclude])
sample_matrix_and_dm6[which(sample_matrix_and_dm6$fRIP==max(sample_matrix_and_dm6$fRIP[-inds_to_exclude])),]
sample_matrix_and_dm6[which(sample_matrix_and_dm6$tot==max(sample_matrix_and_dm6$tot[-inds_to_exclude])),]

norm_peak <- sample_matrix_and_dm6_max_peak/sample_matrix_and_dm6$fRIP
norm_total <- sample_matrix_and_dm6_max_tot/sample_matrix_and_dm6$tot_reads

sample_matrix_and_dm6$norm_peak <- norm_peak
sample_matrix_and_dm6$norm_total <- norm_total

sample_matrix_and_dm6$norm_peak[inds_to_exclude] <- 1
sample_matrix_and_dm6$norm_total[inds_to_exclude] <- 1


write.table(sample_matrix_and_dm6,file='sample_matrix_and_dm6.temp.txt',quote=  F,sep = "\t",col.names = NA)
write.table(t(sample_matrix_and_dm6$norm_peak),'PDTX_old_new_dm6_norm5_peak_recovery.txt',quote=  F,sep = "\t",col.names = F,row.names = F)
write.table(t(sample_matrix_and_dm6$norm_total),'PDTX_old_new_dm6_total_recovery.txt',quote=  F,sep = "\t",col.names = F,row.names = F)
write.table(rep(1,length(sample_matrix_and_dm6$norm_total)),'PDTX_old_new_no_dm6.txt',quote=  F,sep = "\t",col.names = F,row.names = F)
