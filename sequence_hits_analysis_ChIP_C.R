#!/usr/bin/env Rscript --vanilla

#This script takes in input the bed file of interest and 3 randomized cases (for now... todo: extend to variable number of shuffled cases)

# from Giovanni strategy (illustrated in Robert paper)
#he total number of G4 ChIP-seq regions for HaCaT cells in each of the indicated structural classes of G4s are shown.
# - Loop sizes 1–3, 4–5 and 6–7 indicate that at least one loop of this length is present in the G4; 
# - a long loop indicates a G4 with any loop of length > 7 (up to 12 for any loop and 21 for the middle loop); 
# - a simple bulge indicates a G4 with a bulge of 1–7 bases in one G-run or multiple 1-base bulges; 
# - a 2-tetrads/complex bulge indicates G4s with two G-bases per G-run or several bulges of 1–5 bases; 
# - and other indicates other sequences that do not fall into the former categories. 

# #setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/fasta_G4_regions')
args <- commandArgs(trailingOnly = TRUE)
# args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
file1 <- args[2]
file2 <- args[3]
file3 <- args[4]

# setwd('/Users/simeon01/Documents/PDTX/NOVEMBER_2019_12PDTX/consensus_hg19/')
# file <- 'hg19_old_new_q005.all_peaks.25M.over99nt.sorted.fa'
# file1 <- 'hg19_old_new_q005.all_peaks.25M.over99nt.sorted.shuffle_1.fa'
# file2 <- 'hg19_old_new_q005.all_peaks.25M.over99nt.sorted.shuffle_2.fa'
# file3 <- 'hg19_old_new_q005.all_peaks.25M.over99nt.sorted.shuffle_3.fa'

file_out <- gsub('.fa','.results_motifs.Rdata',file)
print(file_out)
# chr_list <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
#               'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
#               'chr21','chr22','chrX','chrY','chrM')
# 
# strands = c('minus','plus')
# # reading information
# nc <- length(chr_list) * 2
nc <- 1

match_me <- function(aseq, a_motif, i) {
  ind <- integer()
  if(length(grep(a_motif,aseq))>0)
    ind <- i
  return(ind)
}

sequences_from_fasta <- function(fasta_file){
  temp_fasta <- read.table(fasta_file, stringsAsFactors=FALSE, sep='\t', header=FALSE)
  fasta_seq <- temp_fasta$V1[seq(2,length(temp_fasta$V1),2)]
  return(fasta_seq)
}

bed_from_fasta <- function(fasta_file){
  temp_fasta <- read.table(fasta_file, stringsAsFactors=FALSE, sep='\t', header=FALSE)
  temp_bed_intervals <- temp_fasta$V1[seq(1,length(temp_fasta$V1),2)]
  chr_temp <- sapply(strsplit(gsub(':','-',gsub('>','',temp_bed_intervals),"-"),split = "-"),"[",1)
  start_temp <- sapply(strsplit(gsub(':','-',gsub('>','',temp_bed_intervals),"-"),split = "-"),"[",2)
  end_temp <- sapply(strsplit(gsub(':','-',gsub('>','',temp_bed_intervals),"-"),split = "-"),"[",3)
  bed_intervals <- cbind(unlist(chr_temp),unlist(start_temp),unlist(end_temp))
  return(bed_intervals)
}


files = rep('',nc)
counts_7 = counts_r1_7 = counts_r2_7 = counts_r3_7 = rep(-1,nc)
counts_long = counts_r1_long = counts_r2_long = counts_r3_long = rep(-1,nc)
counts_bulges = counts_r1_bulges = counts_r2_bulges = counts_r3_bulges = rep(-1,nc)
counts_sbulges = counts_r1_sbulges = counts_r2_sbulges = counts_r3_sbulges = rep(-1,nc)
counts_cbulges = counts_r1_cbulges = counts_r2_cbulges = counts_r3_cbulges = rep(-1,nc)
counts_cbulges_GG = counts_r1_cbulges_GG = counts_r2_cbulges_GG = counts_r3_cbulges_GG = rep(-1,nc)
counts_pqs = counts_r1_pqs = counts_r2_pqs = counts_r3_pqs = rep(-1,nc)
counts_fraction = counts_r1_fraction = counts_r2_fraction = counts_r3_fraction = rep(-1,nc)
counts_5 = counts_r1_5 = counts_r2_5 = counts_r3_5 = rep(-1,nc)
counts_3 = counts_r1_3 = counts_r2_3 = counts_r3_3 = rep(-1,nc)
counts_GG = counts_r1_GG =counts_r2_GG  = counts_r3_GG=rep(-1,nc);
counts_other = counts_r1_other = counts_r2_other = counts_r3_other = rep(-1,nc)
inds_others = rep(-1,nc) # to get the inds of those regions that do not show any queried structure

counts_A = counts_r1_A = counts_r2_A = counts_r3_A = vector("list", nc)
counts_T = counts_r1_T = counts_r2_T = counts_r3_T = vector("list", nc)
counts_C = counts_r1_C = counts_r2_C = counts_r3_C = vector("list", nc)
counts_G = counts_r1_G = counts_r2_G = counts_r3_G = vector("list", nc)
counts_AT = counts_r1_AT = counts_r2_AT = counts_r3_AT = vector("list", nc)
counts_CG = counts_r1_CG = counts_r2_CG = counts_r3_CG = vector("list", nc)
counts_skewGC = counts_r1_skewGC = counts_r2_skewGC = counts_r3_skewGC = vector("list", nc)
counts_skewAT = counts_r1_skewAT = counts_r2_skewAT = counts_r3_skewAT = vector("list", nc)

motif_GG <- '(G{2,}[ATC]{1,7}){3,}G{2,}'
motif_GGATG_1 <- 'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[AT]{1,7}G|G[AT]{1,7}GG)'

lens = rep(-1,nc);
count_f = 0

#for (c in chr_list) {
#  for (s in strands) {
count_f = count_f + 1;
#chr = c
#strand = s
print(file)
files[count_f] = file

chr_intersect_seq <-sequences_from_fasta(file)
bed_intervals_file <-bed_from_fasta(file)
# random suffling should ideally be in open chromatin with OQs for ChIP-seq motif analysis

chr_random1_seq <- sequences_from_fasta(file1)
chr_random2_seq <- sequences_from_fasta(file2)
chr_random3_seq <- sequences_from_fasta(file3)

# explore actual intersect sequence
seqs <- list(chr_intersect_seq, chr_random1_seq, chr_random2_seq, chr_random3_seq);  
seqs_prints <- c('chr_intersect_seq', 'chr_random1_seq', 'chr_random2_seq', 'chr_random3_seq');  

# for each fasta file
for (s1 in 1:length(seqs)){  
  
  print(seqs_prints[s1])
  seq <- toupper(seqs[[s1]])
  
  countsA = countsT = countsC = countsG = countsAT = countsCG = skewGC = skewAT = rep(0, length(seq))
  count_A7 = count_C7 = count_G7 = count_T7 = 0
  
  inds_A7 = inds_C7 = inds_G7 = inds_T7 = rep(1,0)
  inds_12 = inds_7 = inds_5 = inds_3 = inds_C7 = inds_21 = rep(1,0)
  inds_GG = inds_GGATG_1 = inds_GGATG_2 = inds_GGATG_3 = inds_GGATG_4 = rep(1,0)
  inds_GGATG_1_2 = inds_GGATG_1_3 = inds_GGATG_1_4 = inds_GGATG_2_3 = inds_GGATG_2_4 = inds_GGATG_3_4 = rep(1,0)
  inds_GGATG_1_2_3 = inds_GGATG_1_2_4 = inds_GGATG_2_3_4 = inds_GGATG_1_3_4 = rep(1,0)
  inds_GGATG_1_1 = inds_GGATG_2_2 =  inds_GGATG_3_3 =  inds_GGATG_4_4 = rep(1,0)
  inds_GGATG_1n_1n = inds_GGATG_2n_2n = inds_GGATG_3n_3n = inds_GGATG_4n_4n = rep(1,0)
  inds_GGATG_1n_2n = inds_GGATG_1n_3n = inds_GGATG_1n_4n = inds_GGATG_2n_3n = inds_GGATG_2n_4n = inds_GGATG_3n_4n = rep(1,0)
  inds_GGATG_1n_2n_3n = inds_GGATG_1n_2n_4n = inds_GGATG_1n_3n_4n = inds_GGATG_2n_3n_4n = rep(1,0)
  
  # for each sequence in fasta file
  for (i in 1:length(seq))
  {
    flag <- TRUE
    aseq <- seq[i]  
    
    tmp <- strsplit(aseq,'')[[1]]
    len <- length(tmp)/100
    countsT[i] <- length(which(tmp == 'T'))/len
    countsA[i] <- length(which(tmp == 'A'))/len
    countsC[i] <- length(which(tmp == 'C'))/len
    countsG[i] <- length(which(tmp == 'G'))/len
    countsCG[i] <- (countsC[i] + countsG[i])/len
    countsAT[i] <- (countsA[i] + countsT[i])/len
    skewGC[i] <- (countsG[i] - countsC[i])/(countsG[i] + countsC[i])
    skewAT[i] <- (countsA[i] - countsT[i])/(countsA[i] + countsT[i])
    
    if(i%%10000 == 0)
      print(i)
    if (length(grep('A{7,}',aseq))>0) 
    {count_A7 = count_A7 + 1 ; inds_A7= c(inds_A7,i);}
    if (length(grep('T{7,}',aseq))>0) 
    {count_T7 = count_T7 + 1 ; inds_T7= c(inds_T7,i);}
    if (length(grep('C{7,}',aseq))>0) 
    {count_C7 = count_C7 + 1 ; inds_C7= c(inds_C7,i);}
    if (length(grep('G{7,}',aseq))>0) 
    {count_G7 = count_G7 + 1 ; inds_G7= c(inds_G7,i);}
    
    # G-rich strand        
    if (length(grep('(G{3,}[ATGC]{1,3}){3,}G{3,}',aseq))>0) 
    {inds_3= c(inds_3,i); flag<-FALSE}
    if (length(grep('(G{3,}[ATGC]{1,5}){3,}G{3,}',aseq))>0) 
    {inds_5= c(inds_5,i); flag<-FALSE}
    if (length(grep('(G{3,}[ATGC]{1,7}){3,}G{3,}',aseq))>0) 
    {inds_7= c(inds_7,i); flag<-FALSE}
    if (length(grep('(C{3,}[ATGC]{1,7}){3,}C{3,}',aseq))>0) 
    {inds_C7= c(inds_C7,i); }
    if (length(grep('(G{3,}[ATGC]{1,12}){3,}G{3,}',aseq))>0) 
    {inds_12= c(inds_12,i); flag<-FALSE}
    if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{13,21}G{3,}[ATGC]{1,7}G{3,}',aseq))>0) 
    {inds_21= c(inds_21,i); flag<-FALSE}
    
    # C-rich strand        
    if (length(grep('(C{3,}[ATGC]{1,3}){3,}C{3,}',aseq))>0) 
    {inds_3= c(inds_3,i); flag<-FALSE}
    if (length(grep('(C{3,}[ATGC]{1,5}){3,}C{3,}',aseq))>0) 
    {inds_5= c(inds_5,i); flag<-FALSE}
    if (length(grep('(C{3,}[ATGC]{1,7}){3,}C{3,}',aseq))>0) 
    {inds_7= c(inds_7,i); flag<-FALSE}
    if (length(grep('(C{3,}[ATGC]{1,12}){3,}C{3,}',aseq))>0) 
    {inds_12= c(inds_12,i); flag<-FALSE}
    if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{13,21}C{3,}[ATGC]{1,7}C{3,}',aseq))>0) 
    {inds_21= c(inds_21,i); flag<-FALSE}
    
    if (flag)
    {
      # G-rich strand
      # inds_GG = c(inds_GG,match_me(aseq, motif_GG, i))
      # inds_GGATG_1 = c(inds_GGATG_1,match_me(aseq, motif_GGATG_1, i))
      
      if(length(grep('(G{2,}[ATCG]{1,7}){3,}G{2,}',aseq))>0) {
        inds_GG = c(inds_GG,i) }
      if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{1,7}G|G[ATC]{1,7}GG)',aseq))>0){
        inds_GGATG_1 = c(inds_GGATG_1,i) }
      if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{1,7}G|G[ATC]{1,7}GG)[ATGC]{1,7}G{3,}',aseq))>0){
        inds_GGATG_2 = c(inds_GGATG_2,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(GG[ATC]{1,7}G|G[ATC]{1,7}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}',aseq))>0){
        inds_GGATG_3 = c(inds_GGATG_3,i) }
      if (length(grep('(GG[ATC]{1,7}G|G[ATC]{1,7}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_4 = c(inds_GGATG_4,i) }
      if (length(grep('(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_1_2 = c(inds_GGATG_1_2,i) }
      if (length(grep('(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_1_3= c(inds_GGATG_1_3,i) }
      if (length(grep('(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)',aseq))>0) {
        inds_GGATG_1_4= c(inds_GGATG_1_4,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}',aseq))>0){
        inds_GGATG_2_3 = c(inds_GGATG_2_3 ,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)',aseq))>0){
        inds_GGATG_2_4 = c(inds_GGATG_2_4 ,i) }
      if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)',aseq))>0) {
        inds_GGATG_3_4 = c(inds_GGATG_3_4 ,i) }
      if (length(grep('(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_1_2_3 = c(inds_GGATG_1_2_3 ,i) }
      if (length(grep('(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)',aseq))>0) {
        inds_GGATG_1_2_4 = c(inds_GGATG_1_2_4 ,i) }
      if (length(grep('(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)',aseq))>0) {
        inds_GGATG_1_3_4 = c(inds_GGATG_1_3_4 ,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)[ATGC]{1,7}(GG[ATC]G|G[ATC]GG)',aseq))>0) {
        inds_GGATG_2_3_4 = c(inds_GGATG_2_3_4 ,i) }
      if (length(grep('(G[ATC]G[ATC]C)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_1_1 = c(inds_GGATG_1_1 ,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(G[ATC]G[ATC]C)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_2_2 = c(inds_GGATG_2_2,i) }
      if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(G[ATC]G[ATC]C)[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_3_3 = c(inds_GGATG_3_3,i) }
      if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(G[ATC]G[ATC]C)',aseq))>0) {
        inds_GGATG_4_4 = c(inds_GGATG_4_4,i) }
      if (length(grep('(G[ATC]{2,5}G[ATC]{2,5}C)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_1n_1n = c(inds_GGATG_1n_1n,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(G[ATC]{2,5}G[ATC]{2,5}C)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_2n_2n = c(inds_GGATG_2n_2n,i) }
      if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(G[ATC]{2,5}G[ATC]{2,5}C)[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_3n_3n = c(inds_GGATG_3n_3n,i) }
      if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(G[ATC]{2,5}G[ATC]{2,5}C)',aseq))>0) {
        inds_GGATG_4n_4n = c(inds_GGATG_4n_4n,i) }
      if (length(grep('(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_1n_2n = c(inds_GGATG_1n_2n,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_2n_3n = c(inds_GGATG_2n_3n,i) }
      if (length(grep('(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_1n_3n = c(inds_GGATG_1n_3n,i) }
      if (length(grep('(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)',aseq))>0) {
        inds_GGATG_1n_4n = c(inds_GGATG_1n_4n,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)',aseq))>0) {
        inds_GGATG_2n_4n = c(inds_GGATG_2n_4n,i) }
      if (length(grep('G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)',aseq))>0) {
        inds_GGATG_3n_4n = c(inds_GGATG_3n_4n,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)',aseq))>0) {
        inds_GGATG_2n_4n = c(inds_GGATG_2n_4n,i) }
      if (length(grep('(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}',aseq))>0) {
        inds_GGATG_1n_2n_3n = c(inds_GGATG_1n_2n_3n,i) }
      if (length(grep('(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)',aseq))>0) {
        inds_GGATG_1n_2n_4n = c(inds_GGATG_1n_2n_4n,i) }
      if (length(grep('(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)',aseq))>0) {
        inds_GGATG_1n_3n_4n = c(inds_GGATG_1n_3n_4n,i) }
      if (length(grep('G{3,}[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)[ATGC]{1,7}(GG[ATC]{2,5}G|G[ATC]{2,5}GG)',aseq))>0) {
        inds_GGATG_2n_3n_4n = c(inds_GGATG_2n_3n_4n,i) }
      
      # C-rich strand
      if(length(grep('(C{2,}[ATGC]{1,7}){3,}C{2,}',aseq))>0) {
        inds_GG = c(inds_GG,i) }
      if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{1,7}C|C[ATG]{1,7}CC)',aseq))>0){
        inds_GGATG_1 = c(inds_GGATG_1,i) }
      if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{1,7}C|C[ATG]{1,7}CC)[ATGC]{1,7}C{3,}',aseq))>0){
        inds_GGATG_2 = c(inds_GGATG_2,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(CC[ATG]{1,7}C|C[ATG]{1,7}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}',aseq))>0){
        inds_GGATG_3 = c(inds_GGATG_3,i) }
      if (length(grep('(CC[ATG]{1,7}C|C[ATG]{1,7}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_4 = c(inds_GGATG_4,i) }
      if (length(grep('(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_1_2 = c(inds_GGATG_1_2,i) }
      if (length(grep('(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_1_3= c(inds_GGATG_1_3,i) }
      if (length(grep('(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)',aseq))>0) {
        inds_GGATG_1_4= c(inds_GGATG_1_4,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}',aseq))>0){
        inds_GGATG_2_3 = c(inds_GGATG_2_3 ,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)',aseq))>0){
        inds_GGATG_2_4 = c(inds_GGATG_2_4 ,i) }
      if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)',aseq))>0) {
        inds_GGATG_3_4 = c(inds_GGATG_3_4 ,i) }
      if (length(grep('(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_1_2_3 = c(inds_GGATG_1_2_3 ,i) }
      if (length(grep('(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)',aseq))>0) {
        inds_GGATG_1_2_4 = c(inds_GGATG_1_2_4 ,i) }
      if (length(grep('(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)',aseq))>0) {
        inds_GGATG_1_3_4 = c(inds_GGATG_1_3_4 ,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)[ATGC]{1,7}(CC[ATG]C|C[ATG]CC)',aseq))>0) {
        inds_GGATG_2_3_4 = c(inds_GGATG_2_3_4 ,i) }
      if (length(grep('(G[ATG]G[ATG]C)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_1_1 = c(inds_GGATG_1_1 ,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(G[ATG]G[ATG]C)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_2_2 = c(inds_GGATG_2_2,i) }
      if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(G[ATG]G[ATG]C)[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_3_3 = c(inds_GGATG_3_3,i) }
      if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(G[ATG]G[ATG]C)',aseq))>0) {
        inds_GGATG_4_4 = c(inds_GGATG_4_4,i) }
      if (length(grep('(G[ATG]{2,5}G[ATG]{2,5}C)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_1n_1n = c(inds_GGATG_1n_1n,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(G[ATG]{2,5}G[ATG]{2,5}C)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_2n_2n = c(inds_GGATG_2n_2n,i) }
      if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(G[ATG]{2,5}G[ATG]{2,5}C)[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_3n_3n = c(inds_GGATG_3n_3n,i) }
      if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(G[ATG]{2,5}G[ATG]{2,5}C)',aseq))>0) {
        inds_GGATG_4n_4n = c(inds_GGATG_4n_4n,i) }
      if (length(grep('(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_1n_2n = c(inds_GGATG_1n_2n,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_2n_3n = c(inds_GGATG_2n_3n,i) }
      if (length(grep('(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_1n_3n = c(inds_GGATG_1n_3n,i) }
      if (length(grep('(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)',aseq))>0) {
        inds_GGATG_1n_4n = c(inds_GGATG_1n_4n,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)',aseq))>0) {
        inds_GGATG_2n_4n = c(inds_GGATG_2n_4n,i) }
      if (length(grep('C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)',aseq))>0) {
        inds_GGATG_3n_4n = c(inds_GGATG_3n_4n,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)',aseq))>0) {
        inds_GGATG_2n_4n = c(inds_GGATG_2n_4n,i) }
      if (length(grep('(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}',aseq))>0) {
        inds_GGATG_1n_2n_3n = c(inds_GGATG_1n_2n_3n,i) }
      if (length(grep('(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)',aseq))>0) {
        inds_GGATG_1n_2n_4n = c(inds_GGATG_1n_2n_4n,i) }
      if (length(grep('(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)',aseq))>0) {
        inds_GGATG_1n_3n_4n = c(inds_GGATG_1n_3n_4n,i) }
      if (length(grep('C{3,}[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)[ATGC]{1,7}(CC[ATG]{2,5}C|C[ATG]{2,5}CC)',aseq))>0) {
        inds_GGATG_2n_3n_4n = c(inds_GGATG_2n_3n_4n,i) }
    }          
    
  }  
  
  count7 <- length(unique(inds_7))
  count5 <- length(unique(inds_5))
  count3 <- length(unique(inds_3))
  inds_long <- setdiff(c(inds_12,inds_21, inds_7, inds_5, inds_3), c(inds_7, inds_5, inds_3))
  count_long <- length(inds_long)
  inds_bulges <- (unique(c(inds_GGATG_1, inds_GGATG_2, inds_GGATG_3, inds_GGATG_4,
                           inds_GGATG_1_1,inds_GGATG_2_2,inds_GGATG_3_3,inds_GGATG_4_4,
                           inds_GGATG_1_2, inds_GGATG_1_3, inds_GGATG_1_4, inds_GGATG_2_3, inds_GGATG_2_4, inds_GGATG_3_4, 
                                            'long_loops',
                           inds_GGATG_1_2_3, inds_GGATG_1_2_4,inds_GGATG_1_3_4, inds_GGATG_2_3_4,
                           inds_GGATG_1n_1n, inds_GGATG_2n_2n, inds_GGATG_3n_3n, inds_GGATG_4n_4n,
                           inds_GGATG_1n_2n, inds_GGATG_1n_3n, inds_GGATG_1n_4n, inds_GGATG_2n_3n,
                           inds_GGATG_2n_4n, inds_GGATG_3n_4n, inds_GGATG_1n_2n_4n, inds_GGATG_1n_2n_3n,
                           inds_GGATG_1n_3n_4n, inds_GGATG_2n_3n_4n)))
  count_bulges <- length(inds_bulges)
  
  inds_simple_bulges <- unique(c(inds_GGATG_1,inds_GGATG_2,inds_GGATG_3,inds_GGATG_4,
                                 #inds_GGATG_1_2, inds_GGATG_1_3, inds_GGATG_1_4,inds_GGATG_2_3, inds_GGATG_2_4,inds_GGATG_3_4,
                                 inds_GGATG_1_1,inds_GGATG_2_2,inds_GGATG_3_3,inds_GGATG_4_4))
  count_simple_bulges <- length(inds_simple_bulges)
  
  inds_complex_bulges <- setdiff(inds_bulges,inds_simple_bulges)
  count_complex_bulges <- length(inds_complex_bulges)
  inds_GG_nobulge <- setdiff(inds_GG,inds_bulges)
  count_GG_nobulge <- length(inds_GG_nobulge)
  
  # counts of the actual case
  # plot: counts_3, counts_5 - counts_3, counts_7 - counts_5, counts_long, counts_sbulges, counts_cbulges_GG, counts_other
  if (s1 == 1){
    counts_T[[count_f]] <- countsT;counts_A[[count_f]] <- countsA;counts_C[[count_f]] <- countsC;counts_G[[count_f]] <- countsG
    counts_CG[[count_f]] <- countsCG;counts_AT[[count_f]] <- countsAT;counts_skewGC[[count_f]] <- skewGC;counts_skewAT[[count_f]] <- skewAT
    
    lens[count_f] <- length(seq)
    counts_7[count_f] = count7;
    counts_5[count_f] = count5;
    counts_3[count_f] = count3
    counts_long[count_f] = count_long
    counts_bulges[count_f] = count_bulges
    counts_sbulges[count_f] = count_simple_bulges
    counts_cbulges[count_f] = count_complex_bulges
    counts_GG[count_f] = count_GG_nobulge
    counts_cbulges_GG[count_f] = count_complex_bulges + count_GG_nobulge
    counts_pqs[count_f] = count7 + count_long + count_bulges + count_GG_nobulge
    counts_fraction[count_f] = counts_pqs[count_f] / length(seq)
    counts_other[count_f] <- length(seq) - counts_pqs[count_f]
    inds_others <- setdiff(1:length(seq), unique(c(inds_7,inds_long,inds_bulges,inds_GG_nobulge)))
    # save to bed file the coordinates of those peaks not having any match
    bed_to_export <- bed_intervals_file[inds_others,]
    write.table(bed_to_export,file=gsub('.fa','.motifs_others.bed',file),quote = F, sep = "\t", col.names = F, row.names = F)
    
  }
  # counts of the randomly shuffled case number 2
  if (s1 == 2){
    counts_r1_T[[count_f]] <- countsT;counts_r1_A[[count_f]] <- countsA;counts_r1_C[[count_f]] <- countsC;counts_r1_G[[count_f]] <- countsG
    counts_r1_CG[[count_f]] <- countsCG;counts_r1_AT[[count_f]] <- countsAT;counts_r1_skewGC[[count_f]] <- skewGC;counts_r1_skewAT[[count_f]] <- skewAT
    counts_r1_7[count_f] = count7;counts_r1_5[count_f] = count5;counts_r1_3[count_f] = count3
    counts_r1_long[count_f] = count_long
    counts_r1_bulges[count_f] = count_bulges
    counts_r1_sbulges[count_f] = count_simple_bulges
    counts_r1_cbulges[count_f] = count_complex_bulges
    counts_r1_GG[count_f] = count_GG_nobulge
    counts_r1_cbulges_GG[count_f] = count_complex_bulges + count_GG_nobulge
    counts_r1_pqs[count_f] = count7 + count_long + count_bulges + count_GG_nobulge
    counts_r1_fraction[count_f] = counts_r1_pqs[count_f] / length(seq)
    counts_r1_other[count_f] <- length(seq) - counts_r1_pqs[count_f]
  }
  # counts of the randomly shuffled casee number 2
  if (s1 == 3){
    counts_r2_T[[count_f]] <- countsT;counts_r2_A[[count_f]] <- countsA;counts_r2_C[[count_f]] <- countsC;counts_r2_G[[count_f]] <- countsG
    counts_r2_CG[[count_f]] <- countsCG;counts_r2_AT[[count_f]] <- countsAT;counts_r2_skewGC[[count_f]] <- skewGC;counts_r2_skewAT[[count_f]] <- skewAT
    counts_r2_7[count_f] = count7;counts_r2_5[count_f] = count5;counts_r2_3[count_f] = count3
    counts_r2_long[count_f] = count_long
    counts_r2_bulges[count_f] = count_bulges
    counts_r2_sbulges[count_f] = count_simple_bulges
    counts_r2_cbulges[count_f] = count_complex_bulges
    counts_r2_GG[count_f] = count_GG_nobulge
    counts_r2_cbulges_GG[count_f] = count_complex_bulges + count_GG_nobulge
    counts_r2_pqs[count_f] = count7 + count_long + count_bulges + count_GG_nobulge
    counts_r2_fraction[count_f] = counts_r2_pqs[count_f] / length(seq)
    counts_r2_other[count_f] <- length(seq) - counts_r2_pqs[count_f]
  }
  # counts of the randomly shuffled case number 3
  if (s1 == 4){
    counts_r3_T[[count_f]] <- countsT;counts_r3_A[[count_f]] <- countsA;counts_r3_C[[count_f]] <- countsC;counts_r3_G[[count_f]] <- countsG
    counts_r3_CG[[count_f]] <- countsCG;counts_r3_AT[[count_f]] <- countsAT;counts_r3_skewGC[[count_f]] <- skewGC;counts_r3_skewAT[[count_f]] <- skewAT
    counts_r3_7[count_f] = count7;counts_r3_5[count_f] = count5;counts_r3_3[count_f] = count3
    counts_r3_long[count_f] = count_long
    counts_r3_bulges[count_f] = count_bulges
    counts_r3_sbulges[count_f] = count_simple_bulges
    counts_r3_cbulges[count_f] = count_complex_bulges
    counts_r3_GG[count_f] = count_GG_nobulge
    counts_r3_cbulges_GG[count_f] = count_complex_bulges + count_GG_nobulge
    counts_r3_pqs[count_f] = count7 + count_long + count_bulges + count_GG_nobulge
    counts_r3_fraction[count_f] = counts_r3_pqs[count_f] / length(seq)
    counts_r3_other[count_f] <- length(seq) - counts_r3_pqs[count_f]
  }
}

#}
#}

# create data frame and other output results
# fold enrichment will be for instance counts_7 / (average(counts_r1_7,counts_r2_7,counts_r3_7))
#save.image('sequence_analysis_ChIP.RData')

current_actual_count <- c(counts_3,
                          (counts_5-counts_3),
                          (counts_7-counts_5),
                          counts_long,
                          counts_sbulges, 
                          counts_cbulges_GG, 
                          counts_other)
current_radom_average <- c(mean(counts_r1_3,counts_r2_3,counts_r3_3),
                           mean(counts_r1_5,counts_r2_5,counts_r3_5) - mean(counts_r1_3,counts_r2_3,counts_r3_3),
                           mean(counts_r1_7,counts_r2_7,counts_r3_7) - mean(counts_r1_5,counts_r2_5,counts_r3_5),
                           mean(counts_r1_long,counts_r2_long,counts_r3_long),
                           mean(counts_r1_sbulges,counts_r2_sbulges,counts_r3_sbulges),
                           mean(counts_r1_cbulges_GG,counts_r2_cbulges_GG,counts_r3_cbulges_GG),
                           mean(counts_r1_other,counts_r2_other,counts_r3_other))
current_fold_change_over_random <- current_actual_count/current_radom_average

# create a dataframe with the various parameters
results_motifs <- data.frame(actual_count = current_actual_count,
                             radom_average = current_radom_average,
                             fold_change_over_random = current_fold_change_over_random,
                             motif_name = c('loop1-3',
                                            'loop4-5',
                                            'loop6-7',
                                            'long_loops',
                                            'simple_bulges',
                                            '2tetrads_cbulges',
                                            'others'))
results_motifs$motif_name <- factor(results_motifs$motif_name, 
                                    levels=c('loop1-3',
                                             'loop4-5',
                                             'loop6-7',
                                             'long_loops',
                                             'simple_bulges',
                                             '2tetrads_cbulges',
                                             'others'))

#save data.frame in RData
file_out <- gsub('.fa','.results_motifs.Rdata',file)
save(results_motifs,file=file_out)
# use ggplot to plot
#bp <- ggplot(results_motifs, aes(x="",y=actual_count, fill = motif_name)) + geom_bar(width = 1, stat = "identity")

