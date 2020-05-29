#' @title identify CNA
#' @description identify level of copy number alteration, define AMP, GAIN, NEUT, and deletions
#' @param promoters list, expression levels, list of promoters overlapping at least one ΔG4r
#' @return copyNumbersCalled.bed that needs to be processed in bash
#' @author angela simeone  \email{angela.simeone@cruk.cam.ac.uk}



library(QDNAseq)

base_dir <- '/mnt/nfs/nas/group_folders/simeon01/Data/20191216_PDTX_bam_5M/'
setwd(base_dir)

w <- 100
bins <- getBinAnnotations(binSize=w)

#files <- c('139M_181_input.all.hg19.merged.nodup.5M.bam',
#           '139M_284_input.all.hg19.merged.nodup.5M.bam',
#           '143_284_input.all.hg19.merged.nodup.5M.bam',
#           '143_317_input.all.hg19.merged.nodup.5M.bam',
#           '179_181_input.all.hg19.merged.nodup.5M.bam',
#           '179_284_input.all.hg19.merged.nodup.5M.bam',
#           '201_181_input.all.hg19.merged.nodup.5M.bam',
#           '201_284_input.all.hg19.merged.nodup.5M.bam',
#           '316_284_input.all.hg19.merged.nodup.5M.bam',
#           '331_284_input.all.hg19.merged.nodup.5M.bam',
#           '521_284_input.all.hg19.merged.nodup.5M.bam',
#           '5_317_input.all.hg19.merged.nodup.5M.bam',
#           '863_317_input.all.hg19.merged.nodup.5M.bam',
#           '9_317_input.all.hg19.merged.nodup.5M.bam',
#           '98_181_input.all.hg19.merged.nodup.5M.bam',
#           '98_284_input.all.hg19.merged.nodup.5M.bam')

files <- list.files(pattern = 'bam')
for (i in 1:length(files))
{
  print(files[i])
  aname <- gsub('.hg19.*','',files[i])
  readCounts <- binReadCounts(bins, bamfiles=files[i])
  pdf(paste(aname,w,'read_counts.pdf',sep='_'),8,8)
  plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
  dev.off()
  
  readCountsFiltered <- readCounts
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE, chromosomes="Y")
  #pdf(paste(aname,w,'isobarPlot.pdf',sep='_'),8,8)
  #isobarPlot(readCountsFiltered)
  #dev.off()
  
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  pdf(paste(aname,w,'noisePlot.pdf',sep='_'),8,8)
  noisePlot(readCountsFiltered)
  dev.off()
  
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  pdf(paste(aname,w,'copyNumberSmooth.pdf',sep='_'),8,8)
  plot(copyNumbersSmooth)
  dev.off()
  
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="none")
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  pdf(paste(aname,w,'copyNumbersSegmented.pdf',sep='_'),8,8)
  plot(copyNumbersSegmented)
  dev.off()
  
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  pdf(paste(aname,w,'copyNumbersCalled.pdf',sep='_'),8,8)
  plot(copyNumbersCalled)
  dev.off()
  
  exportBins(copyNumbersCalled, file=paste(aname,w,'copyNumbersCalled.bed',sep='_'),format='bed')
  
}



## in bash - from command line generate regions defined as:
# 1)GAIN if(($5 > 0.25) && ($5 <= 0.75))
# 1)HETD if(($5 > -1.4) && ($5 < -0.3))
# 1)NEUT if(($5 >= -0.3) && ($5 <= 0.25))
# 1)HOMD if($5 <= -1.4)
# 1)AMP if($5 > 0.75)

for file in `ls *copyNumbersCalled.bed  | grep -v "STG197*"`  # 197 is the normal tussue (failed) that we have to exclude
do
echo $file
base_name=${file%%_input*}
echo $base_name
tail -n +1 $file | awk -F "\t" '{if(($5 > 0.25) && ($5 <= 0.75)) {print "chr"$1"\t"$2"\t"$3"\t"$5}}' | sortBed -i - |mergeBed -i - > $base_name.CNA_GAIN.100kb.merged.bed
tail -n +1 $file |awk -F "\t" '{if(($5 > -1.4) && ($5 < -0.3)) { print "chr"$1"\t"$2"\t"$3"\t"$5}}'|sortBed -i - |mergeBed -i -  > $base_name.CNA_HETD.100kb.merged.bed
tail -n +1 $file |awk -F "\t" '{if(($5 >= -0.3) && ($5 <= 0.25)) { print "chr"$1"\t"$2"\t"$3"\t"$5}}'|sortBed -i - |mergeBed -i -  > $base_name.CNA_NEUT.100kb.merged.bed
tail -n +2 $file |awk -F "\t" '{if($5 <= -1.4) {print "chr"$1"\t"$2"\t"$3"\t"$5}}'|sortBed -i - |mergeBed -i - > $base_name.CNA_HOMD.100kb.merged.bed
tail -n +1 $file |awk -F "\t" '{if($5 > 0.75) {print "chr"$1"\t"$2"\t"$3"\t"$5}}'|sortBed -i - |mergeBed -i - > $base_name.CNA_AMP.100kb.merged.bed
echo "======"
done