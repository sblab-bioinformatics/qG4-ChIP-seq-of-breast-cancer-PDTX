PDTX processing after peak calling
================
angela simeone

We integrate data from the first and the second phase of the qG4-ChIP-seq study on PDTXs. The total number of models is 22.

In order to integrate everything we followed those steps:

1.  create consensus by merging multi2 of the selected cases (PDTX only);
2.  compute coverages at consensus;
3.  extract normalization factors from drosophila;
4.  normalize human signal;
5.  perform differential a customized differential analysis.

Create consensus regions for drosophila and human
-------------------------------------------------

For drosophila, we have created a general [*drosophila ROI consensus regions*](./input_files/consensus_dm6_over_21_multi2_samples.bed) set based on more than 100 samples. This is includes about 1300 regions across the drosophila genome

For human, the reference set is obtained by merging all PDTX multi2 regions (i.e. regions confirmed in 2 out of 4 technical replicats).

``` bash
### == merge multi2 regions ==

### == hg19 consensus ==
mkdir ~/fastq/trimmed/aligned/merged_lines/consensus_hg19
cd ~/fastq/trimmed/aligned/merged_lines/consensus_hg19
hg19_old_all_merged_peaks=~/RobertPDTX/multi2_hg19_25M/hg19.all_peaks.25M.over99nt.sort.bed
hg19_new_all_merged_peaks_q005=~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs/hg19.q005.all_peaks.25M.over99nt.sorted.bed
cat $hg19_old_all_merged_peaks $hg19_new_all_merged_peaks_q005 | sortBed -i - | mergeBed -i - | awk '{if (($3-$2)>=100) print $0}'| sortBed -i - | mergeBed -i - > hg19_old_new_q005.all_peaks.25M.over99nt.sorted.bed

### == dm6 consensus ==
mkdir ~/fastq/trimmed/aligned/merged_lines/consensus_dm6
cd ~/fastq/trimmed/aligned/merged_lines/consensus_dm6
dm6_old_all_merged_peaks=~/RobertPDTX/multi2_dm6_5M/dm6.all_peaks.5M.bed
dm6_new_all_merged_peaks_q005=~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs/dm6.q005.all_peaks.5M.sorted.bed
cat $dm6_old_all_merged_peaks $dm6_new_all_merged_peaks_q005 | sortBed -i - | mergeBed -i - > dm6_old_new_q005.all_peaks.5M.sorted.bed

## ==  dm6 consensus ROI == 
# ~/consensus_dm6_over_21_multi2_samples.bed

## == define consensus for coverages ==
consensus_PDTX_hg=~/fastq/trimmed/aligned/merged_lines/consensus_hg19/hg19_old_new_q005.all_peaks.25M.over99nt.sorted.bed
consensus_PDTX_dm=/scratchb/sblab/simeon01/all_dm6_multi2_peaks/consensus_dm6_21/consensus_dm6_over_21_multi2_samples.bed

path_newConsensus_hg19=~/consensus_hg19/coverage_hg19_allbam
path_newConsensus_dm6=~/consensus_dm6/coverage_dm6_allbam

mkdir $path_newConsensus_hg19
mkdir $path_newConsensus_dm6

cp $consensus_PDTX_hg ~/consensus_hg19
cp $consensus_PDTX_dm ~/consensus_dm6

## == bams dir ==
path_hg19=~/all_bams
path_dm6=~/all_bams

## == human signal (coverage) at human consensus ==
bam_paths_hg19=($path_hg19)
out_dir_hg19=$path_newConsensus_hg19
for curr_dir in ${bam_paths_hg19[@]}
do
  cd $curr_dir
  #for file in `ls *all.hg19.merged.nodup.bam | grep -v input` #all.hg19.merged.nodup.bam 
  #for file in `ls *all.hg19.merged.nodup.bam | grep input` #all.hg19.merged.nodup.bam 
  for file in *hg19*.bam
  do
    echo $file
    echo $consensus_PDTX_hg
    echo sbatch --mem 8G -o %j.log.out -e %j.log.err --wrap "bedtools coverage -a $consensus_PDTX_hg -b $file  -counts > $out_dir_hg19/${file%%.bam}.hg19.q005.25M.bedgraph"
    sbatch --mem 8G -o %j.log.out -e %j.log.err --wrap "bedtools coverage -a $consensus_PDTX_hg -b $file  -counts > $out_dir_hg19/${file%%.bam}.hg19.q005.25M.bedgraph"
  done
  echo "  =============== "
done


## == drosophila signal (coverage) at drosophila consensus ==
bam_paths_dm6=($path_dm6)
out_dir_dm6=$path_newConsensus_dm6
for curr_dir in ${bam_paths_dm6[@]}
do
  cd $curr_dir
  for file in *dm6*.bam
  do
    echo $file
    echo $consensus_PDTX_dm
    echo sbatch --mem 8G -o %j.log.out -e %j.log.err --wrap "bedtools coverage -a $consensus_PDTX_dm -b $file  -counts > $out_dir_dm6/${file%%.bam}.dm6.q005.5M.bedgraph"
    sbatch --mem 8G -o %j.log.out -e %j.log.err --wrap "bedtools coverage -a $consensus_PDTX_dm -b $file  -counts > $out_dir_dm6/${file%%.bam}.dm6.q005.5M.bedgraph"
  done
  echo "  =============== "
done
```

``` bash

path_stat5=~/all_bams
coverage_hg19=~/consensus_hg19/coverage_hg19_allbam
coverage_dm6=~/consensus_dm6/coverage_dm6_allbam

#rm dm6_stats_with_fnames.txt tot_num_reads_hg19_libraries.txt temp_concat_hg19_libraries.txt 

# cd ~/consensus_dm6/coverage_dm6_allbam # << === this is when all libraries merged peaks are used as DM6 consensus

cd $coverage_dm6

# Generated tRIP (total ReadsInPeaks, .tRIP.txt) files for drosophila

for file in *bedgraph
do
  awk '{sum += $4} END {print sum}' $file > ${file%%.bedgraph}.tRIP.txt
done


touch dm6_stats_with_fnames.txt
touch temp_concat_hg19_libraries.txt
touch tot_num_reads_hg19_libraries.txt

f1=`ls $coverage_hg19/*bedgraph | head -n 1`
cut -f 1,2,3 $f1  > temp_concat_hg19_libraries.txt


for file in *bedgraph
do
  pattern=${file%%dm6.q005*.bedgraph}
  stat5_file=`ls $path_stat5/$pattern*stat5`
  
  if [ -z "$stat5_file" ]
  then
    pattern2=${file%%input*}input*dm6*stat5
    echo $pattern2
    stat5_file=`ls $path_stat5/$pattern2`
    echo "***"
  fi
  echo $stat5_file
  cat $stat5_file
  echo "===== ===== ===== ====="
  
  file_dm6_reads_in_peaks=${file/.bedgraph}.tRIP.txt
  # print 3 col file
  if [[ $file_dm6_reads_in_peaks != *"input"* ]]; then
    tot_reads=`cat $stat5_file`
    rip=`cat $file_dm6_reads_in_peaks`
    norm_factor=`echo "scale=5; $rip / $tot_reads" | bc -l`
    echo $file$'\t'$tot_reads$'\t'$rip>> dm6_stats_with_fnames.txt
    
  else
    tot_reads=`cat $stat5_file`
    rip=`cat $file_dm6_reads_in_peaks`
    echo $file$'\t'$tot_reads$'\t'$rip>> dm6_stats_with_fnames.txt
    
  fi
  
  #check for each file the pattern in the human set
  #pattern for human bedgraph
  semi_pattern_hg19=${file/dm6/hg19}
  bedgraph_hg19=`ls $coverage_hg19/${semi_pattern_hg19%%hg19*}hg19*`
  stat_hg19=`ls $path_stat5/${semi_pattern_hg19%%hg19*}hg19*stat5`
  
  #coverage files
  cut -f 4 $bedgraph_hg19 | paste temp_concat_hg19_libraries.txt - > intermediate_hg19_libraries.txt
  mv intermediate_hg19_libraries.txt temp_concat_hg19_libraries.txt;
  
  # files with total number of reads hg19 and filename
  tot_reads_hg19=`cat $stat_hg19`
  echo ${semi_pattern_hg19%%hg19*}hg19$'\t'$tot_reads_hg19 >> tot_num_reads_hg19_libraries.txt
done


# copy locally the files of interest in ordert to run drosophila normalization and differential analysis
# files of interest:
# - tot_num_reads_hg19_libraries.txt
# - temp_concat_hg19_libraries.txt 
# - dm6_stats_with_fnames.txt
```

original: Analysis\_PDTX\_NO\_cells\_new\_old\_PDTX.\*
