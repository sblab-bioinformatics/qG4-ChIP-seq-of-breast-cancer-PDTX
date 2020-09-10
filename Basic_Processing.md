Processing data from fastq files to peaks and consensus regions (peaks)
================
angela simeone

Data
----

Download and store data in a `fastq` folder.

Processing of the data
----------------------

Steps:

-   trimming
-   alignment
-   spieces splitting
-   merge lanes
-   markduplicates
-   recover stats

``` bash

# === trimming ===
cd ~/fastq
mkdir trimmed
for fq in *.fastq.gz
do
  
  bname=`basename $fq`
  sbatch -o %j.out -e %j.err --mem 16000 --wrap "cutadapt -q 20 -O 3 -a CTGTCTCTTATACACATCT -o trimmed/${bname%%.fastq.gz}.trimmed.fastq.gz $fq"
done

# == align ==
g='hg19_mm10_dm6.fa'
cd ~/fastq/trimmed
mkdir aligned
for f in *.trimmed.fastq.gz
do
    sbatch --mem 16G --wrap "bwa mem -t 8 -M $g $f | samtools view -Sb -F2308 -q 10 > aligned/${f%%.trimmed.fastq.gz}.hg19_mm10_dm6.bam"
done

for f in *bam
do
  sbatch --mem 16G --wrap "samtools sort -@ 8 $f > ${f%%.bam}.sort.bam"
done

# === species splitting ===
cd ~/fastq/trimmed/aligned
hg19_genome=~/reference_genomes/PDTX_genome/hg19_mm10_dm6/hg19.genome.bed
mm10_genome=~/reference_genomes/PDTX_genome/hg19_mm10_dm6/mm10.genome.bed
dm6_genome=~/reference_genomes/PDTX_genome/hg19_mm10_dm6/dm6.genome.bed
for f in *bam
do
  #hg19
  sbatch --mem 8G -o $f.out -e $f.err -J view --wrap "samtools view -b -@8 $f -L $hg19_genome > ${f%%.hg19_mm10_dm6.sort.bam}.hg19.bam"
  #dm6
  sbatch --mem 8G -o $f.out -e $f.err -J view --wrap "samtools view -b -@8 $f -L $dm6_genome > ${f%%.hg19_mm10_dm6.sort.bam}.dm6.bam"
  #mm10
  sbatch --mem 8G -o $f.out -e $f.err -J view --wrap "samtools view -b -@8 $f -L $mm10_genome > ${f%%.hg19_mm10_dm6.sort.bam}.mm10.bam"
done

# === merge files of same condition and same sequencing run === 
cd ~/fastq/trimmed/aligned
for f in *_L001*.dm6.bam
do
  base_name=${f%%.SLX*}
  base_name=${f%%.SLx*}
  echo $base_name
  #reference_name
  #ls ${base_name}*hg19.bam
  
  hg19_f=`ls ${base_name}*hg19.bam`
  dm6_f=`ls ${base_name}*dm6.bam`
  mm10_f=`ls ${base_name}*mm10.bam`
  echo $hg19_f
  echo $dm6_f
  echo $mm10_f
  
  #hg19
  sbatch --mem 8G --wrap "samtools merge -@8 ./merged_lines/${hg19_f%%_L001*}.all.hg19.merged.bam ${hg19_f%%_L001*}*hg19.bam"
  #dm6
  sbatch --mem 8G --wrap "samtools merge -@8 ./merged_lines/${dm6_f%%_L001*}.all.dm6.merged.bam ${hg19_f%%_L001*}*dm6.bam"
  #mm10
  sbatch --mem 8G --wrap "samtools merge -@8 ./merged_lines/${mm10_f%%_L001*}.all.mm10.merged.bam ${mm10_f%%_L001*}*mm10.bam"
  
  echo " ======================================= "
  echo " "
done


cd ~/fastq/trimmed/aligned
mkdir merged_lines
mv *all*bam ./merged_lines
cd merged


# === markduplicates (picard) and recover stats on number of reads == 
cd ~/fastq/trimmed/aligned/merged_lines
for f in *.bam
do
  sbatch -o %j.$bam.tmp.out -e %j.$bam.tmp.err --mem 32000  --wrap "java -Xmx7g -jar ~/applications/picard-2.20.3.jar MarkDuplicates I=$f O=${f%%.bam}.nodup.bam M=${f%%.bam}.md.txt AS=true REMOVE_DUPLICATES=true"
  echo "==========="
  echo "          "
done

# === recover stats ===
cd /scratcha/sblab/simeon01/Data/190902_NEW_PDTX_NEW_normalBreast/SLX-18351/trimmed/aligned
touch all_dm6_counts
for f in *dm6*nodup.bam; do echo $f;echo $f >> all_dm6_counts &&  samtools view -c $f >> all_dm6_counts; done

touch all_mm10_counts
for f in *mm10*nodup.bam; do echo $f;echo $f >> all_mm10_counts &&  samtools view -c $f >> all_mm10_counts; done

touch all_hg19_counts
for f in *hg19*nodup.bam; do echo $f;echo $f >> all_hg19_counts &&  samtools view -c $f >> all_hg19_counts; done


touch all_counts
for f in *.bam ; do echo $f;echo $f >> all_counts &&  samtools view -c $f >> all_counts; done
for f in KZ45*input*bam
do echo $f
echo $f >> all_counts_KZ45 &&  samtools view -c $f >> all_counts_KZ45; done
awk '{printf "%s\t%s",$0,(NR%2?FS:RS)}' all_counts.txt > all_couns.tab.txt
awk '{printf "%s\t%s",$0,(NR%2?FS:RS)}' all_counts_KZ45 > all_counts_KZ45.txt
```

Fix naming problem with one of the samples
------------------------------------------

One run name was wrong. Fixt it : “KZ45\_AB551\_input” -- IS --&gt; “KZ45\_AB790\_input”.

``` bash
ls KZ45_AB551_input*

for file in `ls KZ45_AB551_input*`
do
  echo $file
  echo ${file/KZ45_AB551_input/KZ45_AB790_input}
  mv $file ${file/KZ45_AB551_input/KZ45_AB790_input}
  echo "    "
done
```

Peak calling
------------

Call peaks for each replicate using the paired input.

``` bash


Call peaks on not-subsampled files
```

``` bash
cd ~/fastq/trimmed/aligned/merged_lines


# == recover PDTX model using the PDTX names ===

out_dir_narrow='~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep'
mkdir $out_dir_narrow
pdtx_ids=(AB636 AB577 PAR1022 AB580 STG139 AB555 AB551 STG282 AB790 STG195 PAR1006)
rep_num=(1 2 3 4)

for curr_pdtx in ${pdtx_ids[@]}
do
  echo "curr_pdtx: $curr_pdtx"
  t_hg19=`ls *${curr_pdtx}*ChIP*.hg19.merged.nodup.bam`
  
  c_hg19=`ls *${curr_pdtx}*input*.hg19.merged.nodup.bam`
  c_dm6=`ls *${curr_pdtx}*input*.dm6.merged.nodup.bam`
    
  #ls -lht $t_hg19
  #ls -lht $t_dm6
  
  #ls -lth $c_hg19
  #ls -lth $c_dm6
  
  # additional_loop over 4 reps
  n=1
  for t_hg19_curr in ${t_hg19[@]}
  do
    echo $rep
    curr_t_hg19=$t_hg19_curr
    curr_t_dm6=${curr_t_hg19/hg19/dm6}
    
    echo "CHIP:"
    echo "individual hg19 ChIP: $curr_t_hg19"
    echo "individual dm6 ChIP: $curr_t_dm6"
    
    ls -l $curr_t_hg19
    ls -l $curr_t_dm6
    
    echo "inputs:"
    echo "individual hg19 INPUT: $c_hg19"
    echo "individual dm6 INPUT: $c_dm6"
    
    ls -l $c_hg19
    ls -l $c_dm6
    
    echo "rep number: $n"
    ((n++))
    echo "========="
    
  #actual peak calling step for human (default options)
  sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_hg19 -c $c_hg19 -n $out_dir_narrow/${curr_t_hg19%%.nodup.bam}.q005.all --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05 --keep-dup 'all'"
  sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_dm6 -c $c_dm6 -n $out_dir_narrow/${curr_t_dm6%%.nodup.bam}.q005.all --format=BAM --gsize 'dm' --bw=300 --qvalue 0.05 --keep-dup 'all'"
  
  #actual peak calling step for human (nomodel options)
  #sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_hg19 -c $c_hg19 -n $out_dir_narrow/${curr_t_hg19%%.nodup.bam}.q005.nomodel.all --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05 --keep-dup 'all' --nomodel"
  #sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_dm6 -c $c_dm6 -n $out_dir_narrow/${curr_t_dm6%%.nodup.bam}.q005.nomodel.all --format=BAM --gsize 'dm' --bw=300 --qvalue 0.05 --keep-dup 'all' --nomodel"
  echo "=========================="
  
  #actual peak calling step for human (nomodel option --pvalue 1e-2)
  #sbatch -o %j.${curr_t_hg19%%.nodup.bam}.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_hg19 -c $c_hg19 -n $out_dir_narrow/${curr_t_hg19%%.nodup.bam}.p1e-2.nomodel.all --format=BAM --gsize 'hs' --bw=300 --pvalue 1e-2 --keep-dup 'all' --nomodel"
  #sbatch -o %j.${curr_t_dm6%%.nodup.bam}.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_dm6 -c $c_dm6 -n $out_dir_narrow/${curr_t_dm6%%.nodup.bam}.p1e-2.nomodel.all --format=BAM --gsize 'dm' --bw=300 --pvalue 1e-2 --keep-dup 'all' --nomodel"
  echo "=========================="
  echo " "
  done
done
```

Generate track files
--------------------

Steps:

-   extract library size and print it to file;
-   generate tracks normalized by library size.

``` bash

# === print counts ===
cd ~/fastq/trimmed/aligned/merged_lines
for file in *.hg19.merged.nodup.bam
  do
  bam_hg38=$file
  
  bam_hg38_nodup=$file
  bam_dm6_nodup=${file%%.hg19.merged.nodup.bam}.dm6.merged.nodup.bam
  echo "====="
  ls -lth $bam_hg38_nodup
  ls -lth $bam_dm6_nodup
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${bam_hg38_nodup%%.bam}.stat5"
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_dm6_nodup >  ${bam_dm6_nodup%%.bam}.stat5"
done

# === index bam ===
for bam in *.hg19.merged.nodup.bam
    do
  sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam"
done

# === generate tracks ===
for file in *.hg19.merged.nodup.bam
  do
  bam_hg19_nodup=$file
  bam_dm6_nodup=${file%%.hg19.merged.nodup.bam}.dm6.merged.nodup.bam
  tot_r_hg19=`cat ${bam_hg19_nodup%%.bam}.stat5`
  tot_r_dm6=`cat ${bam_dm6_nodup%%.bam}.stat5`
  scal_factor_hg19=$(awk -v m=$tot_r_hg19 'BEGIN { print 1000000/m }')
  scal_factor_dm6=$(awk -v m=$tot_r_dm6 'BEGIN { print 1000000/m }')
  echo $scal_factor_hg19
  echo $scal_factor_dm6
  sbatch --mem 8G --wrap "bamCoverage --scaleFactor $scal_factor_dm6 -bs 10 -b $bam_dm6_nodup -of \"bigwig\" -o ${bam_dm6_nodup%%.bam}.w10.rpm.bw"
  sbatch --mem 8G --wrap "bamCoverage --scaleFactor $scal_factor_hg19 -bs 10 -b $bam_hg19_nodup -of \"bigwig\" -o ${bam_hg19_nodup%%.bam}.w10.rpm.bw"
done
```

Check stats on number of peaks
------------------------------

Check overlap of narow peaks with the respective organism OQS (obtained by merging K and Li, both stranded).

``` bash

# ===  total number of peaks ===
wc -l  *hg19*narrowPeak | awk '{print $2"\t"$1}'
wc -l **.hg19.merged.q005.all_peaks.narrowPeak | awk '{print $2"\t"$1}'
wc -l *.all.hg19.merged.q005.nomodel.all_peaks.narrowPeak| awk '{print $2"\t"$1}'
# === n peaks overlapping OQs (hg19) ===
oqs=/Users/simeon01/Documents/OQs/OQ_hits.merged.sorted.bed
oqs_dm6=/Users/simeon01/Documents/OQs/all_dm6.sorted.bed

#for file in *.hg19.merged.q005.all_peaks.narrowPeak
for file in *.all.hg19.merged.q005.nomodel.all_peaks.narrowPeak
  do 
  #echo $file
  intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l
  #echo "===="
  done
```

Generate consensus peaks for each PDTX
--------------------------------------

Consensus regions for each PDTX are those observed in 2 out of the 4 technical replicates. We use multiIntersect to extract them.

``` bash

# === generate multi2 for all cases ===
pdtx_ids=(AB636 AB577 PAR1022 AB580 STG139 AB555 AB551 STG282 AB790 STG195 PAR1006 STG197)
mkdir multi2_bed
mkdir multi3_bed
for case in ${pdtx_ids[@]}
do
  echo $case
  bed_files=`ls *$case*.hg19.merged.q005.all_peaks.narrowPeak`
  bed_files_q005_nomodel=`ls *$case*.all.hg19.merged.q005.nomodel.all_peaks.narrowPeak`
  echo $bed_files
  echo $bed_files_q005_nomodel
  multiIntersectBed -i $bed_files | awk '{if($4>=2) print $0}' | sortBed -i | mergeBed -i - > ./multi2_bed/$case.hg19.merged.q005.all_peaksmulti2.bed
  multiIntersectBed -i $bed_files | awk '{if($4>=3) print $0}' | sortBed -i | mergeBed -i - > ./multi3_bed/$case.hg19.merged.q005.all_peaks.multi3.bed
  
  multiIntersectBed -i $bed_files_q005_nomodel | awk '{if($4>=2) print $0}' | sortBed -i | mergeBed -i - > ./multi2_bed/$case.all.hg19.merged.q005.nomodel.all_peaks.multi2.bed
  multiIntersectBed -i $bed_files_q005_nomodel | awk '{if($4>=3) print $0}' | sortBed -i | mergeBed -i - > ./multi3_bed/$case.all.hg19.merged.q005.nomodel.all_peaks.multi3.bed
  
  echo " ================ . ================"
done
  
wc -l *multi2.bed | awk '{print $2"\t"$1}'

# === overlap multi2 hg19 with OQs ===
for file in  *multi2.bed
  do 
  intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l
  done

## == multi3 ==
wc -l *multi3.bed | awk '{print $2"\t"$1}'

# === overlap multi3 hg19 with OQs ===
for file in  *hg19*multi3.bed
  do 
  intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l
  done



# === n peaks (dm6) ===
wc -l *dm6*narrowPeak | awk '{print $2"\t"$1}' 
wc -l *dm6.merged.q005.all_peaks.narrowPeak| awk '{print $2"\t"$1}' 
wc -l *dm6.merged.q005.nomodel.all_peaks.narrowPeak| awk '{print $2"\t"$1}' 

# === n peaks overlapping OQs (dm6) ===
oqs_dm6=/Users/simeon01/Documents/OQs/all_dm6.sorted.bed
rob_dm6=/Users/simeon01/Documents/PDTX/PDTX_ROBERT/macs2_output/dm6.all_peaks.5M.bed
#for file in  *dm6.merged.q005.all_peaks.narrowPeak
for file in  *dm6.merged.q005.nomodel.all_peaks.narrowPeak
  do 
  intersectBed -a $file -b $oqs_dm6 -wa | sort | uniq | wc -l
  done
# === n peaks overlapping rob all peaks (dm6) ===
#for file in  *dm6.merged.q005.all_peaks.narrowPeak
for file in  *dm6.merged.q005.nomodel.all_peaks.narrowPeak
  do 
  intersectBed -a $file -b $rob_dm6 -wa | sort | uniq | wc -l
  done
  
# === generate multi2 dm6 for all cases ===
pdtx_ids=(AB636 AB577 PAR1022 AB580 STG139 AB555 AB551 STG282 AB790 STG195 PAR1006 STG197)
mkdir multi2_bed_dm6
mkdir multi3_bed_dm6
for case in ${pdtx_ids[@]}
do
  echo $case
  bed_files=`ls *$case*.dm6.merged.q005.all_peaks.narrowPeak`
  bed_files_q005_nomodel=`ls *$case*all.dm6.merged.q005.nomodel.all_peaks.narrowPeak`
  echo $bed_files
  echo $bed_files_q005_nomodel  
  multiIntersectBed -i $bed_files | awk '{if($4>=2) print $0}' | sortBed -i | mergeBed -i - > ./multi2_bed_dm6/$case.dm6.merged.q005.all_peaks.multi2.bed
  multiIntersectBed -i $bed_files_q005_nomodel | awk '{if($4>=2) print $0}' | sortBed -i | mergeBed -i - > ./multi2_bed_dm6/$case.all.dm6.merged.q005.nomodel.all_peaks.multi2.bed
  
  multiIntersectBed -i $bed_files | awk '{if($4>=3) print $0}' | sortBed -i | mergeBed -i - > ./multi3_bed_dm6/$case.dm6.merged.q005.all_peaks.multi3.bed
  multiIntersectBed -i $bed_files_q005_nomodel | awk '{if($4>=3) print $0}' | sortBed -i | mergeBed -i - > ./multi3_bed_dm6/$case.all.dm6.merged.q005.nomodel.all_peaks.multi3.bed
  echo " ================ . ================"
  done

# === n peaks multi2  (dm6) ===
wc -l *multi2.bed | awk '{print $2"\t"$1}' 
wc -l *dm6*multi3.bed | awk '{print $2"\t"$1}' 

# === n peaks multi2 overlapping OQs (dm6) ===
oqs_dm6=/Users/simeon01/Documents/OQs/all_dm6.sorted.bed
rob_dm6=/Users/simeon01/Documents/PDTX/PDTX_ROBERT/macs2_output/dm6.all_peaks.5M.bed
for file in  *multi2.bed
  do 
  intersectBed -a $file -b $oqs_dm6 -wa | sort | uniq | wc -l
  done
  for file in  *multi3.bed
  do 
  intersectBed -a $file -b $oqs_dm6 -wa | sort | uniq | wc -l
  done
# === n peaks overlapping rob all peaks (dm6) ===
for file in  *multi2.bed
  do 
  intersectBed -a $file -b $rob_dm6 -wa | sort | uniq | wc -l
  done
  for file in  *multi3.bed
  do 
  intersectBed -a $file -b $rob_dm6 -wa | sort | uniq | wc -l
  done


# comparison of new peaks with old rob peaks
peak143_284=/Users/simeon01/Documents/PDTX/PDTX_ROBERT/macs2_output/143_284_R1.all.hg19.q005.all_peaks.multi2.bed
peak143_317=/Users/simeon01/Documents/PDTX/PDTX_ROBERT/macs2_output/143_317_R1.hg19.q005.all_peaks.multi2.bed
dg4r_143=/Users/simeon01/Documents/PDTX/DBG4_PDTX/STG143_specific_G4.bed
for file in *hg19*multi2.bed
do
  #intersectBed -a $file -b $peak143_284 -wa | sort | uniq | wc -l
  #intersectBed -a $file -b $peak143_317 -wa | sort | uniq | wc -l
  intersectBed -a $file -b $dg4r_143 -wa | sort | uniq | wc -l
done


for file in *hg19*.narrowPeak
do
  #intersectBed -a $file -b $peak143_284 -wa | sort | uniq | wc -l
  #intersectBed -a $file -b $peak143_317 -wa | sort | uniq | wc -l
  intersectBed -a $file -b $dg4r_143 -wa | sort | uniq | wc -l
done
```

Downsample input bam files
--------------------------

Input bam files get downsampled to 25M reads. Definitive peaks have ben identified using inputs of 25M reads depth. As before, the steps are:

-   call peaks
-   generate multi2 (regions confirmed in 2 of the 4 peak replicates)
-   generate consesus of all new PDTX (merging all peaks)

Note: input libraries reported in the pubblication () have been subsampled to 25M (human) and 5M (drosophila). subsupling has been performed specificallyfor the following samples:
AB521M
AB863M
HCI005
HCI009
STG139M_181
STG139M_284
STG143_284
STG143_317
STG201_181
STG201_284
STG316
STG331
VHIO098_181
VHIO098_284
VHIO179_181
VHIO179_284


``` bash

### === DOWNSAMPLE ===

## ==  25M for human == 
cd ~/fastq/trimmed/aligned/merged_lines
for file in *input.all.hg19.merged.nodup.stat5
do
  stat5_file=`cat $file`
  echo $stat5_file
  if [ "$stat5_file" -gt 25000000 ]
  then
      echo "biggern than25M" $file
      factor_down=$(awk -v m=$stat5_file 'BEGIN { print 25000000/m }')
      echo $factor_down
      echo "***"
      sbatch --mem 8G --wrap "samtools view -@ 8 -s  $factor_down -b ${file%%.stat5}.bam > ${file%%.stat5}.25M.bam"
  else
    echo " =================> LESS"
    sbatch --mem 8G --wrap "cp ${file%%.stat5}.bam ${file%%.stat5}.less25M.bam"
  fi
  echo "================="
done

## ==  5M for drosophila == 
for file in *input.all.dm6.merged.nodup.stat5
do
  stat5_file=`cat $file`
  echo $stat5_file
  if [ "$stat5_file" -gt 5000000 ]
  then
      echo "bigger than 5M" $file
      factor_down=$(awk -v m=$stat5_file 'BEGIN { print 5000000/m }')
      echo $factor_down
      echo "***"
      sbatch --mem 8G --wrap "samtools view -@ 8 -s  $factor_down -b ${file%%.stat5}.bam > ${file%%.stat5}.5M.bam"
  else
    echo " =================> LESS"
    sbatch --mem 8G --wrap "cp ${file%%.stat5}.bam ${file%%.stat5}.less5M.bam"
  fi
  echo "================="
done
```

Peak calling using down-sampled input bams
------------------------------------------

Peaks are called with macs2 using default option.

``` bash
## ==  peak calling using downsamples input bam (both hg19 and dm6) == 
out_dir_narrow='~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs'
mkdir $out_dir_narrow
pdtx_ids=(AB636 AB577 PAR1022 AB580 STG139 AB555 AB551 STG282 AB790 STG195 PAR1006 STG197)
rep_num=(1 2 3 4)

for curr_pdtx in ${pdtx_ids[@]}
do
  echo "curr_pdtx: $curr_pdtx"
  t_hg19=`ls *${curr_pdtx}*ChIP*.hg19.merged.nodup.bam`
  
  c_hg19=`ls *${curr_pdtx}*input*.hg19.merged.nodup*25M.bam`
  c_dm6=`ls *${curr_pdtx}*input*.dm6.merged.nodup*5M.bam`
    
  #ls -lht $t_hg19
  #ls -lht $t_dm6
  
  #ls -lth $c_hg19
  #ls -lth $c_dm6
  
  # additional_loop over 4 reps
  n=1
  for t_hg19_curr in ${t_hg19[@]}
  do
    echo $rep
    curr_t_hg19=$t_hg19_curr
    curr_t_dm6=${curr_t_hg19/hg19/dm6}
    
    echo "CHIP:"
    echo "individual hg19 ChIP: $curr_t_hg19"
    echo "individual dm6 ChIP: $curr_t_dm6"
    
    ls -l $curr_t_hg19
    ls -l $curr_t_dm6
    
    echo "inputs:"
    echo "individual hg19 INPUT: $c_hg19"
    echo "individual dm6 INPUT: $c_dm6"
    
    ls -l $c_hg19
    ls -l $c_dm6
    
    echo "rep number: $n"
    ((n++))
    echo "========="
    
  #actual peak calling step for human (default options)
  sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_hg19 -c $c_hg19 -n $out_dir_narrow/${curr_t_hg19%%.nodup.bam}.q005.all --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05 --keep-dup 'all'"
  sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_dm6 -c $c_dm6 -n $out_dir_narrow/${curr_t_dm6%%.nodup.bam}.q005.all --format=BAM --gsize 'dm' --bw=300 --qvalue 0.05 --keep-dup 'all'"
  
  #actual peak calling step for human (nomodel options)
  sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_hg19 -c $c_hg19 -n $out_dir_narrow/${curr_t_hg19%%.nodup.bam}.q005.nomodel.all --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05 --keep-dup 'all' --nomodel"
  sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_dm6 -c $c_dm6 -n $out_dir_narrow/${curr_t_dm6%%.nodup.bam}.q005.nomodel.all --format=BAM --gsize 'dm' --bw=300 --qvalue 0.05 --keep-dup 'all' --nomodel"
  echo "=========================="
  
  #actual peak calling step for human (nomodel option --pvalue 1e-2)
  #sbatch -o %j.${curr_t_hg19%%.nodup.bam}.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_hg19 -c $c_hg19 -n $out_dir_narrow/${curr_t_hg19%%.nodup.bam}.p1e-2.nomodel.all --format=BAM --gsize 'hs' --bw=300 --pvalue 1e-2 --keep-dup 'all' --nomodel"
  #sbatch -o %j.${curr_t_dm6%%.nodup.bam}.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $curr_t_dm6 -c $c_dm6 -n $out_dir_narrow/${curr_t_dm6%%.nodup.bam}.p1e-2.nomodel.all --format=BAM --gsize 'dm' --bw=300 --pvalue 1e-2 --keep-dup 'all' --nomodel"
  echo "=========================="
  echo " "
  done
done
```

Extract multi2 for each PDTX and create consensus by organism
-------------------------------------------------------------

Confimed peaks (we define them as `multi2`) are generated for each PDTX by selecting regions observed in 2 out of 4 replicates. For each specie (human and drosophila), the consensus regions is defined by merging all multi2 PDTX files.

``` bash

## ===  extract multi2 for the different cases === 
# == hg19 cases == 
# hg19.merged.q005.all_peaks.narrowPeak 
# hg19.merged.q005.nomodel.all_peaks.narrowPeak
pdtx_ids=(AB636 AB577 PAR1022 AB580 STG139 AB555 AB551 STG282 AB790 STG195 PAR1006 STG197)
peaks_folder=(~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs ~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep)
cd ~/fastq/trimmed/aligned/merged_lines
for dir in ${peaks_folder[@]}
do
cd $dir
  for curr_pdtx in ${pdtx_ids[@]}
  do
    narrow_peaks_default=`ls *${curr_pdtx}*hg19.merged.q005.all_peaks.narrowPeak`
    narrow_peaks_default_nomodel=`ls *${curr_pdtx}*hg19.merged.q005.nomodel.all_peaks.narrowPeak`
    multiIntersectBed -i ${narrow_peaks_default} | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ${curr_pdtx}.hg19.q005.all_peaks.multi2.bed
    multiIntersectBed -i ${narrow_peaks_default_nomodel} | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ${curr_pdtx}.hg19.q005.nomodel.all_peaks.multi2.bed
  done
  cd ~/fastq/trimmed/aligned/merged_lines
  echo "       ===============                 "
done

## == hg19 consesus bed (merging all multi2 peaks) == 
# peaks on downsampled input bams
cd ~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs
multi2_peaks_default=`ls *hg19.q005.all_peaks.multi2.bed`
multi2_peaks_default_nomodel=`ls *hg19.q005.nomodel.all_peaks.multi2.bed`
cat $multi2_peaks_default | sortBed -i - | mergeBed -i - > hg19.q005.all_peaks.25M.bed
cat $multi2_peaks_default_nomodel | sortBed -i - | mergeBed -i - > hg19.q005.nomodel.all_peaks.25M.bed
for file in hg19*bed
do
  sbatch --mem 4G --wrap "awk '{if ((\$3-\$2) >= 100) print \$0}' $file | sortBed -i - > ${file%%.bed}.over99nt.sorted.bed"
done

# peaks on full input bams
cd ~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep
multi2_peaks_default=`ls *hg19.q005.all_peaks.multi2.bed`
multi2_peaks_default_nomodel=`ls *hg19.q005.nomodel.all_peaks.multi2.bed`
cat $multi2_peaks_default | sortBed -i - | mergeBed -i - > hg19.q005.all_peaks.bed
cat $multi2_peaks_default_nomodel | sortBed -i - | mergeBed -i - > hg19.q005.nomodel.all_peaks.bed
for file in hg19*bed
do
  sbatch --mem 4G --wrap "awk '{if ((\$3-\$2) >= 100) print \$0}' $file | sortBed -i - > ${file%%.bed}.over99nt.sorted.bed"
done

## ===  extract multi2 for the different cases === 
# == dm6 cases == 
# dm6.merged.q005.all_peaks.narrowPeak 
# dm6.merged.q005.nomodel.all_peaks.narrowPeak
pdtx_ids=(AB636 AB577 PAR1022 AB580 STG139 AB555 AB551 STG282 AB790 STG195 PAR1006 STG197)
peaks_folder=(~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs ~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep)
cd ~/fastq/trimmed/aligned/merged_lines
for dir in ${peaks_folder[@]}
do
cd $dir
for curr_pdtx in ${pdtx_ids[@]}
do
narrow_peaks_default=`ls *${curr_pdtx}*dm6.merged.q005.all_peaks.narrowPeak`
narrow_peaks_default_nomodel=`ls *${curr_pdtx}*dm6.merged.q005.nomodel.all_peaks.narrowPeak`
multiIntersectBed -i ${narrow_peaks_default} | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ${curr_pdtx}.dm6.q005.all_peaks.multi2.bed
multiIntersectBed -i ${narrow_peaks_default_nomodel} | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ${curr_pdtx}.dm6.q005.nomodel.all_peaks.multi2.bed
done
cd ~/fastq/trimmed/aligned/merged_lines
echo "       ===============                 "
done

## == dm6 consesus bed (merging all multi2 peaks) == 
# dm6 peaks on downsampled input bams
cd ~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs
multi2_peaks_default=`ls *dm6.q005.all_peaks.multi2.bed`
multi2_peaks_default_nomodel=`ls *dm6.q005.nomodel.all_peaks.multi2.bed`
cat $multi2_peaks_default | sortBed -i - | mergeBed -i - > dm6.q005.all_peaks.5M.bed
cat $multi2_peaks_default_nomodel | sortBed -i - | mergeBed -i - > dm6.q005.nomodel.all_peaks.5M.bed
for file in dm6*bed
do
sbatch --mem 4G --wrap "sortBed -i $file > ${file%%.bed}.sorted.bed"
done

# dm6 peaks on full input bams
cd ~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep
multi2_peaks_default=`ls *dm6.q005.all_peaks.multi2.bed`
multi2_peaks_default_nomodel=`ls *dm6.q005.nomodel.all_peaks.multi2.bed`
cat $multi2_peaks_default | sortBed -i - | mergeBed -i - > dm6.q005.all_peaks.bed
cat $multi2_peaks_default_nomodel | sortBed -i - | mergeBed -i - > dm6.q005.nomodel.all_peaks.bed
for file in dm6*bed
do
sbatch --mem 4G --wrap "sortBed -i $file > ${file%%.bed}.sorted.bed"
done
```

merge old Rob's consesus with new consesus by organism
------------------------------------------------------

Merge general consensus by merging consensus from first phase of the study and consensus from second phase of the study.

``` bash


### == unify old rob merged hg19 peaks with the new ones
mkdir ~/fastq/trimmed/aligned/merged_lines/consensus_hg19
cd ~/fastq/trimmed/aligned/merged_lines/consensus_hg19
hg19_old_all_merged_peaks=/scratchb/sblab/simeon01/RobertPDTX/multi2_hg19_25M/hg19.all_peaks.25M.over99nt.sort.bed
hg19_new_all_merged_peaks_q005=~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs/hg19.q005.all_peaks.25M.over99nt.sorted.bed
cat $hg19_old_all_merged_peaks $hg19_new_all_merged_peaks_q005 | sortBed -i - | mergeBed -i - | awk '{if (($3-$2)>=100) print $0}'| sortBed -i - | mergeBed -i - > hg19_old_new_q005.all_peaks.25M.over99nt.sorted.bed

### == unify old rob merged dm6 peaks with the new ones

mkdir ~/fastq/trimmed/aligned/merged_lines/consensus_dm6
cd ~/fastq/trimmed/aligned/merged_lines/consensus_dm6
dm6_old_all_merged_peaks=/scratchb/sblab/simeon01/RobertPDTX/multi2_dm6_5M/dm6.all_peaks.5M.bed
dm6_new_all_merged_peaks_q005=~/fastq/trimmed/aligned/merged_lines/macs2_individual_rep_downsampled_inputs/dm6.q005.all_peaks.5M.sorted.bed
cat $dm6_old_all_merged_peaks $dm6_new_all_merged_peaks_q005 | sortBed -i - | mergeBed -i - > dm6_old_new_q005.all_peaks.5M.sorted.bed
```
