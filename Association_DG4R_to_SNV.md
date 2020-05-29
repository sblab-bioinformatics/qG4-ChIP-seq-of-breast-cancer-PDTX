Association of G4 signal to SNV
================
angela simeone

ΔG4r association to SNV
-----------------------

Here we assess the association of single nucleotide variants (SNV) studied in each PDTX to the ΔG4r.

The analysis involved different steps.

-   shuffle DG4r acorss genome
-   overlap DG4r with individual SNV the respective PDTXs (i.e. grep model name in the bed file with all SNV)
-   collect the overlaps of the actual data
-   collect the overlaps obtained using randomized data
-   print data to [file](./input_files/res_DG4R_SNV.tab)

``` bash

rm res
for f in p*up.bed
do

  echo $f
  f2=${f/peak*detable./}
  f3=${f2/.*/}
  echo $f3
  grep -e "\t$f3" ../SNSsSamples_into_bed.bed > tmp.bed
  shuffleBed -i $f -g hg19.genome -seed 10 > ${f%%.bed}.shuf_10.bed
  shuffleBed -i $f -g hg19.genome -seed 20 > ${f%%.bed}.shuf_20.bed
  shuffleBed -i $f -g hg19.genome -seed 30 > ${f%%.bed}.shuf_30.bed
  shuffleBed -i $f -g hg19.genome -seed 40 > ${f%%.bed}.shuf_40.bed
  shuffleBed -i $f -g hg19.genome -seed 50 > ${f%%.bed}.shuf_50.bed

  echo $f >> res
  wc -l tmp.bed >> res
  intersectBed -a $f -b tmp.bed -wa | wc -l >> res
  intersectBed -a ${f%%.bed}.shuf_10.bed -b tmp.bed -wa | wc -l >> res
  intersectBed -a ${f%%.bed}.shuf_20.bed -b tmp.bed -wa | wc -l >> res
  intersectBed -a ${f%%.bed}.shuf_30.bed -b tmp.bed -wa | wc -l >> res
  intersectBed -a ${f%%.bed}.shuf_40.bed -b tmp.bed -wa | wc -l >> res
  intersectBed -a ${f%%.bed}.shuf_50.bed -b tmp.bed -wa | wc -l >> res

done

cat res | paste - - - - - - - - > res_DG4R_SNV.tab
```
