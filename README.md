This repository contains data access and computational analyses for the methods developed in our Hänsel-Hertsch et al. 2019 manuscript.

## Data

All the raw sequencing data, qG4-ChIP-seq fastq, bigwig and regions have been deposited in the ArrayExpress database at EMBL-EBI under accession number [E-MTAB-7949](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7949). For clarity purposes, 1) a table [sdrf](sdrf.tsv) explaining the source of all files in "E-MTAB-7949" and 2) a text file [idf](idf.tsv) describing more details of qG4-ChIP-seq protocol.


## Code

- [Mapping, peak calling and peak processing](Mapping_peak.calling_and_peak.processing.txt): Basic analysis to get all BG4-antibody enriched G4 regions.

- [Somatic copy number alteration profiling](Somatic.copy.number.alteration.profiling.txt): For each of 16 PDTX samples the genomic Input data was used to estimate copy number abberations.

- [DG4Rs and CG4Rs.](DG4Rs.and.CG4Rs.txt): Normalisation factor extraction from spike-in/reference coverage variation; R-script to perform differential BG4 enrichment analysis; Extraction of DG4R and CG4R files

- [DG4R PDTX stratification](DG4R.PDTX.stratification.txt): Systematic DG4R overlap with other DG4Rs and similarity analysis.

- [DG4R/CG4R enrichment in CNA classifications relative to random permutation](/DG4R.CG4R.enrichment.in.CNA.classifications.relative.to.random.permutation.txt): DG4R or CG4R enrichments in copy number aberrations.

- [PDTX qG4-ChIP-seq peak annotation and enrichment analysis](PDTX.qG4-ChIP-seq.peak.annotation.and.enrichment.analysis.txt): Genome annotation of PDTX BG4 enriched regions.

- [Promoter - G4 - Gene expression analysis](Promoter.G4.Gene.expression.analysis.txt): Systematic DG4R overlap with other DG4Rs and similarity analysis.

- [Promoter - G4 intensity - Gene expression - CNA analysis](Promoter.G4.intensity.Gene.expression.CNA.analysis.txt): How populated are G4 structures at a given promtoter in relationship to CNA and gene expression status.  

- [DG4R/CG4R - Association - Upregulated genes determining integrative cluster signature](DG4R.CG4R.Association.Upregulated.genes.determining.integrative.cluster.signature.txt): DG4R association with the 10 different IC gene sets.

- [DG4R/CG4R - Association - 45 common driver regions](DG4R.CG4R.Association.45.common.driver.regions.txt): DG4R and CG4R association with 45 common breast cancer driver regions

- [Transcription factor binding site (TFBS) - DG4R enrichment analysis](Transcription.factor.binding.site.(TFBS).DG4R.enrichment.analysis.txt): DG4R enrichment in ChIP-ATLAS breast cancer transcription factor binding sites. 

















