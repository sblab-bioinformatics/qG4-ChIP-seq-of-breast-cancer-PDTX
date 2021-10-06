README.md
================
angela simeone

About this repository
=====================

This repository collects bash and R scripts used for the analysis of the PDTX models by quantitative G4-ChIP-seq (qG4-ChIP-seq). 22 PDTX models have been profiled by qG4-ChIP-seq in 4 technical replicates plus input (5 libraries each).<br />

In our study (PID, DOI:) we process data obteined from qG4-ChIP-seq. We identify regions with differential binding (up-regulated) in each of the 22 breast cancer model and we call these regions as ΔG4.<br /> We then explore and characterize all the 22 ΔG4 regions in respect to expression, to copy number alterations and to transcritpion factor binding profiles.

[*Experimental design*](./Experimental_design_PDTXs).

Listed below are the different steps followed to perform the processing and analysis of qBG4-ChIP-seq

Data
====

Data produced using the qG4-ChIP-seq have been deposited at [GSE152216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152216). The sample file describing the deposited data is [here](sdrf.tsv). Additional information about the qG4-ChIP-seq protocol [here](idf.tsv).

The repository contains fastq files and processed files (peaks files, bw track files).

Data analysis and code
======================


Listed below are the different steps followed to perform the the analysis of qBG4-ChIP-seq.

-   [**Sequencing (Fastq) processing: alignment, peak calling and consensus regions generation from peaks**](./Basic_Processing.md): adapter removal, alignment, species separation, duplicates removal, peak calling, downsample input files and repeat peak calling, extraction of multi2 regions (regions confirmed in 2 out of 4 replicates), merge old consensus (first round of experiments) to new consensus (second round of experiments).

-   [**Analysis of PDTX samples qBG4-ChIP-seq**](./Analysis_first_second_phase_PDTXs.md)<br /> Use human consensus to estimate human signal, use [*Drosophila ROI consensus regions*](./input_files/consensus_dm6_over_21_multi2_samples.bed) to estimate drosophila signal, collect all stats (libraries sizes, coverages for human and drosophila respectively) necessary to proceed to drosophila normalization.<br /> Intermediate files are generated and used in the customized R script to generate sample file and to extract normalization factors.

    -   [**R script to generate sample file and compute drosophila normalization factors**](./Generate_sample_file_PDTX_only.R)
    -   [**R script to perform drosophila normalization**](./ChIP_analysis_normalize_data_input_subtraction_PDTX.R)
    -   [**R script to perform recursive differential analysis and generate DG4 regions (DG4R)**](./ChIP_analysis.diff_bind_PDTX.R)

-   [**Characterization of differential regions ΔG4**](./Characterization_DG4R.md)<br /> Pairwise comparison of PDTX G4 regions (Jaccard index), fold-enrichment of PDTXs at 45 cancer drivers genes, G4-motifs analysis, prepare data to analyze associatio of ΔG4 to expression.

-   [**Association of ΔG4 to expression**](/Compare_expressionValues_DG4R_CG4R.R)<br /> R script to explore association between presence of ΔG4 in promoter and relative expression levels in each individual PDTX model. This analysis shows that The presence of ΔG4 is generally associated to higly expressed gene.

-   [**Association of ΔG4 to CNA**](./Association_DG4R_to_CNA.md)<br /> Explore fold-enrichment of ΔG4 at CNA regions. This analysis is done for each individual PDTX (i.e. indivudual PDTXs CNA have been identified first and then ΔG4 have been compared to them).
    -   [**R script to identify CNA given input library bam (5M subsampling)**](./copy_number_alteration_identification_script.R)<br />
        This script uses the R package "QDNAseq".
-   [**Association of ΔG4 to CNA (copy number alterations) and expression**](./Comparison_Expression_DG4R_CNA.md)
    -   [**R script to perform the comparison between G4 levels, expression and CNA**](./Analysis_promoter_G4levels_expression_CNA.R)
-   [**Association of ΔG4 to SNV**](./Association_DG4R_to_SNV.md)<br />

-   [**Association of ΔG4 to TF binding sites (after downloading data from ChIP-Atlas)**](./Association_DG4R_to_CNA.md)

    -   [**R script to produce correaltions matrix from Fold Enrichment of TF data**](./transcirption_factors_breast_revised_9March2020.R)

Data files
==========

Here is the list of data files used at various stages of the analysis.

Input files to generate the sample file

-   [*file containing total number of reads hg19*](./input_files/tot_num_reads_hg19_libraries.txt)
-   [*file containing total number of reads of dm6*](./input_files/dm6_stats_with_fnames.txt)

Input files to perform drosophila normalization

-   [*Peak counts file*](./input_files/temp_concat_hg19_libraries.txt)
-   [*Sample file*](./input_files/sample_file_PDTX_only.txt)
-   [*Drosophila normalization factors*](./input_files/PDTX_old_new_dm6_norm5_peak_recovery.txt)

Input files (consesus regions) to extract coverages:

-   [*Drosophila ROI consensus regions*](./input_files/consensus_dm6_over_21_multi2_samples.bed)
-   [*Human consensus regions*](./input_files/hg19_old_new_q005.all_peaks.25M.over99nt.sorted.bed)

Input files to perform differential analysis

-   [*Peak counts file after input subtraction-library size norm- drosophila normalization*](peak_counts.norm_filtered.tab)
-   [*Sample file*](./input_files/sample_file_PDTX_only.txt)

Input files to perform the analysis of the association between DG4R and expression

-   [*Coordinates of gene promoters hg19*](./input_files/hg19.gene_name.promoters.bed)
-   [*List of promoter with overlapping DG4R*](./input_files/hg19.gene_name.promoters.DG4r.bed)
-   [*List of promoter with overlapping CG4R*](./input_files/hg19.gene_name.promoters.CG4r.bed)
-   [*Expression values data table*](./input_files/ExpModelsData_all_plus_PARsamples.txt)

Input files to perform analysis of helicases

-   [*List of helicases*](./input_files/helicase_list_angela.csv)

Input files for analysis of the association between DG4R and CNA

-   [*genome file*](hg19_robert_github.size.genome)
-   [*occurances of overlaps between DG4R and CNA\_AMP in the actual and randomized cases*](./input_files/PDTX_vs_all_AMP.sites)
-   [*occurances of overlaps between DG4R and CNA\_GAIN in the actual and randomized cases*](./input_files/PDTX_vs_all_GAIN.sites)
-   [*occurances of overlaps between DG4R and CNA\_NEUT in the actual and randomized cases*](./input_files/PDTX_vs_all_NEUT.sites)
-   [*occurances of overlaps between DG4R and CNA\_HETD in the actual and randomized cases*](./input_files/PDTX_vs_all_HETD.sites)
-   [*occurances of overlaps between DG4R and CNA\_HOMD in the actual and randomized cases*](./input_files/PDTX_vs_all_HOMD.sites)

Processed files
---------------

-   [*FE of DG4R at TF binding regions from ChIP-atlas*](/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction/ChIP_atlas/Fold_enrichments_DG4R_at_TFbindingRegions.processed.csv)
-   [*Corr of TF-FE (DG4R at TF binding regions from ChIP-atlas)*](/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction/ChIP_atlas/Corr_Spearm_TF_Fold_enrichments_DG4R_at_TFbindingRegions.processed.csv)
-   [*Corr of DG4R (FE of DG4R at TF binding regions from ChIP-atlas)*](/Users/simeon01/Documents/PDTX/REVISIONS/PDTX_for_figures_generation/PDTX_only_dm6_inputSubtraction/ChIP_atlas/Corr_Spearm_PDTX_Fold_enrichments_DG4R_at_TFbindingRegions.processed.csv)

Contact us
----------

For any query contact us: [Angela](mailto:angela.simeone@cruk.cam.ac.uk) and [Robert](mailto:robert.haensel-hertsch@uni-koeln.de).
