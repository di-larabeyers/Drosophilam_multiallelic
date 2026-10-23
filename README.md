# Drosophilam_multiallelic

# Purpose

This workflow investigates the evolutionary dynamics of multi-allelic loci in pooled *Drosophila melanogaster* populations across multiple sampling time points. By combining temporal allele-frequency trajectories with fitness matrix estimation, it characterizes the fitness effects of multi-allelic genotypes and evaluates potential heterozygote advantage using both fitness-based metrics and the Lewontin condition. In parallel, multi-allelic sites are annotated to genomic features and genes to investigate their genomic distribution and functional context. Coverage-based analyses assess the relationship between sequencing depth and multi-allelic variation, while a targeted analysis of insecticide-resistance genes provides a more detailed examination of fitness patterns and multi-allelic variation in genes potentially involved in adaptive resistance. Together, these analyses integrate population-genomic, evolutionary, and functional perspectives to investigate the potential role of multiallelic variation in adaptation.

# Main workflow 
1. Locate the input BAM files and write a file list.
2. Fix read groups in BAM files.
3. Run allele extraction with FreeBayes, optionally applying a coverage filter.
4. Convert extracted allele counts to allele-frequency/trajectory tables using either hard filters or coverage filters.
5. Run allele statistics and downstream branches:
   - fitness matrix estimation and heterozygote-advantage analysis;
   - genomic annotation and gene-level summaries;
   - resistance-gene-specific extraction;
   - coverage summaries for multi-allelic sites.
  
# Data
- Pooled sequencing BAM files for the population/time points.
- D. melanogaster reference genome and its indexed reference files e.g. dmel-all-chromosome-r6.12.fasta from Flybase Release 6.12.
- BAM file list (list.txt).
- Gene annotation, e.g. genes_annotations_612.gff3 / FlyBase Release 6.12 annotation.
- Population identifier and estimated generation numbers for the sampling time points.

# Software
- FreeBayes v1.3.2-dirty
- SAMtools 1.13
- BCFtools 1.6
- BEDTools v2.27.1
- BGZIP 1.13
- Tabix 1.13
- MATLAB R2026a
- MATLAB, including Optimization Toolbox and Global Optimization Toolbox for the fitness-estimation branch

# Workflow and data processing

## Input preparation and read-group correction

The pipeline starts with the preparation of pooled sequencing data from *Drosophila melanogaster*. The BAM files are collected in a file list that serves as input for subsequent processing steps. Before variant extraction, the read groups of the BAM files are corrected using `fix_read_groups.sh`. This step ensures that the sequencing files contain the required read-group information for downstream analysis. The corrected BAM files are then used for allele extraction.

## Allele extraction and coverage filtering

Allele extraction is performed using FreeBayes through `run_extraction.sh`. The script processes the corrected BAM files and the reference genome to identify genomic positions with multiple alleles and to extract the corresponding allele counts. The resulting files, including `filtered_sites_allc.tsv` and `allele_counts_snps_all.tsv`, provide the basis for subsequent allele-frequency calculations.

An alternative extraction strategy is implemented in `run_extraction_coveragefilter.sh`. In addition to allele extraction, this script applies coverage-based filtering to reduce the influence of positions with insufficient or unusually high sequencing coverage. The resulting coverage-filtered site tables and allele-count files are used in the subsequent analyses. Both extraction strategies identify positions with at least three alleles, but differ in their treatment of sequencing coverage.

## Conversion into allele-frequency trajectories

The extracted allele counts are converted into allele-frequency tables using a set of Python scripts. These scripts calculate allele frequencies across the sampled time points and organize the data into trajectories that describe temporal changes in allele frequencies within the population. The population identifier and the corresponding generation numbers are incorporated during this step.

The conversion is performed using alternative filtering strategies. `data_conversion_all.py` processes the extraction output using hard-filtering criteria, including a minimum read count and a maximum number of alleles per site. In contrast, `data_conversion_hcf.py` processes coverage-filtered data. The scripts `data_conversion_hcf_mna5_allc_new.py` and `data_conversion_cov_morethanfive.py` apply additional restrictions on the number of alleles retained, allowing subsets with different allelic complexity to be analyzed separately.

The resulting allele-frequency trajectories constitute the central input for the downstream analyses. They are used to characterize multi-allelic variation, associate genomic positions with annotated features, and estimate fitness matrices from changes in allele frequencies over time.

## Characterization of multi-allelic sites

The script `count_alleles.py` summarizes the number of alleles observed at each genomic position and generates a distribution of sites according to their allelic count. This analysis provides an overview of the multi-allelic variation retained after the respective filtering and conversion procedures. The results can be used to compare the number and distribution of multi-allelic sites between datasets generated with different filtering strategies.

A separate coverage analysis examines sequencing depth at multi-allelic positions. The script `alleles_coverage_count.py` calculates coverage-related measures, including total and mean coverage across sampling points, and records the number of alleles per site. This produces a site-level coverage summary that can be used to investigate whether the occurrence of multi-allelic sites is associated with differences in sequencing coverage.

## Genomic annotation of multi-allelic sites

To investigate the genomic context of multi-allelic variation, the allele-frequency data are converted into BED format using `makebed.py`. This representation enables the genomic positions to be compared with gene annotations and coding-sequence coordinates.

The scripts `genefinder.py` and `cdsfinder.py` identify genes and coding sequences overlapping the analyzed positions. This step establishes whether multi-allelic sites are located within annotated genes or coding regions and provides gene-level information for subsequent summaries.

A more detailed annotation is performed by `annotationfinder.py`, which combines the genomic coordinates, the GFF3 gene annotation, and the converted allele-frequency data. The script assigns genomic annotations to multi-allelic positions and generates `annotated_multiallelic_sites.tsv`. This output enables the distribution of multiallelic sites across annotated genomic features to be examined.

For the coverage analysis, `genefinder_for_coverage_allinone.py` combines the site-level coverage summaries with the gene annotation. The resulting `sites_with_genes.tsv` links coverage information to genomic positions and their associated genes. Subsequently, `gene_summary2.py` summarizes the occurrence of multi-allelic sites at the gene level, including site counts, gene length, site density, and coverage statistics. This allows differences in the number of multi-allelic sites between genes to be evaluated in relation to gene length and sequencing coverage.

### Fitness matrix estimation

The temporal allele-frequency trajectories are used to estimate fitness matrices, with the aim of inferring relative genotype-specific fitness effects from changes in allele frequencies across generations. Two alternative implementations are available, based on the methods described by [Vellnow et al. (2026)](https://link.springer.com/article/10.1186/s12864-026-13411-5) and the [bioRxiv preprint (2024)](https://www.biorxiv.org/content/10.1101/2024.11.08.622587v3). The original implementation and the revised implementation differ in their cost function.

The original implementation is initiated by `scan_dros_genome_F.m` and uses `fit_F_sparse_func.m` to estimate fitness matrices, which uses `CostSparseFn.m` to calculate the cost by comparing observed allele frequencies with frequencies predicted by the fitness model across the sampled generations. It accounts for gaps between sampling time points by iteratively predicting the intermediate allele frequencies and sums the squared Euclidean differences between predicted and observed frequencies. The revised implementation uses `scan_dros_genome_F_new.m` and `fit_F_sparse_func_new.m`. The latter incorporates a modified cost function implemented in `CostSparseShahRidge.m`, which combines Shahshahani-weighted optimization with ridge regularization. Both implementations use the allele-frequency trajectories, generation information, and the number of alleles to estimate a fitness matrix and evaluate the associated optimization cost.

The supporting function `SMatVec.m` reconstructs a symmetric fitness matrix from a vector containing its independent elements, thereby reducing the number of parameters that need to be optimized. The supporting function `convert_F_to_F_ij.m` transforms the estimated fitness-effect matrix \(F\) into an evolutionarily equivalent matrix \(F_{ij}\), in which a specified matrix entry \((i,j)\) is set to zero. This transformation provides an alternative representation of the fitness matrix while preserving its evolutionary equivalence.
The function `convert_F_to_w.m` converts the fitness-effect matrix \(F\) into a matrix of relative fitness values \(w\). It normalizes the values so that the genotype with the highest fitness has a relative fitness of 1. These relative fitness values can then be used to interpret and compare the fitness of different genotypes.

The resulting fitness matrices and cost values are saved in MATLAB result files and provide the input for subsequent fitness analysis. Keeping the outputs of the two implementations separate allows their results to be evaluated and compared without mixing the underlying optimization approaches.

### Analysis of heterozygote advantage

The estimated fitness matrices are further processed using `scan_dros_genome_anaylsis.m` or `scan_dros_genome_anaylsis_nice.m`. These scripts evaluate the estimated fitness relationships and calculate additional metrics, including `hav1`, which is used in the workflow to quantify the strength of heterozygote advantage.

The resulting analysis tables contain the derived fitness-related metrics and form the basis for further investigations of the evolutionary properties of multiallelic sites. By linking these estimates to the genomic annotations, the fitness-related results can subsequently be examined in the context of specific genes and genomic regions.

### Analysis of insecticide-resistance genes

Finally, the fitness-analysis results are examined for selected candidate insecticide-resistance genes. The scripts `rdl_filter_F_analysis_Results.m` and the corresponding scripts for `ace` and `cyp` extract results associated with specified genomic regions or candidate genes from the broader analysis output.

This targeted analysis makes it possible to investigate fitness-related metrics in selected resistance genes separately from the genome-wide dataset. The resulting subsets can be compared with the genome-wide results to assess whether the selected regions show patterns that warrant further investigation. The interpretation depends on the genomic coordinates and filtering criteria defined in the respective scripts.
