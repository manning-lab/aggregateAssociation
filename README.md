# Aggregate Association -- SKAT or Burden tests

## Description 

This workflow performs an aggregate association analysis of genotype data with a single phenotype. The primary code is written in R using the GENESIS package for model fitting and association testing. The workflow can either generate a null model from phenotype and relatedness data or use a pregenerated null model. This workflow outputs both the full association statistics (pvalues, test statistics) and summary plots (manhatten and QQ)

### Authors

This workflow is produced and maintained by the [Manning Lab](https://manning-lab.github.io/). Contributing authors include:

* Tim Majarian (tmajaria@broadinstitute.org)
* Alisa Manning (amanning@broadinstitute.org).

## Dependencies

### Workflow execution

This workflow is written with a WDL wrapper for portability. Each task can be run individually from the base script or within the full workflow.

* [WDL](https://software.broadinstitute.org/wdl/documentation/quickstart)
* [Chromwell](http://cromwell.readthedocs.io/en/develop/)

### R packages

All code and necessary packages are available through a [Docker image](https://hub.docker.com/r/manninglab/aggregateassociation/) as well as through the [Github repository](https://github.com/manning-lab/aggregateAssociation).

* [GENESIS](https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html)
* [GWASTools](https://www.bioconductor.org/packages/release/bioc/html/GWASTools.html)
* [SeqArray](https://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)
* [SeqVarTools](https://www.bioconductor.org/packages/release/bioc/html/SeqVarTools.html)
* [Biobase](https://bioconductor.org/packages/release/bioc/html/Biobase.html)
* [dplyr](https://dplyr.tidyverse.org/)
* [tidyr](https://tidyr.tidyverse.org/)
* [stringr](https://cran.r-project.org/web/packages/stringr/vignettes/stringr.html)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [qqman](https://cran.r-project.org/web/packages/qqman/index.html)

## Main Functions

*Italics* indicate optional inputs. Each function also requires the user to specify allocated memory and disk space.

### fitNull

This function generates a null model to be used in association testing in Genesis.

Inputs:
* genotype_file : genotype data for all samples, this is only used to ensure the correct ordering of phenotype data, it need not contain all genotypes to be tested in aggAssocTest (GDS file)
* *phenotype_file* : phenotype data for all samples to be included in analysis (CSV or TSV file)
* *outcome_name* : the outcome to be tested (string)
* *outcome_type* : the type of outcome being tested (dichotomous or continuous)
* *covariates_string* : covariates to condition on in linear mixed modeling (comma separated string, default = None)
* *sample_file* : a file containing a list of sample ids (matching the genotype and phenotype files) to be included, one per line (.txt, optional)
* label : prefix for output filename (string)
* *kinship_matrix* : relatedness measures for all samples (CSV or TSV file)
* *id_col* : column name of id column in phenotype file (string)

Outputs:
* model : generated null model (.RDa)
* log : log file containing paths to inputs and memory and cpu usage

###  aggAssocTest 

This function performs an association test to generate p-values for each aggregation unit included.

Inputs:
* genotype_file : a genotype file containing data for all samples and all variants to be tested (GDS file)
* *null_file* : output of fitNull or a pregenerated null model (.RDa)
* group_file : RData or csv/tsv file with groups to include in analysis, if RData, must be saved as a list with unique names, each entry as a data frame with at least columns for variant.id, position, chromosome, ref, allele, nAlleles, allele.index. If csv, must have at least columns for group_id, position, chromosome, ref, alt (csv or RData)
* label : prefix for output filename (string)
* *test* : SKAT or Burden (string)
* *pval* : if SKAT: davies, kuonen, or liu; if Burden: Score, Wald, or Firth (string)
* *weights* : parameters of beta distribution for variant weights (comma separated string in form: "1,25") 
* *force_maf* : flag to force the minor allele to the least frequent allele (default = True)

Outputs:
* assoc : an RData file of associations results (.RData)
* log : log file containing paths to inputs and memory and cpu usage
* groups : csv file of actualy groups used for analysis. Some reasonable defaults are set during execution: unless force_maf == F, the alternate allele may be changed so that MAF(alt) < MAF(ref). Some variants may also be removed from groups if they do not occuring within the sample being tested. This will be reflected in the groups output.

### summary

Generate a summary of association results including quantile-quantile and manhattan plots for variants, one each for all and only groups with cumulative minor allele count above the specified threshold. Also generates CSV files of all and the top associated variants.

Inputs:
* label : prefix for output filename (string)
* assoc_files : comma separated list of association results, output of aggAssocTest (string)
* *minmac* : cumulative minor allele count threshold for plotting (Int, default = 10)

Outputs:
* plots : manhatten and QQ plots (PNG)
* assoc_res : association results for all groups (CSV)
* assoc_res_variants : association results for all variants in all groups (CSV)
* mac_plots : thresholded by minmac manhatten and QQ plots (PNG)
* mac_assoc_res : association results for groups thresholded by minmac  (CSV)
* mac_assoc_res_variants : association results for all variants in groups thresholded by minmac  (CSV)
* log : log file containing paths to inputs and memory and cpu usage

## Other workflow inputs

* this_fitNull_memory : amount of memory in GB for fitNull task (int)
* this_aggAssocTest_memory : amount of memory in GB for aggAssocTest task (int)
* this_summary_memory : amount of memory in GB for summary task (int)
* this_disk : amount of disk space in GB to allot for each execution of a task (int)



