# Aggregate Association -- SKAT or Burden tests

## Description 

This workflow performs an aggregate association analysis of genotype data with a single phenotype. The primary code is written in R using the GENESIS package for model fitting and association testing. The workflow can either generate a null model from phenotype and relatedness data or use a pregenerated null model.

### Authors

This workflow is produced and maintained by the [Manning Lab](https://manning-lab.github.io/). Contributing authors include:

* Tim Majarian (tmajaria@broadinstitute.org)
* Alisa Manning (amanning@broadinstitute.org).

## Dependencies

### Workflow execution

* [WDL](https://software.broadinstitute.org/wdl/documentation/quickstart)
* [Chromwell](http://cromwell.readthedocs.io/en/develop/)

### R packages

* [GENESIS](https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html)
* [GWASTools](https://www.bioconductor.org/packages/release/bioc/html/GWASTools.html)
* [SeqArray](https://www.bioconductor.org/packages/release/bioc/html/SeqArray.html)
* [SeqVarTools](https://www.bioconductor.org/packages/release/bioc/html/SeqVarTools.html)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [qqman](https://cran.r-project.org/web/packages/qqman/index.html)

## Main Functions

### fitNull

This function generates a null model to be used in association testing in Genesis

Inputs:
* genotype_files : genotype data for all samples (array of VCF or GDS file)
* phenotype_file : phenotype data for all samples to be included in analysis (CSV or TSV file)
* outcome : the outcome to be tested (string)
* outcome_type : the type of outcome being tested (dichotomous or continuous)
* covariates_string : covariates to condition on in linear mixed modeling (comma separated string, default = '')
* conditional_string : comma separated list of chromosome:position in the format: "chr1:100023,chr10:235225" for conditional analysis (string, optional)
* ivars_string : comma separated list of interaction variables (string, optional)
* sample_file : a file containing a list of sample ids (matching the genotype and phenotype files) to be included, one per line (.txt, optional)
* label : prefix for output filename (string)
* kinship_file : relatedness measures for all samples (CSV or TSV file)
* id_col : column name of id column in phenotype file (string)

Outputs:
* model : generated null model (.RDa)

###  aggAssocTest 

This function performs an association test to generate p-values for each variant included.
	File gds_file
	File null_file
	File group_file
	String label
	String? test
	String? pval
	String? weights

Inputs:
* gds_file : a genotype file containing data for all samples are variants to be tested (.gds)
* null_file : output of the function *fitNull* or a pregenerated null model (.RDa)
* group_file : RData or csv file with groups to include in analysis, if RData, must be saved as a list with unique names, each entry as a data frame with at least variant ids and allele. If csv, must have at least columns for group_id, position, alt, ref (csv or RData)
* label : prefix for output filename (string)
* test : SKAT or Burden (string)
* pval : if SKAT: davies, kuonen, or liu; if Burden: Score, Wald, or Firth (string)
* weights : parameters of beta distribution for variant weights (comma separated string in form: "1,25")

Outputs:
* assoc : an RData file of associations results (.RData)

### summary

Generate a summary of association results including quantile-quantile and manhattan plots for variants subseted by minor allele frequency (all variants, maf < 5%, maf >= 5%). Also generates CSV files of all and the top associated variants.

Inputs:
* label : prefix for output filename (string)
* assoc_files : comma separated list of association results, output of aggAssocTest (string)

## Other workflow inputs

* this_memory : amount of memory in GB for each execution of a task (int)
* this_disk : amount of disk space in GB to allot for each execution of a task (int)



