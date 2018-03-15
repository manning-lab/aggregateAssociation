### aggAssociation.R
# Description: This function performs an aggregate association test using the skat or burden methods.
# Inputs:
# gds.file : a genotype file containing data for all samples are variants to be tested (.gds)
# null.file : output of the function *fitNull* or a pregenerated null model (.RDa)
# group.file : RData or csv with aggregation groups. Must contain at least variant id and group id (csv or RData)
# label : prefix for output filename (string)
# test : SKAT or Burden (string)
# pval : if SKAT: davies, kuonen, or liu; if Burden: Score, Wald, or Firth (string)
# weights : parameters of beta distribution for variant weights (comma separated string in form: "1,25")
# Outputs:
# assoc : an RData file of associations results (.RData)

# Load packages
library(GENESIS)
library(SeqArray)
library(SeqVarTools)
library(data.table)
library(Biobase)

# Parse input args
input_args <- commandArgs(trailingOnly=T)
gds.file <- input_args[1] 
null.file <- input_args[2]
group.file <- input_args[3]
label <- input_args[4]
test <- input_args[5]
pval <- input_args[6]
weights <- as.numeric(unlist(strsplit(input_args[7],",")))

# these are from the DCC pipeline, credit -> S. Gogarten
getAggList <- function(gds, variants.orig){
  filtOrig <- seqGetFilter(gds)
  seqSetFilter(gds, variant.id=variants.orig)
  variants.new <- .expandAlleles(gds)
  group <- data.frame(variant.id=variants.new$variant.id, allele.index=variants.new$allele.index)
}

.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             chromosome=seqGetData(gds, "chromosome"),
             position=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             nAlleles=seqNumAllele(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    rename_(allele="alt") %>%
    group_by_("variant.id") %>%
    mutate_(allele.index=~1:n()) %>%
    as.data.frame()
}

# Load nullfile
load(null.file)

# Open gds file
gds.data <- seqOpen(gds.file)

# Filter to samples in null model
seqSetFilter(gds.data,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)

# Genotype data to the correct format
gds.geno.data <- SeqVarData(gds.data)

## load groups
group_ext <- tail(unlist(strsplit(basename(group.file),'\\.')),n=1)
if (group_ext == 'RData'){
  # load if RData file
  load(group.file)  

} else if (group_ext == 'tsv' | group_ext == 'csv') {
  # load with data table
  group.raw <- fread(group.file, data.table=F)
  var.df <- data.frame(variant.id = seqGetData(gds.data, "variant.id"), pos = seqGetData(gds.data, "position"))
  var.df <- var.df[var.df$pos %in% group.raw$position,"variant.id"]
  seqSetFilter(gds.data,variant.id=var.df)
  
  library(SeqVarTools)
  library(dplyr)
  library(tidyr)
  var.df <- .expandAlleles(gds.data)
  var.df <- merge(group.raw, var.df, by.x=c('position','ref','alt'), by.y=c('position','ref','allele'))
  seqSetFilter(gds.data,variant.id=var.df$variant.id)
  var.df$maf <- alleleFrequency(gds.data,n=1)
  
  
  # var.df <- var.df[var.df$pos %in% group.raw$position,]
  # group.raw <- group.raw[group.raw$position %in% var.df$pos,]
  # group.var <- merge(group.raw, var.df, by.x=c('position','ref','alt'), by.y=c('pos','ref','alt'))
  # 
  # seqSetFilter(gds.data,variant.id = group.var$variant.id)
  # library(SeqVarTools)
  # library(dplyr)
  # library(tidyr)
  # gds.df <- .expandAlleles(gds.data)
  # gds.df$maf <- alleleFrequency(gds.data,n=1)
  groups <- list()
  
  for (gid in unique(var.df$group_id)){
    groups[[as.character(gid)]] <- var.df[var.df$group_id == gid,]
  }
  
} else {
  stop("Group file does not have the required extension")
}

# make sure we have no duplicates
groups = groups[!duplicated(names(groups))]

# make sure all groups are in the gds file
# groups.var_id <- do.call(rbind, groups)$variant.id
# if (is.null(gds.geno.data@variantData@data$variant.id)){
#   var.data <- data.frame(variant.id = seqGetData(gds.data, 'variant.id'))
#   var.anno <- AnnotatedDataFrame(var.data)
#   variantData(gds.geno.data) <- var.anno
# }
# 
# if (any(!(groups.var_id %in% gds.geno.data@variantData@data$variant.id))){
#   stop("One or more groups contain variants that are not in the genotype file")
# }

#### run association test
if(tolower(test)=="skat"){
  assoc <- assocTestSeq(gds.geno.data, nullmod, groups, test="SKAT", pval.method=pval, weight.beta = weights)
  assoc$results = assoc$results[order(assoc$results$pval_0),]
  for (group_name in names(assoc$variantInfo)){
    assoc$results[group_name,"MAF"] <- mean(assoc$variantInfo[[group_name]]$freq)
  }
  save(assoc, file=paste(label, ".assoc.RData", sep=""))
} else if (tolower(test) == "burden") {
  assoc <- assocTestSeq(gds.geno.data, nullmod, groups, test="Burden", burden.test=pval, weight.beta = weights)
  assoc$results = assoc$results[order(assoc$results[,paste(pval,".pval",sep="")]),]
  for (group_name in names(assoc$variantInfo)){
    assoc$results[group_name,"MAF"] <- mean(assoc$variantInfo[[group_name]]$freq)
  }
  names(assoc$results)[names(assoc$results) == paste(pval,".pval",sep="")] = "pval_0"
  save(assoc, file=paste(label, ".assoc.RData", sep=""))
  
} else {
  fwrite(list(), file=paste(label, ".assoc.RData", sep=""))
}

seqClose(gds.data)





