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
force.maf <- input_args[8]

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
  # load group file
  group.raw <- fread(group.file, data.table=F)
  
  # need to match with variant data; get variant.id and position
  var.df <- data.frame(variant.id = seqGetData(gds.data, "variant.id"), pos = seqGetData(gds.data, "position"))
  
  # subset by positions in group file
  var.df <- var.df[var.df$pos %in% group.raw$position,"variant.id"]
  seqSetFilter(gds.data,variant.id=var.df)
  
  # get variant info in correct format
  library(SeqVarTools)
  library(dplyr)
  library(tidyr)
  var.df <- .expandAlleles(gds.data)
  
  # merge over pos, ref, and alt
  var.df <- merge(group.raw, var.df, by.x=c('position','ref','alt'), by.y=c('position','ref','allele'))
  
  # filter to those variants in both gds and groups
  seqSetFilter(gds.data,variant.id=var.df$variant.id)
  
  # add the maf to the groups
  ref.freq <- seqAlleleFreq(gds.data, ref.allele=0L, .progress = T)
  
  # force alt to be lower maf allele
  if (force.maf == "False" | force.maf == "F"){
    var.df$maf <- ref.freq
  } else {
    # get the right alt index
    alt.id <- data.frame(variant.id = seqGetData(gds.data,"variant.id"), alt.index = ifelse(ref.freq < 1-ref.freq, 0, 1), maf = pmin(ref.freq, 1-ref.freq))
    
    # merge with group file
    var.df <- merge(var.df,alt.id, by.x = "variant.id", by.y = "variant.id")
    
    # rename cols
    var.df <- var.df[,names(var.df)[names(var.df)!= "allele.index"]]
    names(var.df)[names(var.df)== "alt.index"] <- "allele.index"
  }
  
  groups <- list()
  
  for (gid in unique(var.df$group_id)){
    groups[[as.character(gid)]] <- var.df[var.df$group_id == gid,]
  }
} else {
  stop("Group file does not have the required extension")
}

# make sure we have no duplicates
groups = groups[!duplicated(names(groups))]

# groups to data frame
groups.df <- do.call(rbind,groups)

#### run association test
if(tolower(test)=="skat"){
  assoc <- assocTestSeq(gds.geno.data, nullmod, groups, test="SKAT", pval.method=pval, weight.beta = weights)
  assoc$results = assoc$results[order(assoc$results$pval_0),]
  for (group_name in names(assoc$variantInfo)){
    assoc$results[group_name,"MAF"] <- mean(assoc$variantInfo[[group_name]]$freq)
  }
  
  # add ref/alt to assoc variants
  for (g in names(assoc$variantInfo)){
    assoc$variantInfo[[g]] <- merge(assoc$variantInfo[[g]], groups.df[,c("variant.id","ref", "alt", "position","allele.index")], by.x = "variantID", by.y = "variant.id")
  }
  
  save(assoc, file=paste(label, ".assoc.RData", sep=""))
  save(groups, file=paste(label, "groups.RData", sep=""))
} else if (tolower(test) == "burden") {
  assoc <- assocTestSeq(gds.geno.data, nullmod, groups, test="Burden", burden.test=pval, weight.beta = weights)
  assoc$results = assoc$results[order(assoc$results[,paste(pval,".pval",sep="")]),]
  for (group_name in names(assoc$variantInfo)){
    assoc$results[group_name,"MAF"] <- mean(assoc$variantInfo[[group_name]]$freq)
  }
  names(assoc$results)[names(assoc$results) == paste(pval,".pval",sep="")] = "pval_0"
  
  # add ref/alt to assoc variants
  for (g in names(assoc$variantInfo)){
    assoc$variantInfo[[g]] <- merge(assoc$variantInfo[[g]], groups.df[,c("variant.id","ref", "alt", "position","allele.index")], by.x = "variantID", by.y = "variant.id")
  }
  
  save(assoc, file=paste(label, ".assoc.RData", sep=""))
  save(groups, file=paste(label, "groups.RData", sep=""))
  
} else {
  fwrite(list(), file=paste(label, ".assoc.RData", sep=""))
  save(groups, file=paste(label, "groups.RData", sep=""))
}

seqClose(gds.data)





