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
  seqSetFilter(gds.data, variant.sel = filtOrig$variant.sel, sample.sel = filtOrig$sample.sel)
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
    group_by_("variant.id") %>%
    mutate_(allele.index=~1:n()) %>%
    as.data.frame()
}

.checkHeader <- function(agglist) {
  cols <- colnames(agglist)
  if ( !( any( c( "chr","chromosome","chrom" ) %in% tolower(cols) ) ) ) { return(F) }
  if ( !( any( c( "pos","position" ) %in% tolower(cols) ) ) ) { return(F) }
  if ( !( any( c( "ref","reference","allele1","a1" ) %in% tolower(cols) ) ) ) { return(F) }
  if ( !( any( c( "alt","alternate","allele2","a2" ) %in% tolower(cols) ) ) ) { return(F) }
  if ( !( any( c( "group.id","group_id","group" ) %in% tolower(cols) ) ) ) { return(F) }
  return(T)
}

.fixHeader <- function(agglist) {
  cols <- colnames(agglist)
  cols[tolower(cols) %in% c( "chr","chromosome","chrom" )] <- "chromosome"
  cols[tolower(cols) %in% c( "pos","position" )] <- "position"
  cols[tolower(cols) %in% c( "ref","reference","allele1","a1" )] <- "ref"
  cols[tolower(cols) %in% c( "alt","alternate","allele2","a2" )] <- "alt"
  cols[tolower(cols) %in% c( "group.id","group_id","group" )] <- "group.id"
  colnames(agglist) <- cols
  return( agglist )
}

# Load nullfile
load(null.file)

# Open gds file
gds.data <- seqOpen(gds.file)

# Filter to samples in null model
seqSetFilter(gds.data,sample.id=nullmod$scanID, action="intersect", verbose=TRUE)

# Genotype data to the correct format
gds.geno.data <- SeqVarData(gds.data)

# set flag
has.var <- T

## load groups
group_ext <- tail(unlist(strsplit(basename(group.file),'\\.')),n=1)
if (group_ext == 'RData'){
  # load if RData file
  load(group.file)
  
} else if (group_ext == 'tsv' | group_ext == 'csv') {
  # load group file
  group.raw <- fread(group.file, data.table=F)
  
  # check header of group file
  if ( !( .checkHeader( group.raw ) ) ) {
    stop("Group file must have the minimal column headers of chr,pos,ref,alt,group")
  }
  
  # fix header of group file
  group.raw <- .fixHeader( group.raw )
  
  # check if group file already has required columns
  #  variant.id matching the variant.id in seqData for the variants that should be aggregated, and allele.index
  if (!all(c("group.id","variant.id", "allele.index") %in% names(group.raw))){
    
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
    
    # make sure chromosome in same format
    gds.chr <- seqGetData(gds.data,"chromosome")[1]
    group.chr <- as.character(unique(group.raw$chromosome)[1])
    
    if (startsWith(gds.chr, "chr")){
      if(!(startsWith(group.chr,"chr"))){
        group.raw$chromosome <- sub("^","chr",group.raw$chromosome)
      }
    }
    
    if(!(startsWith(gds.chr,"chr"))){
      if (startsWith(group.chr, "chr")){
        group.raw$chromosome <- sub("chr","", group.raw$chromosome)
      }
    }
    
    # merge over pos, ref, and alt
    var.df <- merge(var.df, group.raw, by.x=c('chromosome','position','ref','alt'), by.y=c('chromosome','position','ref','alt'))
    
    # filter to those variants in both gds and groups
    seqSetFilter(gds.data,variant.id=var.df$variant.id)

    # set has variants flag
    has.var <- T
    
    # stop if we have no variants
    if(length(seqGetData(gds.data, "variant.id")) == 0){
      system(paste("touch ",label, ".assoc.RData", sep=""))
      system(paste("touch ",label, ".groups.RData", sep=""))
      print("No variants were found in genotype file for the input aggregation units, returning empty.")
      has.var <- F
    } else {
    
      # add the maf to the groups
      ref.freq <- seqAlleleFreq(gds.data, ref.allele=0L)
      
      # force alt to be lower maf allele
      if (startsWith(tolower(force.maf), "f")){
        var.df$maf <- ref.freq
      } else {
        # get the right alt index
        alt.id <- data.frame(variant.id = seqGetData(gds.data,"variant.id"), allele.index = ifelse(ref.freq < 1-ref.freq, 0, 1), maf = pmin(ref.freq, 1-ref.freq))
        alt.id$mac <- seqAlleleCount(gds.data, ref.allele = alt.id$allele.index)
        
        # merge with group file
        var.df <- merge(var.df[,names(var.df)[!(names(var.df) %in% c("allele.index","maf"))] ], alt.id, by.x = "variant.id", by.y = "variant.id")
      }
    }
  } else {
    var.df <- group.raw
  }
  
  if (has.var){

    # make the group list for genesis input
    groups <- list()
    
    for (gid in unique(var.df$group.id)){
      groups[[as.character(gid)]] <- var.df[var.df$group.id == gid,]
    }
    
    # remove old df
    rm(var.df)
  }
} else {
  stop("Group file does not have the required extension")
}

if (has.var) {

  # make sure we have no duplicates
  groups <- groups[!duplicated(names(groups))]

  # groups to data frame
  groups.df <- do.call(rbind,groups)

  #### run association test
  if(tolower(test)=="skat"){
    assoc <- assocTestSeq(gds.geno.data, nullmod, groups, test="SKAT", pval.method=pval, weight.beta = weights)
    assoc$results = assoc$results[order(assoc$results$pval_0),]
    for (group_name in names(assoc$variantInfo)){
      assoc$results[group_name,"MAF"] <- mean(assoc$variantInfo[[group_name]]$freq)
      assoc$variantInfo[[group_name]] <- merge(assoc$variantInfo[[group_name]], groups.df[,c("variant.id","ref", "alt", "position","allele.index","mac")], by.x = "variantID", by.y = "variant.id")
      assoc$results[group_name,"MAC"] <- sum(assoc$variantInfo[[group_name]]$mac, na.rm = T)
    }
    
    # add ref/alt to assoc variants
    # for (g in names(assoc$variantInfo)){
      
    # }
    
    save(assoc, file=paste(label, ".assoc.RData", sep=""))
    save(groups, file=paste(label, ".groups.RData", sep=""))
  } else if (tolower(test) == "burden") {
    assoc <- assocTestSeq(gds.geno.data, nullmod, groups, test="Burden", burden.test=pval, weight.beta = weights)
    assoc$results = assoc$results[order(assoc$results[,paste(pval,".pval",sep="")]),]
    for (group_name in names(assoc$variantInfo)){
      assoc$results[group_name,"MAF"] <- mean(assoc$variantInfo[[group_name]]$freq)
    }
    names(assoc$results)[names(assoc$results) == paste(pval,".pval",sep="")] = "pval_0"
    
    # add ref/alt to assoc variants
    for (g in names(assoc$variantInfo)){
      assoc$variantInfo[[g]] <- merge(assoc$variantInfo[[g]], groups.df[,c("variant.id","ref", "alt", "position","allele.index","mac")], by.x = "variantID", by.y = "variant.id")
    }
    
    save(assoc, file=paste(label, ".assoc.RData", sep=""))
    save(groups, file=paste(label, ".groups.RData", sep=""))
    
  } else {
    system(paste("touch ",label, ".assoc.RData", sep=""))
    system(paste("touch ",label, ".groups.RData", sep=""))
  }
}

seqClose(gds.data)





