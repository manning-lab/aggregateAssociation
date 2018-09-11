# General script for generating aggregation units

# Inputs:
## gds.file
## region.file
## variant.file
## filter.terms
## filter.values
## out.pref

# Outputs:
## group.file

args <- commandArgs(trailingOnly=T)
gds.file <- args[1]
region.file <- args[2]
variant.file <- args[3]
max.maf <- args[4]
min.maf <- args[5]
out.pref <- args[6]
index <- args[7]

####################################################################################################################################
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr17.pass_and_fail.gtonly.minDP10.gds"
# region.file <- "/Users/tmajaria/Documents/projects/topmed/data/ferrer/agg_regions/ferrer.islet.enhancer.hub.centric.aggregation.regions.bed"
# variant.file <- "/Users/tmajaria/Documents/projects/topmed/results/rarevar/freeze5b/ptv/ptv_only/freeze5b_dp10_ptv_mask1.tsv"
# max.maf <- 0.01
# min.maf <- 0.000000001
# out.pref <- "testing"
####################################################################################################################################


# Load packages and functions
library(SeqVarTools)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(data.table)


.variantDF <- function(gds) {
  data.frame(variant.id = seqGetData(gds, "variant.id"),
             chr = seqGetData(gds, "chromosome"),
             pos = seqGetData(gds, "position"),
             ref = refChar(gds),
             alt = altChar(gds),
             filter = seqGetData(gds, "annotation/filter"),
             maf = seqAlleleFreq(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    group_by_("variant.id") %>%
    as.data.frame()
}
#########################

# Load regions
region.data <- fread(region.file, data.table = F, stringsAsFactors = F)

if (ncol(region.data) < 5){
  region.data$V5 <- NA
}

# To granges to filter genotypes
region.gr <- GRanges(
  seqnames = region.data$V1,
  ranges = IRanges(
    start = as.numeric(region.data$V2),
    end = as.numeric(region.data$V3)
  ),
  group_id = region.data$V4,
  annotation = region.data$V5
)

# Open genotype file
gds.data <- seqOpen(gds.file)

# Filter by the regions
seqSetFilter(gds.data, region.gr)

# Get data frame of variants
var.region <- .expandAlleles(gds.data)
var.region$maf <- pmin(var.region$maf, 1-var.region$maf)

# Subset by passing variants
var.region <- var.region[var.region$filter == "PASS",]

# Remove mono-allele variants
var.region <- var.region[var.region$maf > 0,]

# Check for maf filter argument and filter if input
if (min.maf != "NA") { var.region <- var.region[var.region$maf > as.numeric(min.maf),] }
if (max.maf != "NA") { var.region <- var.region[var.region$maf < as.numeric(max.maf),] }

# Get the annotations from the region file
# correct chr if nec
if (startsWith(var.region$chr[1],"chr") & !(startsWith(region.data$V1[1],"chr"))){
  var.region$chr <- sub("chr","",var.region$chr)
} else if (!(startsWith(var.region$chr[1],"chr")) & startsWith(region.data$V1[1],"chr")){
  var.region$chr <- sub("^","chr",var.region$chr)
}

var.gr <- GRanges(
  seqnames = var.region$chr,
  ranges = IRanges(
    start = var.region$pos,
    end = var.region$pos
  ),
  ref = var.region$ref,
  alt = var.region$alt,
  maf = var.region$maf
)

var.region.ovp <- findOverlaps(var.gr, region.gr)
var.gr <- var.gr[queryHits(var.region.ovp),]
var.gr$group_id <- region.gr[subjectHits(var.region.ovp),]$group_id
var.gr$annotation <- region.gr[subjectHits(var.region.ovp),]$annotation

var.region <- unique(
  data.frame(
    chr = seqnames(var.gr),
    pos = start(var.gr),
    ref = var.gr$ref,
    alt = var.gr$alt,
    maf = var.gr$maf,
    group_id = var.gr$group_id,
    annotation = var.gr$annotation,
    stringsAsFactors = F
  )
)

# Reset filter
seqResetFilter(gds.data)

# Check if we have a variant file
# if we do, add it to the region data
if (variant.file != "NA"){
  variant.data <- fread(variant.file, data.table = F, stringsAsFactors = F)
  if (all(c("chromosome", "position", "ref", "alt", "group_id") %in% names(variant.data))){
    other.col <- names(variant.data)[!(names(variant.data) %in% c("chromosome", "position", "ref", "alt", "group_id"))]
    variant.data <- variant.data[, c("chromosome", "position", "ref", "alt", "group_id", other.col)]
  }
  names(variant.data) <- paste0("V", seq(1,length(names(variant.data))))
  if (!("V6" %in% names(variant.data))) { variant.data$V6 <- NA }
  variant.gr <- GRanges(
    seqnames = variant.data$V1,
    ranges = IRanges(
      start = as.numeric(variant.data$V2),
      end = as.numeric(variant.data$V2)
    ),
    ref = variant.data$V3,
    alt = variant.data$V4,
    group = variant.data$V5,
    annotation = variant.data$V6
  )
  
  seqSetFilter(gds.data, variant.gr)
  var.variants <- .expandAlleles(gds.data)
  var.variants$maf <- pmin(var.variants$maf, 1-var.variants$maf)
  
  # correct chr if nec
  if (startsWith(var.variants$chr[1],"chr") & !(startsWith(variant.data$V1[1],"chr"))){
    var.variants$chr <- sub("chr","",var.variants$chr)
  } else if (!(startsWith(var.variants$chr[1],"chr")) & startsWith(variant.data$V1[1],"chr")){
    var.variants$chr <- sub("^","chr",var.variants$chr)
  }
  
  var.variants <- merge(variant.data, var.variants, by.x = c("V1", "V2", "V3", "V4"), by.y = c("chr", "pos", "ref", "alt"))
  var.variants <- var.variants[,c(1,2,3,4,9,5,6)]
  names(var.variants) <- c("chr","pos","ref","alt","maf","group_id","annotation")
  var.region <- rbind(var.region, var.variants)
}

# close genotype file
seqClose(gds.data)

# write out group file
write.table(var.region, file = paste0(out.pref,".",index,".csv"), quote = F, sep = ",", row.names = F)
