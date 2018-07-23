# summary.R
# Description: Generate a summary of association results including quantile-quantile and manhattan plots for all groups. Also generates CSV files of all associations.
# Inputs:
# label : prefix for output filename (string)
# assoc.files : comma separated list of association results, output of aggAssocTest (string)

# Check if required packages are installed
packages <- c("qqman","data.table","stringr")
to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install,repos='http://cran.us.r-project.org')

# Load packages
lapply(packages, library, character.only = TRUE)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)
label <- input_args[1]
assoc.files <- unlist(strsplit(input_args[2],","))
minmac <- as.numeric(input_args[3])
var.info <- list()

# remove any empty files because R apparently can't handle loading an empty file.
info = file.info(assoc.files)
assoc.files = rownames(info[!(info$size == 0), ])

if (length(assoc.files) == 0){
  fwrite(list(),paste(label, ".assoc.csv", sep=""),sep=",",row.names=F)
  fwrite(list(), paste(label, ".topassoc.csv", sep=""),row.names=F)
  png(filename = paste(label,"_association_plots.png",sep=""),width = 11, height = 11, units = "in", res=800, type = "cairo")
  dev.off()
  
} else {

  # Prep for association files
  assoc.compilation <- c() 
  
  # Loop through association files
  for (i in seq(1,length(assoc.files))) {
    load(assoc.files[i])
    
    # get group results
    res <- assoc$results
    
    # get variant info for each group
    var.list <- list()
    for (g in names(assoc$variantInfo)){
      df <- assoc$variantInfo[[g]]
      if (NROW(df) == 0){
        next
      }
      df$group <- g
      df$pval <- res[rownames(res) == g,"pval_0"]
      var.list[[g]] <- df
    }
    
    var.df <- do.call(rbind,var.list)
    var.info[[i]] <- var.df

    if (!is.na(res)[1]){
      print(dim(res))
      res <- res[!is.na(res[,"pval_0"]),]
      
      # add chromosome and minimum position to results
      chrs <- c()
      pos <- c()
      for (g in rownames(res)){
        chrs <- c(chrs,var.df[var.df$group == g,]$chr[1])
        pos <- c(pos,min(var.df[var.df$group == g,]$pos))
      }
      
      res$chr <- chrs
      res$pos <- pos
      
      assoc.compilation <- rbind(assoc.compilation, res)
      
      if (i == 1) {
        write.table(res,paste(label, ".assoc.csv", sep=""),sep=",",row.names=F)
      } else {
        write.table(res,paste(label, ".assoc.csv", sep=""),col.names=FALSE,sep=",",row.names=F, append=TRUE)
      } 
    }
  }
}

assoc.compilation <- assoc.compilation[!is.na(assoc.compilation$pval_0),]
assoc.compilation$chr <- ifelse(assoc.compilation$chr == "X", 23, as.numeric(assoc.compilation$chr))

# Calculate genomic control
lam = function(x,p=.5){
  x = x[!is.na(x)]
  chisq <- qchisq(1-x,1)
  round((quantile(chisq,p)/qchisq(p,1)),2)
}

png(filename = paste(label,".association.plots.png",sep=""),width = 11, height = 22, units = "in", res=800)#, type = "cairo")
par(mfrow=c(2,1))

# qq plot
qq(as.numeric(assoc.compilation[,"pval_0"]),main=label)
legend('topleft',c(paste0("GC = ", lam(assoc.compilation[,"pval_0"]))))
manhattan(assoc.compilation,chr="chr",bp="pos",p="pval_0", main=label)
dev.off()

write.csv(assoc.compilation, paste(label, ".groupAssoc.csv", sep=""))
var.df = do.call(rbind,var.info)
fwrite(var.df, file = paste(label, ".all.variants.groupAssoc.csv", sep=""), row.names = F)

# now make the mac thresholded plots
var.df$mac <- 2*var.df$n.obs*var.df$freq

cummac <- data.frame(group = row.names(assoc.compilation), cummac = NA, stringsAsFactors = F)
# get cum mac per group
for (gind in seq(1,nrow(assoc.compilation))){
  g <- row.names(assoc.compilation)[gind]
  cummac[cummac$group == g, ]$cummac <- sum(var.df[var.df$group == g,]$mac)
}
assoc.compilation$group <- row.names(assoc.compilation)
assoc.compilation <- merge(assoc.compilation, cummac, by.x = "group", by.y = "group", all.x = T)
assoc.compilation <- assoc.compilation[assoc.compilation$cummac > minmac,]

png(filename = paste(label,".cummac.",as.character(minmac),".association.plots.png",sep=""),width = 11, height = 22, units = "in", res=800)#, type = "cairo")
par(mfrow=c(2,1))

# qq plot
qq(as.numeric(assoc.compilation[,"pval_0"]),main=label)
legend('topleft',c(paste0("GC = ", lam(assoc.compilation[,"pval_0"]))))
manhattan(assoc.compilation,chr="chr",bp="pos",p="pval_0", main=label)
dev.off()

write.csv(assoc.compilation, paste(label,".cummac.",as.character(minmac), ".groupAssoc.csv", sep=""))
var.df = do.call(rbind,var.info)
fwrite(var.df, file = paste(label,".cummac.",as.character(minmac), ".all.variants.groupAssoc.csv", sep=""), row.names = F)

