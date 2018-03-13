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
all_assoc <- list()

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
    
    res <- assoc$results
    all_assoc[[i]]<-assoc

    if (!is.na(res)[1]){
      print(dim(res))
      res <- res[!is.na(res[,"pval_0"]),]
      
      #add to assoc.compilation
      res <- rbind(res,rep(assoc$variantInfo[[1]]$chr,length(res[,1])))
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

png(filename = paste(label,"_association_plots.png",sep=""),width = 11, height = 11, units = "in", res=800, type = "cairo")
par(mfrow=c(2,1))

# qq plot
qq(as.numeric(assoc.compilation[,"pval_0"]),main=label)

# mh plot
l <- list()
 for (i in seq(1,length(all_assoc))){
   l2 <- list()
     for (j in seq(1,length(all_assoc[[i]]$variantInfo))){
       l2[[length(l2)+1]] <- data.frame(P=rep(all_assoc[[i]]$results$pval_0[j],length(all_assoc[[i]]$variantInfo[[j]][,1])), BP=all_assoc[[i]]$variantInfo[[j]]$pos, CHR=all_assoc[[i]]$variantInfo[[j]]$chr)
     }
   l <- unlist(list(l,l2),recursive=F)
 }

 df <- l[[1]]
 for (i in seq(2,length(l))){
   df <- rbind(df,l[[i]])
 }
df$CHR <- as.numeric(as.vector(df$CHR))

manhattan(df,chr="CHR",bp="BP",p="P", main=label)
dev.off()

write.csv(assoc.compilation, paste(label, ".groupAssoc.csv", sep=""))

