# Author: Colin Tang
#
# CRISPR screen pipeline for data normalization, QC, and basic analysis
# the pipeline requires a plasmid count file for normalization and a libSeqFile containing guide UID and sequence.
# The inputs are the file name for the libSeqFile and the label for the file the screens will be normalized to.
# will process and graph all .results files in the working directory.
# The file names should be according to what you want the labels on graphs to be DesiredLabel.fq.trim.results

library(reshape2)

crispranalysisv1.0<-function(libSeqFile, normalizationFile, sgnontarget = c("luciferase", "EGFP", "LacZ")) {
  ######normalization steps
  tkov3_index <- read.table(libSeqFile, header = TRUE, colClasses = "character", stringsAsFactors = FALSE)
  tkov3_strings <- strsplit(tkov3_index$UID, '_')
  # Convert to data frame and label columns
  metaData <- data.frame(matrix(unlist(tkov3_strings), nrow = nrow(tkov3_index), byrow = TRUE), stringsAsFactors = FALSE)
  metaData <- data.frame(position = metaData[,1], strand = metaData[,3], gene = metaData[,2])
  tkov3_index <- cbind(tkov3_index, metaData)
  
  # Import bowtie results and merge with data frame
  #bowtie <- list.files(pattern = "\\.results$")
  bowtie <- list.files(pattern = "\\.results$")
  for (j in 1:length(bowtie)) {
    # Read in results file for this library
    curCounts <- read.table(bowtie[j], colClasses = c("integer", "character"), header = FALSE, stringsAsFactors = FALSE)
    shortName <- strsplit(bowtie[j],'\\.')[[1]][1]   # Rename column with descriptors
    colnames(curCounts) <- c(shortName, "UID")
    tkov3_index <- merge(tkov3_index, curCounts , by.x = "UID", by.y = "UID", all = TRUE)   # Outer join with library index
  }
  # Discard sgRNA with <30 raw reads in plasmid
  normalizationFile <- strsplit(normalizationFile, split = "\\.")[[1]][1]
  normColumn <- which(colnames(tkov3_index) == normalizationFile)
  tkov3_index <- subset(tkov3_index, tkov3_index[,normColumn] >= 30)
  # If reads in data is "NA", set equal to 0
  tkov3_index[,6:(6+length(bowtie)-1)][is.na(tkov3_index[,6:(6+length(bowtie)-1)])] <- 0
  
  #writes read count file will strand data
  write.table(tkov3_index, "read_strand_counts.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  sgrnareaddata <- tkov3_index[5:(5+length(bowtie))]
  write.table(sgrnareaddata, "readcounts.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  
  #log transforms the data
  tkov3_index[,6:(6+length(bowtie)-1)] <- sweep(tkov3_index[,6:(6+length(bowtie)-1)], 2, apply(tkov3_index[,6:(6+length(bowtie)-1)], 2, sum, na.rm = TRUE), "/") * 1e7
  tkov3_index[,6:(6+length(bowtie)-1)] <- log2(tkov3_index[,6:(6+length(bowtie)-1)] + 1)  #adds 1 to prevent inf after log2 transformation
  sgrnadata <- tkov3_index[5:(5+length(bowtie))]
  
  # Substract plasmid rct from end rct to determine log2-fold change over plasmid
  normColumn <- which(colnames(sgrnadata) == normalizationFile)
  log2.fc <- data.frame(gene = sgrnadata$gene)
  log2.fc <- cbind(log2.fc, sgrnadata[,2:length(sgrnadata)] - sgrnadata[,normColumn])
  log2.fc[,normColumn] <- NULL
  
  # Remove non-targeting genes
  nt_list <- sgnontarget   # non-targeting gene list
  nontargeting <- subset(log2.fc, gene %in% nt_list)
  write.table(nontargeting, "non_targeting.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  targeting <- subset(log2.fc, ! gene %in% nt_list)
  write.table(targeting, "targeting_logFC.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Select second-most enriched guide from each gene
  targeting_2nd_high <- aggregate(. ~ gene, data = targeting, function(x) sort(x, decreasing = TRUE)[2])
  write.table(targeting_2nd_high, "targeting_2nd_high.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Select second-most depleted guide from each gene
  targeting_2nd_low <- aggregate(. ~ gene, data = targeting, function(x) sort(x, decreasing = FALSE)[2])
  write.table(targeting_2nd_low, "targeting_2nd_low.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  
  #average all the guides from each gene records the mean and sd
  sgrnaavg <- aggregate(. ~ gene, data = targeting, mean)
  write.table(sgrnaavg, "sgrna_avg.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  sgrnasd <- aggregate(. ~ gene, data = targeting, sd)
  colnames(sgrnasd) <- paste(colnames(sgrnasd), "SD", sep = "_")
  #merges the average dataframe to get avg next to sd
  sgrna.avg.sd <- subset(sgrnaavg, select = c(1))
  for (i in 2:length(sgrnaavg)) {
    sgrna.avg.sd <- cbind(sgrna.avg.sd, 
                          subset(sgrnaavg, select = c(i)), 
                          subset(sgrnasd, select = c(i)))
  }
  write.table(sgrna.avg.sd, "sgrna_avg_sd.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  
  ######performs graphing for single sample graphs or all sample graphs
  screen.boxplot(sgrnaavg)
  screen.readsboxplot(sgrnareaddata)
  screen.corheatmap(sgrnaavg)
  screen.bargraph(sgrna.avg.sd, sortby = (length(sgrna.avg.sd) - 1), top = TRUE, ngenes = 10)
  #for depletion candidates
  depletioncand <- removecontrols(sgrna.avg.sd, nontargeting = sgnontarget)
  write.table(depletioncand, "cegremoved_sgrna_avg_sd.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  screen.bargraph(depletioncand, sortby = (length(sgrna.avg.sd) - 1), top = FALSE)
  for (i in 2:length(sgrnaavg)) {
    screen.sgrnadist(subset(log2.fc, select = c(1,i)), nontargeting = sgnontarget)
    screen.coverage(subset(sgrnareaddata, select= c(1,i)))
    screen.prc(subset(sgrnaavg, select=c(1,i)))
    screen.rankplot(subset(sgrnaavg, select=c(1,i)))
  }
}

#function to remove controls from dataframe. returns df with controls removed
removecontrols<-function(datatable, nontargeting) {
  nt_list <- nontargeting   # non-targeting gene list
  ceg <- read.table("CoreEssentialGenes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  neg <- read.table("NonEssentialGenes.txt", header = TRUE, stringsAsFactors = FALSE)
  cegindex <- which(datatable[,1] %in% ceg[,1])
  datatable <- datatable[-cegindex,]
  negindex <- which(datatable[,1] %in% neg[,1])
  datatable <- datatable[-negindex,]
  return(datatable)
}


