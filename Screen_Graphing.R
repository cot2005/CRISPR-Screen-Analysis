library("ggplot2")
library("wesanderson")
library("dplyr")
library("ROCR")
library("ggbiplot")
library("rgl")
library("ggrepel")

################################################################################
#Graphing functions called by crispranalysisv1.0
################################################################################
# Script to make grouped bar graph.
# Function will take top boolean, False == bottom
# n = number of genes wanted
# sortby = column number (of desired sample to sort by)
# input dataframe needs to be processed with only targeting genes. Script will take
# top n or bottom n and make a bargraph. data frame needs to be formatted as
# gene, treatment, treatment sd, and so on for any number of samples

screen.bargraph<-function(datatable, ngenes = 20, top = TRUE, sortby = length(datatable)) {
  datatable <- datatable[order(datatable[,sortby], decreasing = top),]
  datatable <- datatable[1:ngenes,]
  datatable[,1] <- factor(datatable[,1], levels = datatable[,1])
  numsamples <- length(datatable)
  ggtable <- data.frame()
  for (i in seq(2,numsamples,2)) {
    temptable <- data.frame(datatable[,1],colnames(datatable)[i], datatable[,i], datatable[,(i+1)] )
    colnames(temptable) <- c("Gene", "Sample", "sgRNA", "SD")
    ggtable <- rbind(ggtable, temptable)
  }
  g <- ggplot(ggtable, aes(x=as.factor(Gene), y=sgRNA, fill=Sample))
  label <- "top"
  if (top == FALSE) {
    label <- "bottom"
    g <- g + scale_y_reverse()
  }
  #scale_fill_manual(values=wes_palette(n=Samples, name="darjeeling")) +
  g + geom_bar(position=position_dodge(), stat="identity", colour='black') +
    geom_errorbar(aes(ymin=sgRNA-SD, ymax=sgRNA+SD), position=position_dodge(.9), width = .5) +
    ylab("sgRNA frequency") + xlab("Gene") + ggtitle("sgRNA Representation in Control") + 
    theme(axis.text.x=element_text(angle=45, hjust=1), 
          panel.background = element_rect(fill = "white", colour = "black", size = 0.7, linetype = "solid"))
  labelname <- paste("bargraph_", label, ngenes, "sortby", sortby,".pdf",sep = "")
  ggsave(labelname, width = 10, height = 7)
}

################################################################################

# Script to make a boxplot for sgRNA represemtation in the samples.

screen.boxplot<-function(datatable) {
  data.reordered <- melt(datatable)
  colnames(data.reordered) <- c("Gene", "Sample", "sgRNA_representation")
  ggdata <- ggplot(data.reordered, aes(x = Sample, y = sgRNA_representation))
  graph <- ggdata + geom_boxplot(outlier.colour = "steelblue", outlier.size = 0.7) + 
    ggtitle("sgRNA Score Boxplot") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.background = element_rect(fill = "white", colour = "black", size = 0.7, linetype = "solid"))
  ggsave("sgRNAscore_boxplot.pdf", width = 10, height = 10)
}


################################################################################

# Script to make a boxplot for sgRNA read distribution in the samples.

screen.readsboxplot<-function(datatable) {
  data.reordered <- melt(datatable)
  colnames(data.reordered) <- c("Gene", "Sample", "Read_coverage")
  ggdata <- ggplot(data.reordered, aes(x = Sample, y = Read_coverage))
  graph <- ggdata + geom_boxplot(outlier.colour = "red", outlier.size = 0.7, outlier.shape = 3) + 
    ggtitle("Read Coverage Boxplot") + 
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "white", colour = "black", size = 0.7, linetype = "solid"))
  ggsave("readcount_boxplot.pdf", width = 10, height = 10)
}


################################################################################

# Script to make read coverage distribution graph.
# Input is dataframe with read count and genes in column 1 then makes graphs for every column in the table.

screen.coverage<-function(datatable) {
  liblength <- (length(datatable[,1]))
  most.reads <- max(datatable[,2])
  readstats <- data.frame(sort(datatable[,2]))
  for (i in 1:most.reads) {
    lessreads <- length(which(datatable[,2] < i))/liblength * 100
    readstats[i,2] <- i
    readstats[i,3] <- lessreads
  }
  pdf(paste(colnames(datatable)[2], "_readcoverage.pdf", sep = ""), 10,7)
  plot(readstats[,2],readstats[,3], log = "x" ,pch = "", xlab = "Number of reads per sgRNA",
       ylab = "Cumulative percentage of Library",
       main = paste("Read coverage per sgRNA in", colnames(datatable)[2]))
  lines(readstats[,2],readstats[,3])
  dev.off()
}


################################################################################

# Script to make a density plot for nontargeting vs essential vs total sgRNA abundance
# Input dataframe needs to have all the guides the format of: sgRNA then normalized score.
# nontargeting input is a list of nontargeting sgRNA gene names in the library (defaults to TKOv3 controls)

screen.sgrnadist<-function(datatable, nontargeting = c("luciferase", "EGFP", "LacZ")) {
  ceg <- read.table("CoreEssentialGenes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  neg <- read.table("NonEssentialGenes.txt", header = TRUE, stringsAsFactors = FALSE)
  cegindex <- which(datatable[,1] %in% ceg[,1])
  ceg.score <- datatable[cegindex,]
  ceg.score <- data.frame(sgRNA = "Essential Genes",Score = ceg.score[,2])
  negindex <- which(datatable[,1] %in% neg[,1])
  neg.score <- datatable[negindex,]
  neg.score <- data.frame(sgRNA = "Non-essential",Score = neg.score[,2])
  total.score <- datatable
  total.score <- total.score[-cegindex,]
  negindex <- which(datatable[,1] %in% neg[,1])
  total.score <- total.score[-negindex,]
  if (is.null(nontargeting) == TRUE) {
    total.score <- data.frame(sgRNA = "Targeting sgRNAs", Score = total.score[,2])
    controlscores <- rbind(ceg.score, neg.score, total.score)
  } else {
    nontargetIndex <- which(datatable[,1] %in% nontargeting)
    nt.score <- datatable[nontargetIndex,]
    nt.score <- data.frame(sgRNA = "Non-targeting", Score = nt.score[,2])
    total.score <- total.score[-nontargetIndex,]   #removes non targeting 
    total.score <- data.frame(sgRNA = "Targeting sgRNAs", Score = total.score[,2])
    controlscores <- rbind(ceg.score, neg.score, total.score, nt.score)
  }
  ggplot(controlscores, aes(x=Score, fill=sgRNA, y = ..density..)) + geom_density(alpha = .25) + 
    ylab("Density") + xlab("Average 2nd Highest sgRNA Frequency") + ggtitle("sgRNA Distribution") + theme_bw()
  ggsave(paste(colnames(datatable)[2], "_sgRNAdist.pdf", sep = ""), height = 5, width = 8)
}


################################################################################

# Script to generate a Precision-Recall curve for screen quality control.
# Generate PR curves with input dataframe with gene colunn and one sgRNA data column only

screen.prc<-function(datatable) {
  # import reference sets: core essential genes, non-essential genes
  ceg <- read.table("CoreEssentialGenes.txt", header = TRUE, stringsAsFactors = FALSE)
  neg <- read.table("NonEssentialGenes.txt", header = TRUE, stringsAsFactors = FALSE)
  ceg$label <- TRUE    # classify and concatenate reference sets
  neg$label <- FALSE
  controls <- data.frame(gene = c(ceg$Gene, neg$Gene), label = c(ceg$label, neg$label))
  
  #trims data to second highest guide
  datatable[,2] <- datatable[,2] * -1   # reverse sign
  geneindex <- match(controls$gene, datatable[,1])   #gets row of controls in data
  controls <- cbind(controls, datatable[geneindex,2])  #puts data in controls df
  controls <- na.omit(controls)   # removes NAs
  
  controlpred <- prediction(controls[,3], controls$label)
  controlperf <- performance(controlpred, "prec", "rec")
  
  pdf(paste(colnames(datatable)[2], "_precision_recall.pdf", sep = ""), height = 7, width = 7)
  controlauc <- performance(controlpred, "auc")@y.values[[1]]
  controlauc <- round(controlauc, digits = 3)
  plot(controlperf, col = "#00C48F", lwd = 3, xaxs = "i", yaxs = "i")
  text(0.2, 0.6,  labels = paste("AUC = ", controlauc), adj = c(0,0),cex =1.5)
  dev.off()
}


################################################################################

# Script to make correlation heatmap (uses cor.test and defaults to spearman).
# Input is in the column format of: gene name then sample sgRNA data.

screen.corheatmap<-function(datatable, method = "spearman", width = 7, height = 5) {
  datatable$gene <- NULL
  numsamples <- length(datatable)
  cortable <- data.frame(colnames(datatable))
  #makes correlation table
  for (i in 1:numsamples) {
    for(j in 1:numsamples) {
      cortable[(j),(i+1)] <- cor.test(datatable[,i],datatable[,j], method = method)$estimate
    }
  }
  colnames(cortable) <- c("Sample", colnames(datatable))
  cortable.reordered <- melt(cortable)
  colnames(cortable.reordered) <- c("SampleX", "SampleY", "Rank_correlation")
  ggdata <- ggplot(cortable.reordered, aes(x = SampleX, y = SampleY))
  graph <- ggdata + geom_tile(aes(fill = Rank_correlation), color ="black") + 
    geom_text(aes(label=round(Rank_correlation, digits = 2))) + 
    scale_fill_gradientn(colors = c("steelblue","white") ,
                         values = scales::rescale(c(0.5,1))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  ggsave("corrleation_heatmap.pdf", width = width, height = height)
}


################################################################################

# Script to make a rank plot. Datatable must have gene names in column 1 and sgRNA data in column 2

screen.rankplot<-function(datatable) {
  filename <- colnames(datatable)[2]
  values <- datatable[order(datatable[,2]),]
  values$rank <- 1:length(values[,1])
  pdf(paste(filename, "_rankplot.pdf", sep = ""), 10,7)
  par(mar=c(5,6,4,1)+.1)
  plot(values$rank, values[,2], pch = 20, 
       main = paste("Rank Plot of Avg sgRNA Scores", filename),
       xlab="sgRNA enrichment rank", 
       ylab="sgRNA frequency")
  #(log2 normalized sgRNA counts)")
  abline(h = 0, col = 2, lty=2)
  dev.off()
}


################################################################################
#Extra graphing Functions
################################################################################

# updated version uses ggrepel and uses ngenes is used to define number of top genes if top is TRUE, 
# otherwise it can accept a vector containing a list of specific gene names desired to label. 
# blackLabels boolean sets label color to black, otherwise it is color of the point.

screen.scatterplot<-function(datatable, ngenes = 10, top = TRUE, blackLabels = FALSE, width = 8, height = 5) {
  datatable <- datatable[order(datatable[,3], decreasing = top),]
  datatable$Effect <- datatable[,3] - datatable[,2]
  datatable <- na.omit(datatable)
  if (top == TRUE && ngenes > 0) {
    labeltext <- rep("", length(datatable[,1]))
    labeltext[1:ngenes] <- datatable[1:ngenes,1]
  } else if (top == FALSE && ngenes > 0) {   #enters if wanting a custom label of genes
    labeltext <- rep("", length(datatable[,1]))
    genesofInterest <- match(ngenes, datatable[,1])
    labeltext[genesofInterest] <- datatable[genesofInterest,1]
  }
  if (blackLabels == FALSE && ngenes > 0) {
    ggdata <- ggplot(datatable, aes(x = datatable[,2], y = datatable[,3], color=Effect, label = labeltext)) + 
      geom_point(size = 1.5) + geom_text_repel(size = 5)
  } else if (blackLabels == TRUE && ngenes > 0) {
    ggdata <- ggplot(datatable, aes(x = datatable[,2], y = datatable[,3], color=Effect, label = labeltext)) + 
      geom_point(size = 1.5) + geom_text_repel(color = "black", box.padding = 0.31)
  } else {
    ggdata <- ggplot(datatable, aes(x = datatable[,2], y = datatable[,3], color=Effect)) + geom_point(size = 1.5)
  }
  graph <- ggdata + 
    scale_color_gradientn(colors = c("steelblue","grey","red"), 
                          values = scales::rescale(c( min(datatable$Effect), 0, max(datatable$Effect)))) + 
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.5) +
    ylab(paste(colnames(datatable)[3], "sgRNA Frequency", sep = " ")) + xlab(paste(colnames(datatable)[2], "sgRNA Frequency", sep = " ")) + 
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.7, linetype = "solid"), 
          axis.text = element_text(size=14), axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold")) +
        ylim(c(min(datatable[,3]), max(datatable[,3]))) #+   #for automatic lims
      #coord_cartesian(xlim = c(-4,2.5), ylim = c(-1.25,7.25))  #for manual lims   
  print(graph)
  ggsave("sgRNA_scatterplot.pdf", width = width, height = height)
}


################################################################################

# Script to make labeled rank plot from the drugZ output file

drugz.rankplot<-function(dataFile, genelabels = NULL, width = 6, height = 6) {
  datatable <- read.table(dataFile, sep = "\t", header = T, stringsAsFactors = F)
  datatable <- subset(datatable, select = c(1,4,6))
  datatable <- datatable[order(datatable$normZ),]
  datatable <- na.omit(datatable)
  labeltext <- rep("", length(datatable[,1]))
  if (is.null(genelabels) == FALSE) {
    generows <- data.frame(gene = genelabels, index = match(genelabels, datatable$GENE))
    generows <- na.omit(generows)
    labeltext[generows[,2]] <- as.character(generows[,1])
  }
  ggdata <- ggplot(datatable, aes(x = rank_synth, y = normZ, color=normZ))
  graph <- ggdata + geom_hline(yintercept = 0, color="black", linetype="dashed", size=0.5) +
    geom_point(size = 2) + 
    scale_color_gradientn(colors = c("steelblue","grey","red"), 
                          values = scales::rescale(c(min(datatable$normZ),-2.5, 0, 2.5,max(datatable$normZ)))) + 
    ylab("Normalized Z score") + xlab("sgRNA Rank") + theme_bw() + 
    theme(axis.text = element_text(size=18), axis.title = element_text(size=18, face="bold"), legend.position = "none") +
    #coord_cartesian(xlim = c(-4,2.5), ylim = c(-1.25,7.25))   #for manual lims
    ylim(c(min(datatable$normZ), max(datatable$normZ))) + #for automatic lims
    scale_y_continuous(breaks = seq(-5,20,5), labels = seq(-5,20,5)) + 
    geom_text_repel(label = labeltext, box.padding = 1, min.segment.length = 0.1, size = 7)
  ggsave("drugZ_rankplot.pdf", width = width, height = height)
}

################################################################################

#Script to perform PCA on screen results. Will make a 2D graph.
#needs the screen datatable file with format from the screen analysis
#also will need a df defining sample names and conditions in the order of the columns
#in the data table minus the gene column with headers Sample, Condition

screen.pca<-function(datatable, conditions) {
  datatable <- read.table(datatable, sep = "\t", header = TRUE,check.names = FALSE)
  pcgroups <- read.table(conditions, sep = "\t", header = TRUE)
  datatable <- na.omit(datatable)
  #adds a log transformation of the data but not needed if a logFC already
  screenmatrix <- log(datatable[,-1])
  screenmatrix <- screenmatrix[is.finite(rowSums(screenmatrix)),]
  
  screenmatrix <- t(screenmatrix)
  pca.input <- screenmatrix[,apply(screenmatrix,2,var,na.rm=TRUE) != 0]
  pca.output <- prcomp(pca.input) # computes variance
  ggbiplot(pca.output, cex.lab=3, scale = .75, choices = c(1,2), ellipse = F,
           groups = pcgroups[,2], labels.size = 2,
           var.axes = 'F',
           labels = pcgroups[,1]) + 
    geom_point(aes(colour= pcgroups[,2]), size=7, pch=19, alpha=0.5)+
    theme(title = element_text(size = 20),axis.title = element_text(size=20),axis.text = element_text(size=20))+
    ggtitle('Unsupervised PCA', subtitle = "Using Genes from Panel")
  ggsave("PCAnalysis.pdf", width = 12, height = 8)
}

