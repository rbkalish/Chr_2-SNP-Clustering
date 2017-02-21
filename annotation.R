matchseq<-function(vcf, db = TxDb.Hsapiens.UCSC.hg19.knownGene, noncoding = TRUE){
  intersect(seqlevels(vcf), seqlevels(db))
  
  ## Rename seqlevels with renameSeqlevesl().
  vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf))) 
  
  ## Confirm.
  intersect(seqlevels(vcf), seqlevels(db))
  
  ## Overlaps for all possible variant locations.
  loc <- locateVariants(vcf, db, AllVariants())
  names(loc) <- NULL
  out <- as.data.frame(loc)
  out$names <- names(vcf)[ out$QUERYID ]
  out <- out[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "PRECEDEID", "FOLLOWID")]
  out <- unique(out)
  
  #If noncoding is false, the dataframe is subsetted to remove noncoding SNPS and the preceding and following gene ids
  if (noncoding == FALSE){
    out<- subset(out, select = -c(PRECEDEID, FOLLOWID))
    out <-subset(out, !is.na(out$GENEID))
  }
  out
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

getmode <- function(v) {
    uniqv <- unique(v[!is.na(v)])
    counter <- lapply(uniqv, function(id){length(grep(id, v))})
    do.call(rbind, Map(data.frame,ID=uniqv,COUNT=counter))
}

#Creation of RSID files for NCBI batch Query
rsid<-snp_id$SNPs
rssplit<-split(rsid, ceiling(seq_along(rsid)/20000))
length(rssplit$`1`)
length(rssplit$`5`)
write(rssplit$`1`, "rsid1") #rsid files put through dbSNP batch query to obtain VCF
write(rssplit$`2`, "rsid2")
write(rssplit$`3`, "rsid3")
write(rssplit$`4`, "rsid4")
write(rssplit$`5`, "rsid5")

#Batch Query returns vcf files which are reinputted
vcf1 <- readVcf("VCF\\VCF1.vcf", genome="hg19")
vcf2<- readVcf("VCF\\VCF2.vcf", genome="hg19")
vcf3 <- readVcf("VCF\\VCF3.vcf", genome="hg19")
vcf4 <- readVcf("VCF\\VCF4.vcf", genome="hg19")
vcf5 <- readVcf("VCF\\VCF5.vcf", genome="hg19")



#files are read into a list and put into matchseq function
#matchseq outputs are bound together and merged with original SNP id dataframe
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
colnames(snp_id) <- c("graph", "names")
files <- list(vcf1, vcf2, vcf3, vcf4, vcf5)
out1<-do.call(rbind, lapply(files, function(x){matchseq(x, noncoding = FALSE)}))
output<-merge(snp_id, out1, by= "names")


library(KEGGREST)
library(plyr)     

# subsetting genes of SNPS in cluster 5
#Cluster 5 isolated as outlier in cluster_analysis.R

`%!in%` = Negate(`%in%`)
snp5_ann<-subset(output, names %in% snp_5$SNPs)
notsnp5_ann <- subset(output, names %!in% snp_5$SNPs)


#Submission of cluster4 vs non-cluster4 to snp pathway

snp5_ann$PATHWAY<-unname(get_kegginfo(snp5_ann$GENEID))
notsnp5_ann$PATHWAY<-get_kegginfo(notsnp5_ann$GENEID)

get_kegginfo<-function(..., id_list = NULL, org_code = "hsa"){
  `%!in%` = Negate(`%in%`)
  
  #adds indvidual gene ids and input id_list together to make new list for submission
  id_list = c(..., id_list)
  
  #attaches org_code to ids so they are ready to submit
  id_list<-lapply(id_list, function(x){paste(org_code, x, sep = ":")})
  info = list()
  info2 = c()
  
  
  for (id in id_list){
    if (id %!in% info2){
      
      #if id has not already been submitted, program will try to fetch results from kegg
      info1<-tryCatch(
        {  keggGet(id) 
        }, error = function(cond) {
          message(paste(cond, "ID: ", id)) #error message with id will return for each gene not found
          return(NULL) #program will continue with null input
        }
      )
    }
    
    # if id has already been submitted, index will be found and comparaed to find pathway
    else {  
      ind<-match(id, info2)
      info <- append(info, info[ind])
      info2 <- append(info2, id)
      next
    }
    
    #if pathway exists, it will be added to pathway list otherwise NA will be appended
    if (!is.null(info1[[1]]$PATHWAY)){
      info <- append(info, info1[[1]]$PATHWAY[1])
      info2 <- append(info2, id)
    }
    else{
      info<- append(info, NA)
      info2 <- append(info2, id)
    }
  }
  #returns list of pathways in same order as gene ID list inputed
  info
}


# counts of pathways sorted by cluster

#cluster 5
snp5_pathway<-getmode(unname(snp5_ann$PATHWAY))
snp5_pathway$CLUSTER <- "Outliers"
snp5_pathway<-subset(snp5_pathway, COUNT > 30)


#other clusters
notsnp5_pathway<-getmode(notsnp5_ann$PATHWAY)
notsnp5_pathway$CLUSTER <- "Normal Distribution"
notsnp5_pathway <- subset(notsnp5_pathway, COUNT > 400)


#Preparing table for plotting
pathway_count <- rbind(snp5_pathway,notsnp5_pathway)
pathway_count$ID <- factor(pathway_count$ID, levels = pathway_count$ID[order(pathway_count$COUNT)])


#Plotting of relationships
library(ggplot2)
ggplot(pathway_count, aes(x = CLUSTER, y = COUNT, fill= ID)) + 
  geom_bar(position = "fill",stat = "identity") + 
  labs(title="Pathway Distribution of Outliers Relative to Normal",y = "PROPORTION", x = NULL)



snp5_ann <- data.frame(lapply(snp5_ann, as.character), stringsAsFactors=FALSE)
write.csv(pathway_count[order(pathway_count$COUNT),], "PathwayDistribution.csv")


