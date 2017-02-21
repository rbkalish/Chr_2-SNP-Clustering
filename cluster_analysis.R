library(igraph)
covariance_calc<-function(..., fileList = NULL, fileType = "gml"){
  
  #addition of individual graph file and file list
  fileList = c(..., fileList)
  
  #declaration of parameter lists
  graph = c()
  top_betweenness = c()
  top_degree = c()
  top_close = c()
  top_edge = c()
  top_dm = c()
  
  #calc uses igraph functions to return data frame of covariance parameters
  for (file in fileList){
    if (grepl(fileType, file)){
      graph <- c(graph, file)
      g1<-read_graph(file, format = fileType)
      btw<-betweenness(g1)
      top_betweenness <- c(top_betweenness, sum(btw)/length(btw))
      deg <- degree(g1)
      top_degree <- c(top_degree, sum(deg)/length(deg))
      close <- closeness(g1)
      top_close <- c(top_close, sum(close)/length(close))
      edge<-edge.betweenness(g1)
      top_edge <- c(top_edge, sum(edge)/length(edge))
      top_dm <- c(top_dm, diameter(g1))
    }
  } 
  
  #Creation of data.frame
  estimate <- data.frame(graph)
  estimate$betweenness <- top_betweenness
  estimate$degree <- top_degree
  estimate$closeness <- top_close
  estimate$edge_betweenness <- top_edge
  estimate$diameter <- top_dm
  estimate
}


#read in estimate of covariance
data <- covariance_calc(list.files())
data$graph <-lapply(data$graph, function(x){as.integer(substr(x, 13,16))})



#scale data and find euclidean distance
data.stand <- scale(data[-1])
distance <- dist(data.stand, method = "euclidean")

#cluster using wards method and squared eclidean distance
H.fit <- hclust(distance^2, method="ward.D")

plot(H.fit)


# define some clusters
cluster <- cutree(H.fit, h = 1000)

#4 is the smallest cluster, will examine cluster further
table(cluster)

data[["cluster"]]<-(cluster)

cluster5 <-subset(data, cluster == "5")

#add clusters to estimate csv file
data <- data.frame(lapply(data, as.character), stringsAsFactors=FALSE)
write.csv(data,"estimateR.csv")


#analyzation of SNP reference data frame

load("C:/Users/Robert/Documents/topologies/RSID files/snp_id.rda")
snp_id$Graph<-lapply(snp_id$Graph, function(x){as.integer(substr(x, 13,16))})

#returns data frame of SNPs that are included in cluster5
snp_5 <-subset(snp_id, Graph %in% cluster5$graph)



