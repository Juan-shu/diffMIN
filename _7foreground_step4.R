# Set compute parameter
class <- "g1"

# Loading package and data
load("./fgnet_inter4.RData")

# Add node information
library(dplyr)
node_node <- rownames(net.list[[1]])
node_node <- strsplit(node_node,"_") %>% unlist() %>% matrix(nrow = 2) %>% t() %>% as.data.frame()
colnames(node_node) <- c("nodeA","nodeB")
node_node$edge <- paste0(node_node$nodeA,"_",node_node$nodeB)

# Data pre-process
for (i in 1:length(net.list)) {
  colnames(net.list[[i]]) <- c(paste0("net",i),paste0("ran",1:100))
  net <- net.list[[i]]
  net <- net[apply(net,1,function(x){x[1] > quantile(x[-1],0.9) | x[1] < quantile(x[-1],0.1)}),]
  net.list[i] <- list(net)
  print(i)
}

# Compute p value
library(parallel)
cl <- makeCluster(40)
result_net <- parLapply(cl,net.list,function(x){
  # Define permutation test function
  permutation_test <- function(x, y, B = 1000) {
    set.seed(1)
    orig_diff <- abs(mean(x) - mean(y))
    diffs <- replicate(B, {
      perm <- sample(c(x, y), length(x) + length(y), replace = FALSE)
      abs(mean(perm[1:length(x)]) - mean(perm[(length(x)+1):(length(x)+length(y))]))
    })
    pval <- sum(diffs >= orig_diff) / B
    return(pval)
  }
  result <- apply(x,1,function(y){
    a <- y[1]
    b <- y[-1]
    result <- permutation_test(a,b)
    return(result)
  })
  x <- x[result < 0.05,]
  return(x)
})

# Abstract net
all.net <- data.frame(edge = node_node$edge)
for (i in 1:length(result_net)) {
  net <- result_net[[i]]
  net <- as.data.frame(net)
  net <- net[node_node$edge,]
  all.net <- cbind(all.net,net[,1])
}
all.net <- all.net[,-1]
colnames(all.net) <- c(paste0("fg_cell",1:length(result_net)))
rownames(all.net) <- node_node$edge

save(all.net,file = paste0("./GSE95140_fg_",class,".RData"))
