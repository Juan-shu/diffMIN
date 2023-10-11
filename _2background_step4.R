# Loading data
load(file = "./bg1.RData")

# Add node information
library(dplyr)
node_node <- rownames(bg.net[[1]])
node_node <- strsplit(node_node,"_") %>% unlist() %>% matrix(nrow = 2) %>% t() %>% as.data.frame()
colnames(node_node) <- c("nodeA","nodeB")
node_node$edge <- paste0(node_node$nodeA,"_",node_node$nodeB)

# Data preprocess
# The processing here is to reduce the amount of computation,
# and the removed part almost always has a p value greater than 0.05 in the replacement test
library(parallel)
cl <- makeCluster(40)

processed.bg.net <- parLapply(cl,bg.net,function(x){
  x <- x[apply(x,1,function(y){y[1] > quantile(y[-1],0.9) |y[1] < quantile(y[-1],0.1)}),]
})
stopCluster(cl)

# Compute p value
cl <- makeCluster(40)
result_net <- parLapply(cl,processed.bg.net,function(x){
  # Permutation test function
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

# Save result
for (i in 1:length(result_net)) {
  sample_net <- data.frame(result_net[[i]])
  result_net[i] <- list(sample_net[node_node$edge,])
}

net <- data.frame(bg1 = result_net[[1]][,1],row.names = node_node$edge)
for (i in 2:length(result_net)) {
  net <- cbind(net,result_net[[i]][,1])
}
colnames(net) <- paste0("bg",1:length(result_net))


save(net,node_node,file = paste0("./GSE95140_bg.RData"))
