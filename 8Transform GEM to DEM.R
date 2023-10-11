# Integrate all the result
library(stringr)
library(dplyr)
EM_to_DM <- function(EM,fd){
  library(dplyr)
  gene <- rownames(fd)
  dim1 <- dim(fd)[1]
  dim2 <- dim(fd)[2]
  matrix0 <- matrix(0,dim1,dim1)
  colnames(matrix0) <- rownames(fd)
  rownames(matrix0) <- rownames(fd)
  degree.matrix <- data.frame(gene = rownames(matrix0))
  for (i in 1:dim2) {
    edge <- EM[,i]
    names(edge) <- rownames(EM)
    edge <- na.omit(edge)
    edge <- names(edge)
    node <- strsplit(edge,"_") %>% unlist() %>% matrix(nrow = 2) %>% t()
    node <- data.frame(V1 = match(node[,1],rownames(matrix0)),V2 = match(node[,2],rownames(matrix0)))
    dir1 <- node$V1+(node$V2-1)*dim1
    dir2 <- node$V2+(node$V1-1)*dim1
    sample <- matrix0
    sample[c(dir1,dir2)] <- 1
    degree.matrix <- cbind(degree.matrix,apply(sample,1,sum))
    cat(paste0(i,"\n"))
  }
  degree.matrix <- degree.matrix[,-1]
  DM <- degree.matrix
  colnames(DM) <- colnames(fd)
  return(DM)
}


files <- list.files()[str_detect(list.files(),"GSE95140_fg_g")]
load("./GSE95140_bg.RData")
for (i in files) {
  load(paste0("./",i))
  all.net <- all.net[rownames(net),]
  net <- cbind(net,all.net)
  print(i)
}

load("./GSE95140.RData")
colnames(net) <- colnames(fd)
DM <- EM_to_DM(net,fd)
