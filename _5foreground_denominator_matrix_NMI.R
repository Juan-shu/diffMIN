load("./GSE95140.RData")
load(file = "./fgnet_inter1.RData")
gene_count <- ncol(dat)

# Compute non-perturb data NMI
library(parallel)
cl <- makeCluster(40)
fg.NMI.down <- parLapply(cl,fg.denominator.list,function(x){
  NMI <- function(dis){
    library(infotheo)
    library(reshape2)
    library(dplyr)
    # Compute mutual information 
    var_id <- colnames(dis)
    dis <- data.matrix(dis)
    res <- .Call("buildMIM", dis, NROW(dis), NCOL(dis), 0, PACKAGE = "infotheo")
    dim(res) <- c(NCOL(dis), NCOL(dis))
    res <- as.matrix(res)
    rownames(res) <- var_id
    colnames(res) <- var_id
    #Compute entrophy
    library(entropy)
    H_gene <- apply(dis,2,function(x){table(x) %>% entropy})
    Entrophy <- outer(H_gene,H_gene,FUN = function(x,y){(x+y)/2})
    NMI_value <- melt(res/(Entrophy+.Machine$double.eps))[melt(upper.tri(res/Entrophy))[,3],]
    name <- paste0(NMI_value[,1],"_",NMI_value[,2])
    NMI_value <- NMI_value[,3]
    # Eliminates floating point arithmetic errors  
    NMI_value[NMI_value < 0 ] <- 0
    # Add a minimum value to prevent the denominator from being 0
    NMI_value <- NMI_value+.Machine$double.eps
    names(NMI_value) <- name
    return(NMI_value)
  }
  result <- lapply(x,NMI)
  result <- do.call(cbind,result)
  colnames(result) <- c("Cell",paste0("Cell_random",1:100))
  return(result)
})

save(fg.NMI.down,file = "./fgnet_inter3.RData")
