# The algorithm is computationally intensive, and we recommend using multithreaded operations
# 48 CPU and 200G memory is enough in this example
# foreground network are recommended to compute divided into some batches so we recommend running the work in the background of the system

# Loading data and library packages
load("./GSE95140.RData")
library(dplyr)
library(parallel)

# Define NMI function
NMI <- function(dis){
  # library packages
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
  #Compute entrophy and normalized mutual information
  library(entropy)
  H_gene <- apply(dis,2,function(x){table(x) %>% entropy})
  Entrophy <- outer(H_gene,H_gene,FUN = function(x,y){(x+y)/2})
  NMI_value <- melt(res/(Entrophy+.Machine$double.eps))[melt(upper.tri(res/Entrophy))[,3],]
  name <- paste0(NMI_value[,1],"_",NMI_value[,2])
  NMI_value <- NMI_value[,3]
  NMI_value[NMI_value < 0 ] <- 0# Eliminates floating point arithmetic errors
  names(NMI_value) <- name
  return(NMI_value)
}



### 1 Construct perturbed matrix (PM) and random perturbed matrix (RPM)

# In this example, group sham is considered as a background group
dat.bg <- dat[1:sum(group == group[1]),]

# 100 RPMs are built, More quantity means more computation
for (i in 1:100) { 
  disorder_bg_exp <- apply(dat[1:sum(group == group[1]),],2,
                           function(x){x = sample(x,length(x))}) %>% as.data.frame()
  dat.bg <- cbind(dat.bg,disorder_bg_exp)
  cat("Random background", i, "has been constructed\n")
}

# Changing the format of the data store to facilitate multithreading
gene_number <- ncol(dat.bg)/101
dat.bg.perturb.list <- list()
for (i in 1:101) {
  dat.bg.perturb.list[i] <- list(dat.bg[,(1+gene_number*i-gene_number):(gene_number*i)])
}


### 2 Construct denominator matrix (DeM) and random denominator matrix (RDeM)
dat.bg.denominator.list <- list()
for (i in 1:nrow(dat.bg)) {
  def <- dat.bg[-i,]
  def.list <- list()
  for (j in 1:101){
    def.list[j] <- list(def[,(1+gene_number*j-gene_number):(gene_number*j)])
  }
  dat.bg.denominator.list[i] <- list(def.list)
}


### 3 Compute NMI of PMs and RPMs
cl <- makeCluster(40)
dat.bg.NMI <- parLapply(cl,dat.bg.perturb.list,NMI)
dat.bg.NMI <- do.call(cbind,dat.bg.NMI)
colnames(dat.bg.NMI) <- c("Cell",paste0("Cell_random",1:100))

### 4 Compute NMI of DeMs and RDeM
dat.bg.denominator.NMI <- parLapply(cl,dat.bg.denominator.list,function(x){
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

### 5 Compute delta NMI
bg.net <- lapply(dat.bg.denominator.NMI, function(x){
  x <- dat.bg.NMI/x
})

### Save result
save(bg.net,file = "./bg1.RData")
