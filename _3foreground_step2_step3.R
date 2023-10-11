class <- "g1"

# Load packages and data
load("./GSE95140.RData")
library(dplyr)

# Divide foreground data into some small groups
con_length <- sum(group == group[1])
all_length <- length(group)
group <- rep(c("Con",paste0("g",1:ceiling((all_length-con_length)/50))),
             c(con_length,rep(50,ceiling((all_length-con_length)/50) -1),(all_length-(ceiling((all_length-con_length)/50) -1) * 50 -con_length)))
cat("Computing",class,"...\n")
group.list <- list()
fg.dat <- dat[group == "Con" | group == class,]


# Construct perturbed matrix (PM) and random perturbed matrix (RPM)
# Construct denominator matrix (DeM) and random denominator matrix (RDeM)
fg.perturb.list <- list()
fg.denominator.list <- list()
gene_number <- ncol(fg.dat)

fg.perturb.list <- list()
fg.denominator.list <- list()
for (i in (con_length+1):nrow(fg.dat)) {
  fg.ran <- fg.dat[c(1:con_length,i),]
  for (j in 1:100) {
    fg.ran.dat <- apply(fg.dat[c(1:con_length,i),],2,function(x){x = sample(x,length(x))}) %>% as.data.frame()
    fg.ran <- cbind(fg.ran,fg.ran.dat)
  }
  bg.ran <- fg.ran[1:con_length,]
  
  # Changing the format of the data store to facilitate multithreading
  bg.ran.list <- list()
  fg.ran.list <- list()
  for (j in 1:101) {
    bg.ran.list[j] <- list(bg.ran[,(1+gene_number*j-gene_number):(gene_number*j)])
    fg.ran.list[j] <- list(fg.ran[,(1+gene_number*j-gene_number):(gene_number*j)])
  }
  fg.perturb.list[i-(con_length)] <- list(fg.ran.list)
  fg.denominator.list[i-(con_length)] <- list(bg.ran.list)
}

# Saving result
save(fg.perturb.list,fg.denominator.list,file = "./fgnet_inter1.RData")
