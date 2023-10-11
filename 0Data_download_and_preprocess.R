### Data download and preprocess
### Dataset GSE95140
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95140/suppl/GSE95140_single-cardiomyocyte_RNA-seq.txt.gz",
              destfile = "./GSE95140.txt.gz")
# Load data
exp <- read.table(file = "./GSE95140.txt.gz",header = T)
rownames(exp) <- exp$geneName
exp <- exp[,-1]
exp <- exp[-22135,]

# Remove bulk sample
library(stringr)
exp <- exp[,!(str_detect(colnames(exp),"Bulk"))]

# Remove low expressed samples and probes
filter.matrix <- exp > 0
exp <- exp[rowSums(filter.matrix) > 0,]
filter.matrix <- exp > 0.1
exp <- exp[,colSums(filter.matrix) > 5000]

# Select sub-data to compute network
gene_count <- rowSums(exp > 0)
select.data <- exp[names(gene_count[gene_count > 20]),] # gene count
var.select <- apply(select.data, 1, var) # gene varible
select.gene <- names(sort(var.select,decreasing = T))
fd <- select.data[select.gene[1:400],] #larger matrix means longer time, we select 400 genes in this example

# Discrete
library(infotheo)
dat <- discretize(t(log2(fd+1)))
rownames(dat) <- colnames(fd)

# Group information
library(dplyr)
group <- rep(c("sham","d3","w1","w2","w4","w8"),c(64,58,82,61,73,58))

# Save data
save(group,fd,dat,file = "./GSE95140.RData")
