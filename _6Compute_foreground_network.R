# Loading package and data
load("./fgnet_inter2.RData")
load("./fgnet_inter3.RData")

# Compute foreground net
net.list <- list()
for(i in 1:length(fg.NMI.up)){
  net <- fg.NMI.up[[i]]/fg.NMI.down[[i]]
  net.list[i] <- list(net)
  print(i)
}

save(net.list,file = "./fgnet_inter4.RData")
