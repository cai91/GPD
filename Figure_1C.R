# load libraries
library(UpSetR)
library(ggplot2)

# load files
clst.data = read.delim("matrix_virome_DBs.txt", row.names=1, check.names = FALSE)
species.multiple = rownames(clst.data)[which(rowSums(clst.data) > 1)]
clst.data[clst.data > 0] = 1

# plot
upset(clst.data, sets.x.label = "Total VCs",text.scale = c(2.5,2.2, 2.5, 1.5, 2.0, 2.3),queries = list(list(query = intersects, params = list("GPD"), color = "darkred", active = T)),nintersects = 10, sets = colnames(clst.data), order.by="freq",  
      sets.bar.color = rev(c("black", "black","black","darkred")))