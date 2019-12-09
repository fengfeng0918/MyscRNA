library(scran)
# setwd("C:\\Users\\fdong\\OneDrive\\python_remote\\single_cell_seq")

data_mat=read.csv("X.txt",header = T,sep ="\t")
obs = read.csv("obs.txt",header = T,sep ="\t")
data_mat = t(data_mat)
input_groups = obs$groups
size_factors = computeSumFactors(data_mat, clusters=input_groups, min.mean=0.1)
write.csv(size_factors,file="size_factors.txt",row.names=F, col.names = FALSE)