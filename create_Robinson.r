#! /usr/bin/Rscript

library('seriation')

args = commandArgs(trailingOnly=TRUE)

n = as.numeric(args[1])
m <- random.robinson(n, anti=FALSE)
write.table(m, file="random_robinson_matrix.csv", sep=",",row.names=FALSE, col.names=FALSE)
