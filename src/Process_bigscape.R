library(plyr)
library(dplyr)
library(tidyverse)
library(reshape2)
library(data.table)

args = commandArgs(trailingOnly=TRUE)


##Create binary table BiG-SCAPE for fungilife



network_annotations = read.delim(args[2], header = T)
network_annotations$origin = substr(network_annotations$BGC, 1, 3)
network_annotations = network_annotations[network_annotations$origin != 'BGC',]
network_annotations$origin = NULL


clusters = read.delim(args[1], header = T)
clusters$origin = substr(clusters[,1], 1, 3)
clusters = clusters[clusters$origin != 'BGC',]
clusters$origin = NULL

clusters$genome = sapply(strsplit(clusters[,1] ,"_"), `[`, 1)


# Generate abs/prese table

clusters$count = 1
df = dcast(clusters, Family.Number~genome,
      value.var='count', sum, fill= 0)
df$Family.Number = paste0('GCF', df$Family.Number)
rownames(df) = df$Family.Number
df$Family.Number = NULL


##Fix column names

names_equivalence = read.table(args[6], header = T)
names_equivalence$id = sapply(strsplit(names_equivalence[,2] ,"_"), `[`, 3)

names_equivalence[,2] = str_remove(names_equivalence[,2], '_O.fna')

M= data.frame(old = colnames(df))
M = merge(M, names_equivalence, by.x = 'old', by.y = 'id')
colnames(df) = M[,3]
write.table(df, args[4], row.names = T)

#Generate annotation file

annotation = data.frame(bgc = clusters$X.BGC.Name ,GCF_No = clusters$Family.Number)
annotation = merge(annotation, network_annotations, by.x = 'bgc', by.y = 'BGC')
annotation2save = annotation[,c(2,6)]
colnames(annotation2save) = c("GCF No","BGC Class")
write.csv(annotation2save, args[3], row.names = F)



#Generate BGC_description files

network_annotation2save = annotation
network_annotation2save$GCF_No <- paste0('GCF', network_annotation2save$GCF_No)
network_annotation2save = network_annotation2save[,c('GCF_No',	'bgc',	'Description',	'Product.Prediction',	'BiG.SCAPE.class',	'Organism')]
colnames(network_annotation2save)[1:2] = c('GCF', 'Accession.ID')
write.table(network_annotation2save, args[5], col.names = T, row.names = F, quote = F, sep = '\t')
