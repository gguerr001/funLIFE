
library(ape)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

genes = read.FASTA(args[1], type = "AA")

df = data.frame(augustus_name = names(genes), microLife_name = str_remove(args[1], '_O.aa') )
df$microLife_name = paste0(df$microLife_name, '|', sprintf('%0.6d', 1:nrow(df)))
df$microLife_name = sapply(strsplit(df$microLife_name ,"/"), `[`, 4)
#rename genes
names(genes) = df$microLife_name
write.FASTA(genes, args[2])

#write tags
write.table(df, args[3], row.names = F, quote = F, sep = '\t')

