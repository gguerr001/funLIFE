
library(tidyr)

##I = 2
args = commandArgs(trailingOnly=TRUE)

out_mcl = read.csv(args[1], header = F)
out_mcl$V2 = seq(1:nrow(out_mcl))

out_mcl = separate_rows(out_mcl,V1,sep = "\t", convert = T)

colnames(out_mcl) = c('#BGC Name', 'Family Number')

write.table(out_mcl, args[2], sep = '\t', row.names = F, quote = F)
