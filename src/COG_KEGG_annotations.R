library(data.table)
args = commandArgs(trailingOnly=TRUE)
##arg[1] is the annotation file from emapper
##arg[2] is the cog descriptions
##arg[3] is the kegg descriptions
##arg[4] is the cog output
##arg[5] is the kegg output

annotation <- read.delim2(args[1], header = F, comment.char = '#')
annotation <-annotation[,c(1,5,7,8,12)]
colnames(annotation) <- c('id', 'cogid', 'cog_category', 'cog_description', 'keggid')



#Extract COGs
cog_annotation <- annotation[,c(1,2,3,4)]
cog_annotation$cogid =  sapply(strsplit(cog_annotation$cogid,"@"), `[`, 1)
cog_annotation <- data.table(cog_annotation)
cog_annotation[cog_annotation$cog_category == '-']$cog_category = 'S'
cog_annotation$cog_category <- substr(cog_annotation$cog_category, 1, 1)

cog_categories = cog_annotation[,2:3]
cog_categories = unique(cog_categories)
COG_general= read.csv(args[2], header = F)
COG_general = COG_general[,c(2,7)]
COG_general = unique(COG_general)
COG_general = COG_general[!is.na(COG_general$V7),]

cog_categories = merge(cog_categories, COG_general, by.x = 'cog_category', by.y = 'V2')
cog_categories = cog_categories[,c(2,1,3)]
new_row = data.frame(cogid = 'Unknown', cog_category = 'S', V7 = 'Function unknown')
cog_categories = rbind(cog_categories, new_row)
write.table(cog_categories, 'intermediate_files/cog_annotation/cog_custom_groups.txt', row.names = F, col.names = F, sep= '\t')


cog_annotation = cog_annotation[,c(1,2,4)]
colnames(cog_annotation)[1] = 'gene'
write.table(cog_annotation, args[4], row.names = F)


#Extract KEGGs
kegg_annotation <- annotation[,c(1,5)]
keep <- grep("ko", kegg_annotation$keggid)
kegg_annotation <-kegg_annotation[keep,]
kegg_annotation$keggid <- sapply(strsplit(kegg_annotation$keggid,","), `[`, 1)
kegg_annotation$keggid <- sapply(strsplit(kegg_annotation$keggid,":"), `[`, 2)

##KO descriptions
ko_descriptions <- read.delim2(args[3], header = F)
ko_descriptions$V1 <- sapply(strsplit(ko_descriptions$V1,":"), `[`, 2)
colnames(ko_descriptions)<- c('keggid', 'kegg_description')

kegg_annotation <- merge(kegg_annotation, ko_descriptions, by= 'keggid', all.x = T)
kegg_annotation <-kegg_annotation[,c(2,1,3)]
colnames(kegg_annotation)[1]<- 'gene'

write.table(kegg_annotation, args[5], row.names = F)



