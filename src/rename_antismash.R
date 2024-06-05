
library(stringr)

folders = list.files('intermediate_files/antismash/')
ids = sapply(strsplit(folders ,"_"), `[`, 3)
dir.create('intermediate_files/antismash_renamed')
for (i in 1:length(folders)){
  folder = folders[i]
  gbks = list.files( paste0('intermediate_files/antismash/', folder, '/'  ), pattern = '.gbk'  )
  gbks = gbks[-c(1)]
  if (length(gbks) == 0){
    next
  }else{
    gbks_renamed = paste0( ids[i], '_', seq(1, length(gbks)), '.' ,sapply(strsplit(str_remove( gbks, '.gbk'), ".", fixed=TRUE), tail, 1), '.gbk'   )
    dir.create(paste0('intermediate_files/antismash_renamed/', folder))
    file.copy(paste0('intermediate_files/antismash/', folder, '/', gbks  ),  paste0('intermediate_files/antismash_renamed/', folder)  )
    file.rename( paste0('intermediate_files/antismash_renamed/', folder, '/', gbks), paste0('intermediate_files/antismash_renamed/', folder, '/', gbks_renamed) )
    
  }
  
  
}

checkpoint = data.frame(Checkpoint= 'antismash', value = 'Yes')
write.table(checkpoint, 'intermediate_files/antismash_checkpoint.txt')
