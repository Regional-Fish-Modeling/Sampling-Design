outDir<-"~/output"
files<-list.files(outDir)

for(f in files){
  assign(f,readRDS(file.path(outDir,f)))
}

results<-do.call(rbind,mget(files))
rownames(results)<-NULL
saveRDS(results,"~/trendSimResults.rds")