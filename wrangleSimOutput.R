require(dplyr)
#combine results from the simulations to get single data.frames of gelman and 
resultsDir<-"simTest2"
files<-list.files(resultsDir)
resultsFiles<-files[grep("res",files)]
gelmanFiles<-files[grep("gelman",files)]
rm(files)

getGelman<-function(files){
  output<-list()
  for(f in files){
    fileNum<-which(f==files)
    result<-readRDS(file.path(resultsDir,f)) %>% data.frame()
  
    result$nSites<-strsplit(f,"nSites")[[1]][2] %>%
                 strsplit("nYears") %>% .[[1]] %>% .[1] %>%
                 as.numeric()
  
    result$nYears<-strsplit(f,"nYears")[[1]][2] %>%
                 strsplit("r") %>% .[[1]] %>% .[1] %>%
                 as.numeric()
  
  result$trueTrend<-strsplit(f,"nYears")[[1]][2] %>%
                    strsplit(".rds") %>% .[[1]] %>% .[1] %>%
                    strsplit("r") %>% .[[1]] %>% .[2] %>%
                    as.numeric()
  output[[fileNum]]<-result
  }
  output<-bind_rows(output) %>% data.frame()
  return(output)
}

gelman<-getGelman(gelmanFiles)
gelman[which(apply(gelman[,1:9],1,function(x){any(x>1.1)})),c("nSites","nYears","trueTrend")]

