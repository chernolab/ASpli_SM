
library(ASpli)

### Load simulation data
load("SimuData/AT1_ES-IR-Alt3ss-Alt5ss.Rdata")

### Running ASpli ####
contrast         <- c(1, -1) #Two conditions, A-B
bin_genome       <- TRUE
summarize_counts <- TRUE
find_du_bins     <- TRUE

#Bams for analysis 
bamPath     <- "BAMs/" 
bamFiles    <- c("A1.s.bam",
                 "A2.s.bam",
                 "A3.s.bam",
                 "B1.s.bam",
                 "B2.s.bam",
                 "B3.s.bam")

targets <- data.frame(bams = paste0(bamPath, bamFiles), 
                      phenotype = rep(c("A", "B"), each=3), 
                      stringsAsFactors = FALSE)

mergedBams <- data.frame(bams = paste0(bamPath, c("A.bam", "B.bam")), 
                         condition = c("A", "B"), stringsAsFactors = F)

#bin genome?
if(bin_genome){
  require(GenomicFeatures)
  if(FALSE){
   txdb        <- makeTxDbFromGFF(file = "BAMs/genes.gtf")
   saveDb(txdb,file="BAMs/tair10.sqlite")
  }else{
   txdb<-loadDb("BAMs/tair10.sqlite")  
  } 
  features    <- binGenome(txdb)
  save(features, file="ASpliRun/features.RData")
}else{
  (load("ASpliRun/features.RData")) 
}

#summarize bin and junction counts
if(summarize_counts){
  counts      <- gbCounts(features = features, targets = targets, minReadLength = 100, maxISize = 5000)
  asd         <- jCounts(counts = counts, features = features, minReadLength = 100)
  save(counts, asd, file="ASpliRun/counts.RData")
}else{
  (load("ASpliRun/counts.RData"))
}


#Find differentialy used bins
if(find_du_bins){
  j           <- jDUreport(asd, contrast = contrast, filterWithContrasted = T, 
                           mergedBams = mergedBams$bams, strongFilter = T, runUniformityTest = T)
  gb          <- gbDUreport(counts, contrast = contrast)
  sr          <- splicingReport(gb, j, counts)

  #integrateSingals to export
  is          <- integrateSignals(sr, asd )
  
  #export integrated signal results
  exportIntegratedSignals(is,sr=sr,counts = counts, features = features, asd = asd,
                          mergedBams = mergedBams, output.dir = "ASpliRun") 
  save(j, gb, sr, is,file="ASpliRun/du.RData")
}
