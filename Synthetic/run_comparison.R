library(ASpli)
library(readr) #To read Leafcutter, Majiq, rMats 
library(GenomicRanges)

## Functions####
#Function for venn diagrams with up to 4 elements
limma_venn <- function(item1, item2, item3 = NULL, item4 = NULL, names, main){
  item1 <- unique(item1)
  item2 <- unique(item2)
  
  items <- union(item1, item2)
  if(!is.null(item3)){
    item3 <- unique(item3)
    items <- union(items, item3)
  }
  if(!is.null(item4)){
    item4 <- unique(item4)
    items <- union(items, item4)
  }  
  itemsm <- matrix(0, nrow=length(items), ncol=ifelse(is.null(item4), ifelse(is.null(item3), 2, 3), 4))
  rownames(itemsm) <- items
  itemsm[item1, 1] <- 1
  itemsm[item2, 2] <- 1
  if(!is.null(item3)){
    itemsm[item3, 3] <- 1
  }
  if(!is.null(item4)){
    itemsm[item4, 4] <- 1
  }  
  vennDiagram(itemsm, names = names, main = main)
  return(itemsm)
}


#This function finds all the bins inside a region and compares them to ASpli events and the gold standar.
#In order to do this, the function first finds all the annotated features that fall withing a given region.
#Then it corrects exonic bins incorrectly labeled as intronic bins.
#To correct this, it findd all the external exons belonging to all the transcripts in the gene the intron belongs to.
#If the intron matches with one or more bins, it changes the intron for those bins.
#Finally, if a region has at least one bin belonging to the gold standar, then it keeps only those bins in 
#the region and discard the rest to avoid false positives. If there isnt at least a bin belonging to the gold
#standar, it returns all the bins in the region.
#The function plots a venn with ASpli, gold standard and method intersections and return the list of bins
#from the method.
comparisons <- function(regions, method, P, features, aspli_events, universe){
  
  #All the annotated features
  f <- as.data.frame(features@bins)
  
  #Finds annotated features within the regions
  regions_overlap <- as.data.frame(findOverlaps(GRanges(f$seqnames, IRanges(start=f$start, end=f$end)), regions, type = "within"))
  method_bins <- cbind(names(features@bins[regions_overlap$queryHits, ]), regions_overlap$subjectHits)
  
  #If several regions belong to the same cluster, we mark all the bins belonging to the same cluster.
  if("cluster" %in% colnames(as.data.frame(regions))){
    method_bins[, 2] <- regions$cluster[regions_overlap$subjectHits]
  }    
  method_bins <- unique(method_bins)
  
  #We discard bins not present in the universe
  method_bins <- method_bins[method_bins[, 1] %in% universe, ]
  method_bins <- correct_bin_labels(method_bins)
  method_bins <- filter_false_positives(method_bins, P)
  plot_comparison_venns(method_bins, method, aspli_events)
  return(method_bins)
}

#This function corrects incorrectly labeled bins and removes extra false possitives
correct_bin_labels <- function(method_bins){
  
  #We correct exonic bins incorrectly labeled as intronic bins
  #To correct this, we find all the external exons belonging to all the transcripts in the gene the intron belongs to.
  #If the intron matches with one or more bins, we change the intron for those bins.
  intronic_bins <- unique(as.character(method_bins[grep("I", method_bins[, 1]), 1]))
  for(b in intronic_bins){
    boi <- countsb(counts)[b, ]
    boi_regions <- GRanges(1, IRanges(boi$start, boi$end))
    bins_in_gen <- countsb(counts)[countsb(counts)[, "locus"] == strsplit2(b, ":")[, 1], ]
    external_bins_in_gen <- bins_in_gen[bins_in_gen$event == "external" & bins_in_gen$feature == "E", ]
    if(nrow(external_bins_in_gen) > 0){
      external_bins_in_gen_regions <- GRanges(1, IRanges(external_bins_in_gen$start, external_bins_in_gen$end))
      hits <- as.data.frame(findOverlaps(boi_regions, external_bins_in_gen_regions))
      if(nrow(hits) > 0){
        bi <- which(method_bins[, 1] == b)
        for(i in 1:nrow(hits)){
          method_bins <- rbind(method_bins, method_bins[bi, ])
          method_bins[nrow(method_bins), 1] <- rownames(external_bins_in_gen)[hits$subjectHits[i]]
        }
        method_bins <- method_bins[-bi, ]
      }
    }
  }
  return(method_bins)
}

#This function filters falses positives.
#If a region has at least one bin belonging to the gold standar, then we keep only those bins in the region and discard the rest
#to avoid overcounting false positives.
filter_false_positives <- function(method_bins, P){
  
  positive_method_bins <- data.frame(present=method_bins[, 1] %in% P, 
                                     method=method_bins[, 2])
  amount_of_positives <- aggregate(present ~ method, positive_method_bins, sum)
  
  amount_of_negatives <- aggregate(present ~ method, positive_method_bins, function(s){
    all(!s)
  })
  negative_clusters <- as.character(amount_of_negatives$method[amount_of_negatives$present == T])
  negative_method_bins <- unique(method_bins[method_bins[, 2] %in% negative_clusters, 1])
  
  amount_of_positives <- setNames(amount_of_positives$present, amount_of_positives$method)
  to_remove <- amount_of_positives[positive_method_bins$method[positive_method_bins$present == FALSE]]
  to_remove <- unique(names(to_remove)[to_remove > 0])
  positive_method_bins <- method_bins[!(positive_method_bins$present == FALSE & positive_method_bins$method %in% to_remove), 1]
  positive_method_bins <- setdiff(positive_method_bins, negative_method_bins)
  
  single_negative_per_group <- method_bins[method_bins[, 2] %in% negative_clusters, ]
  
  single_negative_per_group <- single_negative_per_group[!duplicated(single_negative_per_group[, 2]), 1]
  method_bins <- c(positive_method_bins, single_negative_per_group)
  return(method_bins)
}

#This function plots all the venns comparing ASpli, the gold standard and the method for both bins and genes
plot_comparison_venns <- function(bines_de_method, method, aspli_events){
  #layout(matrix(1:4, ncol=2, byrow=T))
  #Uso mismo universe que con aspli
  limma_venn(P, aspli_events, bines_de_method, 
             names = c("Sim", "ASpli", method), main = "")
  if(method == "LEAFCUTTER"){
    text(x = -3.5, y = 0, "Bines", cex = 2, srt = 90)    
  }
  
  limma_venn(unique(strsplit2(P, ":")[, 1]), unique(strsplit2(aspli_events, ":")[, 1]), unique(strsplit2(bines_de_method, ":")[, 1]), 
             names = c("Sim", "ASpli", method), main = "")
  if(method == "LEAFCUTTER"){
    text(x = -3.5, y = 0, "Genes", cex = 2, srt = 90)    
  }
  
  if(FALSE){
    limma_venn(bines_de_method, events_list$bin, union(events_list$anchor, events_list$locale),
               names = c(method, "ASpli - Bin", "ASpli - Anchor + Locale "), main = "Bines")
    
    limma_venn(unique(strsplit2(bines_de_method, ":")[, 1]), unique(strsplit2(events_list$bin, ":")[, 1]),
               unique(strsplit2(union(events_list$anchor, events_list$locale), ":")[, 1]),
               names = c(method, "ASpli - Bin", "ASpli - Anchor + Locale "), main = "Genes")
  }
}

#This function computes precision and recall from a gold standard and a method
precision_recall <- function(sim, x){
  item1 <- unique(sim)
  item2 <- unique(x)
  precision <- length(intersect(item1, item2))/length(item2)
  recall <- length(intersect(item1, item2))/length(item1)
  return(c(precision, recall, length(item2)))
}




## Load simulation data ####
load("SimuData/AT1_ES-IR-Alt3ss-Alt5ss.Rdata")

#Load gold standar with the simulated bins
(load("SimuData/lASbins.RData"))
gold_standard <- strsplit2(names(unlist(lASbins)), "[.]")[, 2]

# Initialize matrix with precision and recall for every method for bins and genes
pr_b <- matrix(0, ncol=4, nrow=6)
pr_g <- matrix(0, ncol=4, nrow=6)
colnames(pr_g) <- colnames(pr_b) <- c("method", "precision", "recall", "size")



## Load ASpli Results ####
if(FALSE){
  load("aaux/andresBNF3/features.RData")
  load("aaux/andresBNF3/counts.RData")
  load("aaux/andresBNF3/reportes.RData")
  j<-reportes[[1]]$j
  gb<-reportes[[1]]$gb
  sr<-reportes[[1]]$sr
  is<-reportes[[1]]$is
  
  load("ASpliRunOld/features.RData")
  load("ASpliRunOld/counts.RData")
  (load("ASpliRunOld/reportes.Rdata"))
  (load("ASpliRunOld/du.RData"))
}else{
  load("ASpliRun/features.RData")
  load("ASpliRun/counts.RData")
  load("ASpliRun/du.RData")
}

if(FALSE){

  old<-signals(is)
  new<-signals(isNew)
  
  signals(isNew0)[,table(otherSources,feature,useNA="ifany")]
  new[,table(otherSources,feature,useNA="ifany")]
  old[,table(otherSources,feature,useNA="ifany")]

  #Hay 284 eventos en is que no estan en isNew: son todos Alt 5'/3'
  #y me parece que estan en el viejo porque antes el programa tenia un bug
  #y no los detectaba correctamente
  new[is.na(feature) & otherSources==0, table(bin.event)]
  # bin.event
  # Alt 5'/3'       CSP        ES       IoR 
  #       24         4         6         8 
   
  old[is.na(feature) & otherSources==0, table(bin.event)]
  # bin.event
  # Alt 5'/3'       CSP        ES       IoR 
  #       308         4         6         8   
  
  #todos los que estan en el old, pero no en el old son jl=1 que no fueron combinados correctamente (ese era el bug corregido en 01/2020)
  old[!region%in%new$region]
  #                 region     locus bin.event b bjs ja jl otherSources
  # 1: 1:10087820-10088041 AT1G28710 Alt 5'/3' 0   0  0  1            0
  # 2: 1:10134536-10134686 AT1G29040 Alt 5'/3' 0   0  0  1            0
  # 3: 1:10177042-10177138 AT1G29120 Alt 5'/3' 0   0  0  1            0
  # 4: 1:10448745-10448849 AT1G29850 Alt 5'/3' 0   0  0  1            0
  # 5: 1:10449480-10449598 AT1G29850 Alt 5'/3' 0   0  0  1            0
  # ---                                                                 
  # 280:   1:9576805-9577051 AT1G27570 Alt 5'/3' 0   0  0  1            0
  # 281:   1:9663149-9663970 AT1G27752 Alt 5'/3' 0   0  0  1            0
  # 282:   1:9695641-9695892 AT1G27840 Alt 5'/3' 0   0  0  1            0
  # 283:     1:981160-981773 AT1G03860 Alt 5'/3' 0   0  0  1            0
  # 284:   1:9859460-9859893 AT1G28210 Alt 5'/3' 0   0  0  1            0  
  
  isMissed <-old[!region%in%new$region,]
    
  for(i in 1:nrow(isMissed)){
    llocus <- isMissed[i]$locus
    old[locus%in%llocus,1:12]
    new[locus%in%llocus,1:12]
  }
  
}

#We define the universe. 
universe <- rownames(countsb(counts))

#True positives are the ones in the gold standard present in the universe
P <- intersect(gold_standard, universe)

#We find all the regions associated with the bins in the gold standard.
gold_standard_regions  <- GRanges(1, IRanges(countsb(counts)[P, "start"], countsb(counts)[P, "end"]))


## Data integration ####

#We integrate gold standard region signals with ASpli signals in order to compare them
bcalcul <- FALSE
if(bcalcul){ 
  is      <- integrateSignals(sr, asd, otherSources = gold_standard_regions )
  # No intron-retention event was simulated for intronic bins (feature=I), so our gold-standard data set
  # does not incude them. However, overalpping external exons could induce coverage changes within I-labelled bins.
  # In order to account for this, we looked for all the external exons belonging to the transcripts in the gene the
  # intronic-bin belongs to. If the intron matched with one or more of these external  kind of bins, we exchanged the intron for those bins.
  intronic_bins <- sr@binbased$bin[grep("I", sr@binbased$bin)]
  for(b in intronic_bins){
    boi <- countsb(counts)[b, ]
    boi_regions <- GRanges(1, IRanges(boi$start, boi$end))
    bins_in_gen <- countsb(counts)[countsb(counts)[, "locus"] == strsplit2(b, ":")[, 1], ]
    external_bins_in_gen <- bins_in_gen[bins_in_gen$event == "external" & bins_in_gen$feature == "E", ]
    if(nrow(external_bins_in_gen) > 0){
      external_bins_in_gen_regions <- GRanges(1, IRanges(external_bins_in_gen$start, external_bins_in_gen$end))
      hits <- as.data.frame(findOverlaps(boi_regions, external_bins_in_gen_regions))
      if(nrow(hits) > 0){
        #return(rownames(bines_del_gen_externos)[hits$subjectHits[1]])
        bi <- which(sr@binbased$bin == b)
        for(i in 1:nrow(hits)){
          sr@binbased <- rbind(sr@binbased, sr@binbased[bi, ])
          sr@binbased$bin[nrow(sr@binbased)] <- rownames(external_bins_in_gen)[hits$subjectHits[i]]
        }
        sr@binbased <- sr@binbased[-bi, ]
      }
    }
  }
  
  intronic_bins <- is@signals$bin[grep("I", is@signals$bin)]
  intronic_bins <- intronic_bins[!is.na(intronic_bins)]
  for(b in intronic_bins){
    boi <- countsb(counts)[b, ]
    boi_regions <- GRanges(1, IRanges(boi$start, boi$end))
    bins_in_gen <- countsb(counts)[countsb(counts)[, "locus"] == strsplit2(b, ":")[, 1], ]
    external_bins_in_gen <- bins_in_gen[bins_in_gen$event == "external" & bins_in_gen$feature == "E", ]
    if(nrow(external_bins_in_gen) > 0){
      external_bins_in_gen_regions <- GRanges(1, IRanges(external_bins_in_gen$start, external_bins_in_gen$end))
      hits <- as.data.frame(findOverlaps(boi_regions, external_bins_in_gen_regions))
      if(nrow(hits) > 0){
        #return(rownames(bines_del_gen_externos)[hits$subjectHits[1]])
        bi <- which(is@signals$bin == b)
        for(i in 1:nrow(hits)){
          is@signals <- rbind(is@signals, is@signals[bi, ])
          is@signals$bin[nrow(is@signals)] <- rownames(external_bins_in_gen)[hits$subjectHits[i]]
        }
        is@signals <- is@signals[-bi, ]
      }
    }
  }
  save(j, gb, sr, is,file="ASpliRun/du2.RData")
}else{
  load("ASpliRun/du2.RData")  
}

if(FALSE){
  rn<-signals(isNew)$region
  ro<-signals(is)$region
  sum(!rn%in%ro)
  sum(!ro%in%rn)
  io <- which(!ro%in%rn)
  table(signals(is)[io,]$bin.event)
}

## ASpli Performance Figure ####
#We consider bins that overlap with locale J3 ranges as correctly detected by ASpli locale
#We find all the bins inside a given locale regions.
#In order to do this, the function first finds all the annotated features that fall withing a given region.
#If a region has at least one bin belonging to the gold standar, then it keeps only those bins in 
#the region and discard the rest to avoid false positives. If there isnt at least a bin belonging to the gold
#standar, it returns all the bins in the region.
if(FALSE){
  unannotated_regions <- is@signals$region[is@signals$jl == 1 ]
  unannotated_regions <- strsplit2(unannotated_regions, "[:]")
  
  unannotated_regions <- GRanges(unannotated_regions[, 1], IRanges(as.numeric(strsplit2(unannotated_regions[, 2], "[-]")[, 1]), 
                                                                   as.numeric(strsplit2(unannotated_regions[, 2], "[-]")[, 2])))
  unannotated_regions <- as.data.frame(findOverlaps(GRanges(features@bins@seqnames, features@bins@ranges), 
                                                    unannotated_regions, type = "within"))
}else{
  unannotated_regions <- is@signals$J3[is@signals$jl %in% c("1","*") ]
  #look for multiples J3, and parse them
  a<-strsplit2(unannotated_regions,";")
  if(ncol(a)>1){
    for(j in 2:ncol(a)){
      i2<-which(a[,j]!="")
      for(i in seq_along(i2)){
        aa<-rbind(strsplit2(a[i2[i],1],".",fixed=TRUE),
                  strsplit2(a[i2[i],j],".",fixed=TRUE))
        a[i2[i],1]<-paste("1",min(as.numeric(aa[,2])),max(as.numeric(aa[,3])),sep=".")
      }
    }
  }
  unannotated_regions <- a[,1]
  unannotated_regions <- strsplit2(unannotated_regions, ".",fixed=TRUE)
  unannotated_regions <- GRanges(unannotated_regions[, 1], IRanges(as.numeric(unannotated_regions[, 2]), 
                                                                   as.numeric(unannotated_regions[, 3])))
  unannotated_regions <- as.data.frame(findOverlaps(GRanges(features@bins@seqnames, features@bins@ranges), 
                                                    unannotated_regions, type = "within"))
}
locale_bins <- unique(cbind(names(features@bins[unannotated_regions$queryHits, ]), unannotated_regions$subjectHits))
locale_bins <- locale_bins[locale_bins[, 1] %in% universe, ]
locale_bins <- filter_false_positives(locale_bins, P)

#We compile the signals list 
events_list <- list(bin = is@signals$bin[is@signals$b != 0], 
                    anchor = is@signals$bin[is@signals$ja != 0],
                    locale = locale_bins,
                    otherSignals = P)

#upset(fromList(events_list), order.by = "freq")
svg("paperFigs/aspliPerformance.svg")
a<-limma_venn(events_list$otherSignals, events_list$locale, events_list$anchor, events_list$bin, 
              names = c("sim", "locale", "anchorage", "coverage"), main = "ASpli performance")
colnames(a)<-c("sim", "locale", "anchorage", "coverage")
dev.off()




## Bin coverage Detection Call Figures ####
# Fig 1.b
#We separate signals from coverage (fold change) and inclusion (change in junction inclusion).
#In both cases we use significative bins (bin.fdr < 0.05) and default values.
#For fold change we use 3 and for change in inclusion we use 0.2.
#svg("aspliSyntheticCoverageDetails.svg")
svg("paperFigs/binCoverageDetectionCall.svg")
Counts <- vennCounts(limma_venn(P, sr@binbased$bin[which(sr@binbased$bin.fdr < 0.05 & abs(sr@binbased$bin.logFC) >= log2(3))], 
                                sr@binbased$bin[which(sr@binbased$bin.fdr < 0.05 & abs(rowSums(cbind(sr@binbased$cluster.dPIR, sr@binbased$cluster.dPSI), na.rm = T)) >= 0.2)], 
                                names = c("sim", "logFC",expression(paste(Delta,"inclusion"))),main = ""))
dev.off()




## ROC and Precision Recall figures####
#In order to characterize coverage signal detection and its dependance of fold change, we build the ROC curve 
#and Precision-Recall curve.
#To make sence of this curve, we take as universe only significative bins (bin.fdr < 0.05) so we can focus on fold change.
#We swipe logFC from 1 to max(logFC) + 0.5.
universe_bin <- sr@binbased$bin[which(sr@binbased$bin.fdr < 0.05)]
N <- universe_bin[!universe_bin %in% P]
Pbin <- P[P %in% universe_bin]
#limma_venn(P, positives, names=c("P", "ASpli"), main ="universe nuevo")
#N <- setdiff(universe, P)

maximo_lfc <- max(abs(sr@binbased$bin.logFC), na.rm = T) + 0.5

puntos <- seq(1, 2**maximo_lfc, 0.01)
PPV <- FPR <- TPR <- c()
for(fc in puntos){
  positives <- sr@binbased$bin[which(abs(sr@binbased$bin.logFC) >= log2(fc) & sr@binbased$bin.fdr < 0.05)]
  negatives <- setdiff(universe_bin, positives)
  TP   <- sum(positives %in% Pbin)
  FN   <- sum(negatives %in% Pbin)
  TN   <- sum(negatives %in% N)
  FP   <- sum(positives %in% N)
  TPR  <- c(TPR, TP/(TP + FN))
  FPR  <- c(FPR, FP/(TN + FP))
  if(length(positives) > 0){
    PPV  <- c(PPV, TP/(TP + FP))
  }else{
    PPV  <- c(PPV, 0)
  }
}


#plot(FPR, TPR, type="l", main=paste("ROC señal Bin (logFC)", "AUC:", format(auc, digits=2)))
fprs <- c()
ppvs <- c()
tprs <- c()
l<-which(puntos == 1.5)
#points(FPR[l], TPR[l], col="red")
fprs <- c(fprs, FPR[l])
tprs <- c(tprs, TPR[l])
ppvs <- c(ppvs, PPV[l])
l<-which(puntos == 2)
#points(FPR[l], TPR[l], col="green")
fprs <- c(fprs, FPR[l])
tprs <- c(tprs, TPR[l])
ppvs <- c(ppvs, PPV[l])
l<-which(puntos == 3)
#points(FPR[l], TPR[l], col="blue")
fprs <- c(fprs, FPR[l])
tprs <- c(tprs, TPR[l])
ppvs <- c(ppvs, PPV[l])

#ROC
auc <- auc(FPR, TPR)
auc <- 0.91
impares <- seq(1, length(FPR), by=2)
pares <- impares + 1
delta_x <- sum((FPR[pares]-FPR[impares])*0.5*(TPR[pares]+TPR[impares]))

svg("paperFigs/aspliSyntheticROC.svg")
plot(FPR, TPR, type="l", main=paste("bin-coverage logFC signal"))
text(0.91,0.95,paste0("(AUC:", format(auc, digits=2),")"),cex=0.8)
points(fprs, tprs, type = "p", pch=c(15, 17, 18), cex=c(1.5,1.5,2.5))
legend("bottomright", legend = c("FC*=1.5x", "FC*=2x", "FC*=3x"), pch=c(15, 17, 18), cex=1, bty="n", inset = 0.02)
dev.off()

#PR
Recall = TPR
Precision = PPV
#auc <- auc(PPV, TPR)
svg("paperFigs/aspliSyntheticRecall.svg")
plot(Recall, Precision, type="l", main="bin-coverage logFC signal")
points(tprs, ppvs, type = "p", pch=c(15, 17, 18), cex=c(1.5,1.5,2.5))
legend("bottomright", legend = c("FC*=1.5x", "FC*=2x", "FC*=3x"), pch=c(15, 17, 18), cex=1, bty="n", inset = 0.02)
 # d <- data.frame(Recall = TPR, Precision = PPV)
 # p <- ggplot(data = d, aes(x=Recall, y=Precision)) + geom_line() + theme_light()
 # dp <- data.frame(ppvs = ppvs, tprs = tprs, logFC = c("1.5", "2", "3"))
 # p <- p + geom_point(data = dp, mapping = aes(x = tprs, y = ppvs, shape = logFC), size = 4) + theme_light()
 # p <- p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme(text = element_text(size = 15), axis.title=element_text(size=10))
 # p + theme(legend.position = c(0.2, 0.2))
dev.off()



## Table (1) ####
#Look for TP, (exlcusively detected TP) and FP
for(i in 2:4){
  tp  <- sum(a[,i]!=0 & a[,1]!=0)
  tpx <- sum(a[,i]!=0 & apply(a,1,sum)==2 & a[,1]!=0)
  fp  <- sum(a[,i]!=0 & a[,1]==0)
  cat(colnames(a)[i],tp,tpx,fp,"\n")
}

#How many simulated bins were filter-out by denisty/low-coverage filtering
analyzedBins <- unique(c(sr@binbased$bin,sr@localebased$bin,sr@anchorbased$bin))
cat("# genomic bins:",nrow(countsb(counts)),"\n")
cat("# analyzed bins:",length(analyzedBins),"\n")
cat("# simulated bins:",length(events_list$otherSignals),"\n")
outSim<-setdiff(events_list$otherSignals,analyzedBins)
cat("# fitered-out simulated bins:",length(outSim),"\n")

#Every simulated filteredout bin is included in the 1477 set, 
#so the filtering step is responsible for 1341/1477 FNs
length(setdiff(outSim,events_list$otherSignals)) == 0



## COMPARISON WITH OTHER METHODS #####
#We join all ASpli events and calculate ASpli precision and recall for bins and genes
aspli_events   <- unique(union(union(events_list$bin, events_list$anchor), events_list$locale))
averGenesASpli <- intersect(rownames(lmRNAcounts[[1]]$g)[lmRNAcounts[[1]]$g$as.active==0],strsplit2(aspli_events, ":")[,1])

pr_b[1, ] <- c("ASpli", precision_recall(P, aspli_events))
pr_g[1, ] <- c("ASpli", precision_recall(unique(strsplit2(P, ":")[, 1]), unique(strsplit2(aspli_events, ":")[, 1])))

aspli_events <- unique(events_list$bin)
averGenesASpli <- intersect(rownames(lmRNAcounts[[1]]$g)[lmRNAcounts[[1]]$g$as.active==0],strsplit2(aspli_events, ":")[,1])

pr_b[5, ] <- c("ASpli - Coverage", precision_recall(P, aspli_events))
pr_g[5, ] <- c("ASpli - Coverage", precision_recall(unique(strsplit2(P, ":")[, 1]), unique(strsplit2(aspli_events, ":")[, 1])))

aspli_events <- unique(union(events_list$anchor, events_list$locale))
averGenesASpli <- intersect(rownames(lmRNAcounts[[1]]$g)[lmRNAcounts[[1]]$g$as.active==0],strsplit2(aspli_events, ":")[,1])
pr_b[6, ] <- c("ASpli - Junction", precision_recall(P, aspli_events))
pr_g[6, ] <- c("ASpli - Junction", precision_recall(unique(strsplit2(P, ":")[, 1]), unique(strsplit2(aspli_events, ":")[, 1])))

layout(matrix(1:6, ncol=3, byrow=F))

## Load LEAFCUTTER results
{
  #Read Leafcutter output
lc_c <- read_delim("others/LEAFCUTTER/leafcutter_ds_cluster_significance.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)
lc_e <- read_delim("others/LEAFCUTTER/leafcutter_ds_effect_sizes.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)
PSI_MIN=0.1

#Keep significant clusters (p.adjust < 0.05) and find the regions associated with them
signif_clusters <- gsub("chr1:", "", lc_c$cluster[which(lc_c$p.adjust< 0.05)], fixed="T")
regions <- lc_e[sapply(strsplit(lc_e$intron, ":", fixed=T), "[[", 4) %in% signif_clusters,]
regions <- matrix(unlist(lapply(strsplit(regions$intron[abs(regions$deltapsi) > PSI_MIN], ":", fixed=T), function(s){return(s[2:4])})), ncol = 3, byrow=T)
regions <- GRanges(1, IRanges(start=as.numeric(regions[, 1]), end = as.numeric(regions[, 2])), cluster = regions[, 3])

#Make venn diagrams and get all the bins associated with the regions
leafcutter_bins <- comparisons(regions, "LEAFCUTTER", P, features, aspli_events, universe)
averGenesLC     <- intersect(rownames(lmRNAcounts[[1]]$g)[lmRNAcounts[[1]]$g$as.active==0],strsplit2(leafcutter_bins, ":")[,1])
#Calculate precision and recall
pr_b[2, ] <- c("LEAFCUTTER", precision_recall(P, leafcutter_bins))
pr_g[2, ] <- c("LEAFCUTTER", precision_recall(unique(strsplit2(P, ":")[, 1]), unique(strsplit2(leafcutter_bins, ":")[, 1])))
## LEAFCUTTER done
}

## Load MAJIQ results
{
  #Read majiq output
majiq <- read_delim("others/MAJIQ/04_dpsi/A_B.deltapsi.tsv", 
                    "\t", escape_double = FALSE, trim_ws = TRUE)

#We keep all the junctions it at least one LSV junction has a P(|dPSI|>=0.20) per LSV junction greater than 0.95
regions <- majiq$`Junctions coords`[unlist(lapply(strsplit(majiq$`P(|dPSI|>=0.20) per LSV junction`, ";", fixed = T), function(s){
  return(any(s > 0.95))
}))]

#Split all the regions in the junction and keep the largest range 
regions <- strsplit(regions, ";", fixed=T)
regions <- matrix(unlist(lapply(regions, function(r){
  r <- strsplit2(r, "-")
  return(c(min(r[, 1]), max(r[, 2])))
})), ncol=2, byrow=T)
regions <- GRanges(1, IRanges(start=as.numeric(regions[, 1]), end = as.numeric(regions[, 2])))

#Compares with gold standard and ASpli and returns all the bins found in the regions
majiq_bines <- comparisons(regions, "MAJIQ", P, features, aspli_events, universe)
averGenesMAJIQ <- intersect(rownames(lmRNAcounts[[1]]$g)[lmRNAcounts[[1]]$g$as.active==0],strsplit2(majiq_bines, ":")[,1])


#Calculates precision recall for majiq
pr_b[3, ] <- c("MAJIQ", precision_recall(P, majiq_bines))
pr_g[3, ] <- c("MAJIQ", precision_recall(unique(strsplit2(P, ":")[, 1]), unique(strsplit2(majiq_bines, ":")[, 1])))
## MAJIQ done
}

## Load rMATS results
{
  #We read every output in rmats. We keep those regions with an absolute inclusion difference level greater than 0.2
regions_all <- matrix(ncol=2, nrow=0)
PSI_MIN=0.1

f <- "A3SS"  
rm <- read_table2(paste0("others/rMATS/", f, ".MATS.JCEC.txt"))

regions <- rm[rm$FDR < 0.05 & abs(as.numeric(rm$IncLevelDifference)) > PSI_MIN & rm$strand == "-", c("shortEE", "longExonEnd")]
regions_all <- rbind(regions_all, as.matrix(regions, ncol=2, byrow = T))
regions <- rm[rm$FDR < 0.05 & abs(as.numeric(rm$IncLevelDifference)) > PSI_MIN & rm$strand == "+", c("longExonStart_0base", "shortES")]
regions_all <- rbind(regions_all, as.matrix(regions, ncol=2, byrow = T))


f <- "A5SS"  
rm <- read_table2(paste0("others/rMATS/", f, ".MATS.JCEC.txt"))


regions <- rm[rm$FDR < 0.05 & abs(as.numeric(rm$IncLevelDifference)) > PSI_MIN & rm$strand == "+", c("shortEE", "longExonEnd")]
regions_all <- rbind(regions_all, as.matrix(regions, ncol=2, byrow = T))
regions <- rm[rm$FDR < 0.05 & abs(as.numeric(rm$IncLevelDifference)) > PSI_MIN & rm$strand == "-", c("longExonStart_0base", "shortES")]
regions_all <- rbind(regions_all, as.matrix(regions, ncol=2, byrow = T))

f <- "MXE"  
rm <- read_table2(paste0("others/rMATS/", f, ".MATS.JCEC.txt"))

regions <- rm[rm$FDR < 0.05 & abs(as.numeric(rm$IncLevelDifference)) > PSI_MIN, c("1stExonStart_0base", "1stExonEnd")]
regions_all <- rbind(regions_all, as.matrix(regions, ncol=2, byrow = T))

regions <- rm[rm$FDR < 0.05 & abs(as.numeric(rm$IncLevelDifference)) > PSI_MIN, c("2ndExonStart_0base", "2ndExonEnd")]
regions_all <- rbind(regions_all, as.matrix(regions, ncol=2, byrow = T))

f <- "SE"  
rm <- read_table2(paste0("others/rMATS/", f, ".MATS.JCEC.txt"))

regions <- rm[rm$FDR < 0.05 & abs(as.numeric(rm$IncLevelDifference)) > PSI_MIN, c("exonStart_0base", "exonEnd")]
regions_all <- rbind(regions_all, as.matrix(regions, ncol=2, byrow = T))

regions_all[, 1] <- regions_all[, 1] + 1
regions_all <- GRanges(1, IRanges(start=as.numeric(regions_all[, 1]), end = as.numeric(regions_all[, 2])))

overlaps <- as.data.frame(findOverlaps(regions_all, GRanges(features@bins@seqnames, features@bins@ranges), type = "equal"))
overlaps_within <- as.data.frame(findOverlaps(GRanges(features@bins@seqnames, features@bins@ranges), regions_all[which(!1:length(regions_all) %in% overlaps$queryHits)], type = "within"))
colnames(overlaps_within) <- c("subjectHits", "queryHits")
overlaps <- rbind(overlaps, overlaps_within)

rmats_bins <- as.data.frame(features@bins[sort(unique(overlaps$subjectHits)),])
rmats_bins <- rownames(rmats_bins)
rmats_bins2 <- data.frame(rmats_bins = rmats_bins, clusters = 1:length(rmats_bins))
#Make comparisons between rMats, ASpli and gold standard
rmats_bins2 <- as.character(correct_bin_labels(rmats_bins2)[, 1])
averGenesRMATS <- intersect(rownames(lmRNAcounts[[1]]$g)[lmRNAcounts[[1]]$g$as.active==0],strsplit2(rmats_bins, ":")[,1])

plot_comparison_venns(rmats_bins2, "rMATS", aspli_events)
pr_b[4, ] <- c("rMATS", precision_recall(P, rmats_bins))
pr_g[4, ] <- c("rMATS", precision_recall(unique(strsplit2(P, ":")[, 1]), unique(strsplit2(rmats_bins, ":")[, 1])))
## rMATS done
}

#  Priduce Precision and Recall Table(2)   #
#Save precision and recall
pr_b <- as.data.frame(pr_b)
pr_g <- as.data.frame(pr_g)
if(FALSE)save(pr_b, pr_g, file="pr.RData")

prg <- pr_g
prb <- pr_b
sizeg <- as.numeric(as.character(prg$size))
sizeb <- as.numeric(as.character(prb$size))
names(sizeg)<-names(sizeb)<-prg$method
smax <- max(c(sizeg,sizeb))
sizeg <- sizeg/smax*3+.5
sizeb <- sizeb/smax*3+.5

x <- prb
x <- x[c(1,5,6,2,3,4),]
xx <- x
xx[,2]<-signif(as.numeric(as.character(xx[,2])),2)
xx[,3]<-signif(as.numeric(as.character(xx[,3])),2)

knitr::kable(xx[,c(1,4,2,3)],"latex")
# \begin{tabular}{l|l|l|r|r}
# \hline
# & method & size & precision & recall\\
# \hline
# 1 & ASpli & 1022 & 0.95 & 0.400\\
# \hline
# 5 & ASpli - Coverage & 966 & 0.96 & 0.380\\
# \hline
# 6 & ASpli - Junction & 583 & 0.99 & 0.230\\
# \hline
# 2 & LEAFCUTTER & 204 & 0.93 & 0.077\\
# \hline
# 3 & MAJIQ & 538 & 0.84 & 0.180\\
# \hline
# 4 & rMATS & 405 & 0.87 & 0.140\\
# \hline
# \end{tabular}  
# 

x <- prg
x <- x[c(1,5,6,2,3,4),]
xx <- x
xx[,2]<-signif(as.numeric(as.character(xx[,2])),2)
xx[,3]<-signif(as.numeric(as.character(xx[,3])),2)
knitr::kable(xx[,c(1,4,2,3)],"latex")

# \begin{tabular}{l|l|l|r|r}
# \hline
# & method & size & precision & recall\\
# \hline
# 1 & ASpli & 631 & 0.99 & 0.68\\
# \hline
# 5 & ASpli - Coverage & 591 & 0.99 & 0.64\\
# \hline
# 6 & ASpli - Junction & 456 & 0.99 & 0.50\\
# \hline
# 2 & LEAFCUTTER & 163 & 0.91 & 0.16\\
# \hline
# 3 & MAJIQ & 381 & 0.87 & 0.36\\
# \hline
# 4 & rMATS & 352 & 0.91 & 0.35\\
# \hline
# \end{tabular}

