library(ASpli)
library(SuperExactTest)
library(data.table)
library(GenomicRanges)

## Load results ####

## ASpli 
{
  universeASpli <- c() 

  ## Data Set A
  load("data/A/aspli/A_report.Rdata")

  sr.A            <- sr                 # splicing report object for experiment A 
  is.A            <- signals(is)        # ASpli signal table for experiment A
  as.A            <- unique(is.A$locus) # gene-level numbers
  universeASpli.A <- binsDU(gb)$locus   # background list statistically tested

  # estimation of the number of examined genomic regions: #bines + #locales - #annotatedES - #annotatedAlt
  nregsA <- nrow(binsDU(gb)) + nrow(localej(jdur)) - nrow(jes(jdur)) -nrow(jalt(jdur))


  ## Data Set B
  load("data/B/aspli/B_report.Rdata")

  sr.B            <- sr                  #splicing report for dataset B 
  is.B            <- signals(is)         #splicing signals for dataset B
  as.B            <-unique(is.B$locus)
  universeASpli.B <- binsDU(gb)$locus

  universeASpli<-unique(c(universeASpli.A,universeASpli.B)) #complete ASpli background list 

  # estimation of the number of examined genomic regions: #bines + #locales - #annotatedES - #annotatedAlt
  nregsB <- nrow(binsDU(gb)) + nrow(localej(jdur)) - nrow(jes(jdur)) -nrow(jalt(jdur))
  nregs  <- max(nregsA,nregsB)

  cat("ASpli results...loaded\n")
}

## Load other algorithm analysis
## LeafCutter
{
 saux    <- "data/geneName_geneID.txt"
 name2id <- fread(saux,header=TRUE,sep="\t")
 
 lc.deltapsi <- 0.1    # minimum deltapsi level to filter out splicing events


 ## Load Data Set A
 ppath <- "data/A/leafCutter/"
 lc1 <- read.table(paste0(ppath,"leafcutter_ds_cluster_significance.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="")
 lc2 <- read.table(paste0(ppath,"leafcutter_ds_effect_sizes.txt"),sep="\t",header=TRUE,quote="")
 
 universeLeafCutter.A     <- unique(lc1$genes)                  # how many genes were tested
 universeLeafCutter.reg.A <- strsplit2(lc2$intron,":clu")[,1]   #how many regions where tested


 # identify significative changes
 lc1 <- lc1[lc1$p.adjust < 0.05 & !is.na(lc1$p.adjust),]
 clusters1 <- strsplit2(lc1$cluster,":")[,2]
 clusters2 <- strsplit2(lc2$intron,":")[,4]
 i2 <- match(clusters2,clusters1)
 i2a <- which(!is.na(i2)) 
 i2b <- which(abs(lc2[,"deltapsi"])>lc.deltapsi)
 lc.A.reg <- strsplit2(lc2$intron[intersect(i2a,i2b)],":clu")[,1]  # splicing-altered regions 

 # convert gene_name to gene_id
 a<-name2id[gene_name%in%strsplit2(lc1$genes,",")[,1],]
 imultiple <- which(a[,table(gene_id)]>1)
 if(length(imultiple)>0){
  cat("hay gene_names con mas de un gene_id asociado...cual elijo?\n")
  a[gene_id%in%names(imultiple)] 
  name2id[name2id$gene_name%in%lc1$genes,]
 }
 lc.A<-a[,unique(gene_id)]            #genes hosting splicing changes

 
 
 ## Load Data set B
 ppath <- "data/B/leafCutter/"
 lc1 <- read.table(paste0(ppath,"leafcutter_ds_cluster_significance.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="")
 lc2 <- read.table(paste0(ppath,"leafcutter_ds_effect_sizes.txt"),sep="\t",header=TRUE,quote="")


 universeLeafCutter.B <- lc1$genes                              # how many genes were tested            
 universeLeafCutter   <- unique(c(universeLeafCutter.A,universeLeafCutter.B))

 universeLeafCutter.reg.B <- strsplit2(lc2$intron,":clu")[,1]   # how many regions were tested
 universeLeafCutter.reg <- unique(c(universeLeafCutter.reg.A,universeLeafCutter.reg.B))

 #identify significative changes

 lc1 <- lc1[lc1$p.adjust < 0.05 & !is.na(lc1$p.adjust),]

 # convert gene_name a gene_id
 a<-name2id[gene_name%in%strsplit2(lc1$genes,",")[,1],]
 imultiple <- which(a[,table(gene_id)]>1)
 if(length(imultiple)>0){
  cat("hay gene_names con mas de un gene_id asociado...cual elijo?\n")
  a[gene_id%in%names(imultiple)] 
  name2id[name2id$gene_name%in%lc1$genes,] 
 }
 lc.B<-a[,unique(gene_id)]

 lc1 <- lc1[lc1$p.adjust < 0.05 & !is.na(lc1$p.adjust),]
 clusters1 <- strsplit2(lc1$cluster,":")[,2]
 clusters2 <- strsplit2(lc2$intron,":")[,4]
 i2 <- match(clusters2,clusters1)
 i2a <- which(!is.na(i2)) 
 i2b <- which(abs(lc2[,"deltapsi"])>lc.deltapsi)
 lc.B.reg <- strsplit2(lc2$intron[intersect(i2a,i2b)],":clu")[,1]


 cat("Leaf Cutter results...loaded\n")
}

## MAJIQ
{
  posteriorP<-0.95    #threshold for a posteriori probability
  
  ## Data set A
  ppath <- "data/A/majiq/"
  m1  <- read.table(paste0(ppath,"col0_prmt5.deltapsi.tsv"),sep="\t",stringsAsFactors=FALSE,header=TRUE,comment.char="@")
  
  universeMAJIQ.A     <- unique(m1$Gene.ID)   #number of genes statistically tested
  universeMAJIQ.reg.A <- unique(m1$LSV.ID)    #number of regions statistically tested
  
  #apply significance filter
  icol <- grep(".0.20.",colnames(m1))
  a<-m1[,icol]#$P..dPSI...0.20..per.LSV.junction
  a<-strsplit2(a,";")
  res<-c()
  for(i in 1:nrow(a)){
    n <- as.numeric(a[i,])
    res<-c(res,any(n>posteriorP,na.rm=TRUE))
  }
  mGenes<-m1[res,1]
  mGenes<-unique(mGenes)
  
  m.A     <- mGenes     #altered genes 
  m.A.reg <- m1[res,2]  #altered regions
  
  
  ## Data set B
  ppath <- "data/B/majiq/"
  m1  <- read.table(paste0(ppath,"col0_prmt5.deltapsi.tsv"),sep="\t",stringsAsFactors=FALSE,header=TRUE,comment.char="@")
  
  universeMAJIQ.B     <-unique(m1$Gene.ID)
  universeMAJIQ.reg.B <- unique(m1$LSV.ID)
  
  universeMAJIQ     <-unique(c(universeMAJIQ.A,universeMAJIQ.B))
  universeMAJIQ.reg <- unique(c(universeMAJIQ.reg.A,universeMAJIQ.reg.B))
  
  # significant changes
  icol <- grep(".0.20.",colnames(m1))
  a<-m1[,icol]#$P..dPSI...0.20..per.LSV.junction
  a<-strsplit2(a,";")
  res<-c()
  for(i in 1:nrow(a)){
    n <- as.numeric(a[i,])
    res<-c(res,any(n>posteriorP,na.rm=TRUE))
  }
  mGenes<-m1[res,1]
  mGenes<-unique(mGenes)
  
  m.B <- mGenes           #altered genes
  m.B.reg <- m1[res,2]    #altered regions
  
  cat("MAJIQ results...loaded\n")
}

## rMATS
{
# set parameters to identify significative splicing alterations
dPSI <- 0.1
qv   <- 0.05

#internal function that retunrs genomic regions satisfying  qv<qv_value y dpsi>dpsi_value
regionesRMATS<-function(rm,class=c("A3SS", "A5SS", "MXE", "RI", "SE")[4],qv=0.05,dpsi=0.2){
  if(class=="A3SS"){  
    regiones <- rm[rm$FDR < qv & abs(as.numeric(rm$IncLevelDifference)) > dpsi & rm$strand == "-", c("chr","shortEE", "longExonEnd")]
    maux <- rm[rm$FDR < qv & abs(as.numeric(rm$IncLevelDifference)) > dpsi & rm$strand == "+", c("chr","longExonStart_0base", "shortES")]
    colnames(maux)<-colnames(regiones)
    regiones <- rbind(regiones,maux)
    regiones <- apply(regiones,1,paste,collapse=".")
  }
  
  if(class=="A5SS"){  
    regiones <- rm[rm$FDR < qv & abs(as.numeric(rm$IncLevelDifference)) > dpsi & rm$strand == "+", c("chr","shortEE", "longExonEnd")]
    maux <- rm[rm$FDR < qv & abs(as.numeric(rm$IncLevelDifference)) > dpsi & rm$strand == "-", c("chr","longExonStart_0base", "shortES")]
    colnames(maux)<-colnames(regiones)
    regiones <- rbind(regiones,maux)
    regiones <- apply(regiones,1,paste,collapse=".")
  }
  
  if(class=="MXE"){  
    regiones <- rm[rm$FDR < qv & abs(as.numeric(rm$IncLevelDifference)) > dpsi, c("chr","X1stExonStart_0base", "X1stExonEnd")]
    maux <- rm[rm$FDR < qv & abs(as.numeric(rm$IncLevelDifference)) > dpsi, c("chr","X2ndExonStart_0base","X2ndExonEnd")]
    colnames(maux)<-colnames(regiones)
    regiones <- rbind(regiones,maux)
    regiones <- apply(regiones,1,paste,collapse=".")
  }
  
  if(class=="SE"){  
    regiones <- rm[rm$FDR < qv & abs(as.numeric(rm$IncLevelDifference)) > dpsi, c("chr","exonStart_0base", "exonEnd")]
    regiones <- apply(regiones,1,paste,collapse=".")
  }
  
  if(class=="RI"){  
    regiones <- rm[rm$FDR < qv & abs(as.numeric(rm$IncLevelDifference)) > dpsi, c("chr","riExonStart_0base", "riExonEnd")]
    regiones <- apply(regiones,1,paste,collapse=".")
  }
  
  return(gsub(" ","",regiones))
}


## Data set A
ppath <- "data/A/rMATS/"

#process and unify rMATS discoveries
asClass <- c("SE","RI","MXE","A5SS","A3SS")
rmad<-rmadEC<-rmadEC.reg<-universerMATS.reg<-universerMATS.reg1<-universerMATS.reg2<-c()
for(i in seq_along(asClass)){
  a  <- read.table(paste0(ppath,asClass[i],".MATS.JCEC.txt"),sep="\t",header=TRUE,stringsAsFactors = FALSE,quote="")
  a[,2]<-sub("\"","",sub("\"","",a[,2]))
  a[,3]<-sub("\"","",sub("\"","",a[,3]))  
    
  rmadEC     <- rbind(rmadEC,cbind(event=asClass[i],a[,c("ID","GeneID","FDR","IncLevelDifference")]))
  rmadEC.reg <- c(rmadEC.reg,regionesRMATS(a,asClass[i],qv=qv,dpsi=dPSI))
  
  universerMATS.reg1<-c(universerMATS.reg1,regionesRMATS(a,asClass[i],qv=2,dpsi =0))
}

universerMATS.A     <- unique(rmadEC$GeneID)      #tested genes in A
universerMATS.reg.A <- unique(universerMATS.reg1) #tested genomic regions in A

rmadEC<-rmadEC[!is.na(rmadEC$IncLevelDifference),]
rmats.A<-unique(rmadEC[rmadEC$FDR<qv& abs(rmadEC$IncLevelDifference)>dPSI,"GeneID"])  # genes hosting significative changes in A
rmats.A.reg <- unique(rmadEC.reg[rmadEC.reg!="NA.NA.NA"])                             # significative genomic regions in A


## Data set B
ppath <- "data/B/rMATS/"

asClass <- c("SE","RI","MXE","A5SS","A3SS")
for(i in seq_along(asClass)){
  a  <- read.table(paste0(ppath,asClass[i],".MATS.JCEC.txt"),sep="\t",header=TRUE,stringsAsFactors = FALSE,quote="")
  a[,2]<-sub("\"","",sub("\"","",a[,2]))
  a[,3]<-sub("\"","",sub("\"","",a[,3]))  
    
  rmadEC <- rbind(rmadEC,cbind(event=asClass[i],a[,c("ID","GeneID","FDR","IncLevelDifference")]))
  rmadEC.reg <- c(rmadEC.reg,regionesRMATS(a,asClass[i],qv=qv,dpsi=dPSI))
  
  universerMATS.reg2<-c(universerMATS.reg2,regionesRMATS(a,asClass[i],qv=2,dpsi =0))
}

universerMATS.B     <-unique(rmadEC$GeneID)       #tested genes in B
universerMATS.reg.B <- unique(universerMATS.reg2) #tested genomic regions in B

universerMATS     <- unique(c(universerMATS.A,universerMATS.B))         #overall tested genes
universerMATS.reg <- unique(c(universerMATS.reg.A,universerMATS.reg.B)) #overall tested regions

rmadEC<-rmadEC[!is.na(rmadEC$IncLevelDifference),]
rmats.B<- unique(rmadEC[rmadEC$FDR<qv & rmadEC$IncLevelDifference>dPSI,"GeneID"]) # genes hosting significative changes in B
rmats.B.reg <- unique(rmadEC.reg[rmadEC.reg!="NA.NA.NA"])                         # significative regions in B  

rmats.B <-rmats.B[!is.na(rmats.B)]
rmats.A <-rmats.A[!is.na(rmats.A)]


cat("rMATs results...loaded\n")
}

## PCR data
{
  pcr <- read.table("data/pcr/TablitaRTPCR_PRMT5.csv",sep="\t",header=TRUE,stringsAsFactors = FALSE,quote="")
  colnames(pcr)<-c("locus","event","isoformX","isoformsOthers","frac.wt.1","frac.wt.2","frac.wt.3",
                   "frac.prmt5.1","frac.prmt5.2","frac.prmt5.3","delta.pct","pv","description")
  pcr[,1] <- toupper(pcr[,1])
  pcr     <- pcr[,-ncol(pcr)]
  
  #  Xiang paper PCR validated splicing event
  xiang   <- c("AT5G05010","AT1G75670","AT1G68190","AT1G78630","AT1G69620","AT5G14740","AT5G55125","AT4G26570","AT2G17340","AT3G48360","AT3G52720","AT5G44640")
  
  cat("PCR results...loaded.\n")
}


## Reproducibility analysis ####

# Lists of discoveries at gene and genomic-region levels
lRNAseq     <- list(as.A=as.A,as.B=as.B,lc.A=lc.A,lc.B=lc.B,m.A=m.A,m.B=m.B,rmats.A=rmats.A,rmats.B=rmats.B)
lRNAseq.reg <- list(as.A=as.character(is.A$region),as.B=as.character(is.B$region),
                    lc.A=lc.A.reg,lc.B=lc.B.reg,
                    m.A=m.A.reg,m.B=m.B.reg,
                    rmats.A=rmats.A.reg,rmats.B=rmats.B.reg)


# Agreement analysis of discoveries between experiments (per methodology, GENE-level)  
{
  
  lUniverse<-list(as=universeASpli,lc=universeLeafCutter,m=universeMAJIQ,rm=universerMATS)
  ul <- unlist(lapply(lUniverse,length))
  lres<-list()
  
  # Agreement at gene-level
  mres.g<-c()
  for(i in seq_along(ul)){
    lres[[i]] <- supertest(lRNAseq[2*(i-1)+(1:2)],n=ul[i])
    a <- summary(lres[[i]])
    mres.g <- rbind(mres.g,c(universe=a$n,a$otab,eo=a$etab[3],pv=a$P.value[3],FE=a$otab[3]/a$etab[3]))
  }
  rownames(mres.g)<-names(lres)<-names(ul)
  colnames(mres.g)<-c("universe","B","A","A&B","EO","pv","FE")
  print(mres.g)
  
  #Union & different normalized intersection statistics (we used thhe last one in the paper)
  mres2.g<-t(apply(mres.g,1,function(x){
    aUb<-unname(sum(x[2:3])-x[4])
    a<-c(aub=aUb,a=x[4]/x[2],b=x[4]/x[3],agreement=x[4]/aUb,overlap_coef=x[4]/min(x[2:3]))
    names(a)<-c("AUB","AnB_A","AnB_B","AnB_AUB","AnB_minAB")
    return(a)
  }))  
  print(mres2.g)
  
  #    universe    B    A  A&B         EO            pv        FE
  # as    13764 2957 3544 2115 761.378088  0.000000e+00  2.777858
  # lc     3342  490  450  261  65.978456 2.403989e-126  3.955837
  # m      8000  264  258  150   8.514000 1.552294e-171 17.618041
  # rm     1797   63  108   52   3.786311  9.850303e-59 13.733686
  #     AUB     AnB_A     AnB_B   AnB_AUB AnB_minAB
  # as 4386 0.7152519 0.5967833 0.4822161 0.7152519
  # lc  679 0.5326531 0.5800000 0.3843888 0.5800000
  # m   372 0.5681818 0.5813953 0.4032258 0.5813953
  # rm  119 0.8253968 0.4814815 0.4369748 0.8253968
}

#Agreement analysis of discoveries between experiments (per methodology, GENOMI-REGION-level) paper's Table(3)  
{
  #number of examined regions  
  ul <- c(as=nregs,lc=length(universeLeafCutter.reg),m=length(universeMAJIQ.reg),rm=length(universerMATS.reg))
  
  #region-level agreement (Table-3 of the paper)
  lres<-list()
  mres<-c()
  for(i in seq_along(ul)){
    lres[[i]] <- supertest(lRNAseq.reg[2*(i-1)+(1:2)],n=ul[i])
    a <- summary(lres[[i]])
    mres <- rbind(mres,c(universe=a$n,a$otab,eo=a$etab[3],pv=a$P.value[3],FE=a$otab[3]/a$etab[3]))
  }
  rownames(mres)<-names(lres)<-names(ul)
  colnames(mres)<-c("universe","B","A","A&B","EO","pv","FE")
  print(mres)
  
  #estimation of overlap coefficient value (last column of mres2)
  mres2<-t(apply(mres,1,function(x){
    aUb<-unname(sum(x[2:3])-x[4])
    a<-c(aub=aUb,a=x[4]/x[2],b=x[4]/x[3],agreement=x[4]/aUb,overlap_coef=x[4]/min(x[2:3]))
    names(a)<-c("AUB","AnB_A","AnB_B","AnB_AUB","AnB_minAB")
    return(a)
  }))  
  print(mres2)
  #    universe    B    A  A&B         EO            pv        FE
  # as   140191 3904 4687 2241 130.522273  0.000000e+00 17.169483
  # lc     8113  675  603  327  50.169481 3.649021e-219  6.517907
  # m     16441  284  277  149   4.784867 9.488903e-203 31.139841
  # rm     2405  158  119  119   7.817879 7.849869e-168 15.221519
  #     AUB     AnB_A     AnB_B   AnB_AUB AnB_minAB
  # as 6350 0.5740266 0.4781310 0.3529134 0.5740266
  # lc  951 0.4844444 0.5422886 0.3438486 0.5422886
  # m   412 0.5246479 0.5379061 0.3616505 0.5379061
  # rm  158 0.7531646 1.0000000 0.7531646 1.0000000
}

#Graphical representation of number of discoveries and between-experiment agreement
{
  layout(matrix(1:2,1,2))
  ndiscoveries<-unlist(lapply(lRNAseq,length))
  b<-barplot(ndiscoveries,log="y",ylim=c(10,5000),main="genes")
  text(b,ndiscoveries*1.20,paste0("(",as.character(ndiscoveries),")"),cex=0.8)
  for(i in 1:4){
    ii<-2*(i-1)+(1:2)
    if(FALSE){
      rect(b[ii[1]],mres.g[i,4]*1.3,b[ii[2]],mres.g[i,4]*0.8,border=NA,col="lightgrey")
      lines(b[ii],rep(mres.g[i,4],2))
      text(mean(b[ii]),mres.g[i,4],mres.g[i,4],cex=0.8)
    }else{
      lines(b[ii],rep(mres.g[i,4],2),lwd=3)
    }
  }
  
  
  ndiscoveries<-unlist(lapply(lRNAseq.reg,length))
  b<-barplot(ndiscoveries,log="y",ylim=c(10,10000),main="regions")
  text(b,ndiscoveries*1.20,paste0("(",as.character(ndiscoveries),")"),cex=0.8)
  for(i in 1:4){
    ii<-2*(i-1)+(1:2)
    if(FALSE){
      rect(b[ii[1]],mres[i,4]*1.3,b[ii[2]],mres[i,4]*0.8,border=NA,col="lightgrey")
      lines(b[ii],rep(mres[i,4],2))
      text(mean(b[ii]),mres[i,4],mres[i,4],cex=0.8)
    }else{
      lines(b[ii],rep(mres[i,4],2),lwd=3)
    }
  }
}

#Graphical representation of expected and reported overlap
{
  ccex <- (mres2[,"AnB_minAB"]-min(mres2[,"AnB_minAB"]))/(max(mres2[,"AnB_minAB"])-min(mres2[,"AnB_minAB"])) 
  ccex <- 1 + 1.5*ccex
  layout(1)
  plot(mres[,5],mres[,4],pch=20,xlab="expected overlap",ylab="observed overlap",log="xy",cex=ccex,ylim=c(10,3000))
  text(mres[,5],mres[,4],c("ASpli","   LeafCutter","MAJIQ","rMATS"),pos = c(1,1,1,1),cex=0.8)
  text(mres[,5],mres[,4],paste0("(",signif(mres[,"FE"],2),"x)"),pos = 3,cex=0.6)
  
  abline(0,1,col="gray",lty=2)
}

#Between methods dicovery comparisons
{
  
  
  #Standarize genomic-region IDs: 
  laux <- lRNAseq.reg
  laux[1:2]<-lapply(laux[1:2],function(x){paste0("chr",x)})
  laux[3:4]<-lapply(laux[3:4],function(x){
    res<-c()
    for(i in seq_along(x)){
      ipos<-gregexpr(":",x[i])[[1]][2]
      res <- c(res,paste0(substr(x[i],1,(ipos-1)),"-",substr(x[i],ipos+1,nchar(x[i]))))
    } 
    return(res)
  })
  laux[5:6]<-lapply(laux[5:6],function(x){
    res<-c()
    for(i in seq_along(x)){
      ipos<-gregexpr(":",x[i])[[1]][2]
      chr<-substr(x[i],3,3)
      res<-c(res,paste0("chr",chr,":",substr(x[i],ipos+1,nchar(x[i]))))
    }
    return(res)
  })
  laux[7:8]<-lapply(laux[7:8],function(x){
    res<-c()
    for(i in seq_along(x)){
      
      ipos<-gregexpr(".",x[i],fixed=TRUE)[[1]]
      res<-c(res,paste0(substr(x[i],1,(ipos[1]-1)),":",substr(x[i],ipos[1]+1,ipos[2]-1),
                        "-",substr(x[i],ipos[2]+1,nchar(x[i]))))
      
    }                    
    return(res)
  })
  
  
  #Assign genomic ranges to region IDs
  lranges <-  lapply(laux,function(x){
    res<-c()
    for(i in seq_along(x)){
      ipos<-regexpr(":",x[i],fixed=TRUE)
      jpos<-regexpr("-",x[i],fixed=TRUE)
      res<-rbind(res,data.frame(seqname=substr(x[i],1,ipos-1),
                                start=substr(x[i],ipos+1,jpos-1),
                                end=substr(x[i],jpos+1,nchar(x[i])),stringsAsFactors = FALSE))
    }
    res$start <- as.numeric(res$start)
    res$end   <- as.numeric(res$end)
    ineg <- which(res$start>res$end)
    if(length(ineg)>0) res[ineg,2:3] <- res[ineg,3:2]
    
    
    a<-GRanges(res$seqname,IRanges(res$start,res$end))
    return(a)
  })
  
  # Loop to estimate different kind of region overlaps
  #  any    : any kind of partial overlap
  #  within : region-i within region-j
  #  withinS: symmetrized within overlap (used in the paper)
  mOver.within <- mOver.withinS <- mOver.any <- matrix(0,nrow=length(lranges),ncol=length(lranges))
  for(i in seq_along(lranges)){
    for(j in (i+1):length(lranges)){
      if(j>length(lranges)) next
      overlapij      <- data.frame(findOverlaps(lranges[[i]], lranges[[j]],type="within"))
      overlapji      <- data.frame(findOverlaps(lranges[[j]], lranges[[i]],type="within"))
      overlapeq      <- data.frame(findOverlaps(lranges[[j]], lranges[[i]],type="equal"))
      overlapany     <- data.frame(findOverlaps(lranges[[j]], lranges[[i]],type="any"))
      
      mOver.within[i,j]<-nrow(overlapij)
      mOver.within[j,i]<-nrow(overlapji)
      mOver.withinS[i,j]<-mOver.withinS[j,i] <- nrow(overlapij)+nrow(overlapji)-nrow(overlapeq)
      mOver.any[i,j]<-mOver.any[j,i] <- nrow(overlapany)
    }
  }  
  colnames(mOver.within) <- colnames(mOver.withinS) <- colnames(mOver.any) <- rownames(mOver.within) <- rownames(mOver.withinS) <- rownames(mOver.any) <- names(lranges)
  diag(mOver.within) <- diag(mOver.withinS) <- diag(mOver.any) <- unlist(lapply(lranges,length)) 
}

#Significance estimation and labeledHeatmap production (SupFig S10)
{

library(WGCNA)

#Function to estimate different type of agrrement coefficients:
# x: list of sets
# type:
#   intersection: |x_i & x_j|
#   union       : |x_i & x_j|/|x_i U x_j|
#   max         : |x_i & x_j|/max(x_i,x_j)
#   min         : |x_i & x_j|/min(x_i,x_j)
mysim <- function (x,type=c("intersect","union","max","min")[1]) {
  if (!is.list(x)) 
    stop("Input x must be list\n")
  nL = length(x)
  if (nL < 2) 
    stop("Input x should have at least two entries\n")
  x = lapply(x, unique)
  Mat = matrix(NA, nL, nL)
  colnames(Mat) = rownames(Mat) = names(x)
  diag(Mat) = 1
  if(type=="max"){
   for (i in 1:(nL - 1)) {
    for (j in (i + 1):nL) {
      Mat[i, j] = Mat[j, i] = sum(x[[i]] %in% x[[j]])/max(length(unique(x[[i]])), length(unique(x[[j]])))
    }
   }
  }else{
    if(type=="min"){
      for (i in 1:(nL - 1)) {
        for (j in (i + 1):nL) {
          Mat[i, j] = Mat[j, i] = sum(x[[i]] %in% x[[j]])/min(length(unique(x[[i]])), length(unique(x[[j]])))
        }
      }
    }else{
      if(type=="union"){
      for (i in 1:(nL - 1)) {
        for (j in (i + 1):nL) {
          Mat[i, j] = Mat[j, i] = sum(x[[i]] %in% x[[j]])/length(union(x[[i]], x[[j]]))
        }
      }
      }else{
      for (i in 1:(nL)) {
        for (j in i:nL) {
          Mat[i, j] = Mat[j, i] = sum(x[[i]] %in% x[[j]])
        }
      }
        
      }
    }   
  }  
  Mat
}

#Gene-level overlap matrix
tm1 <- mysim(lRNAseq,"intersect")
tm2 <- mysim(lRNAseq,"min")
tm<-apply(cbind(tm1,tm2),1,function(x){paste0(x[1:ncol(tm1)],"\n(",signif(x[ncol(tm2)+(1:ncol(tm1))],2),")")})
labeledHeatmap(log(tm1+0.1),xLabels=colnames(tm1),textMatrix = tm,plotLegend = FALSE,main="gene agreement",cex.text=1.2)

#Lets start the genomic-region level analysis
universeNumbers.reg <- c(as.A=nregsA,as.B=nregsB,
                         lc.A=length(universeLeafCutter.reg.A),lc.B=length(universeLeafCutter.reg.B),
                         m.A =length(universeMAJIQ.reg.A),m.B =length(universeMAJIQ.reg.B),
                         rmats.A=length(universerMATS.reg.A),rmats.B=length(universerMATS.reg.B))

# Funcion para estimar significancia de overlaps entre descubrimientos realizados por diferentes metodos 
# En la hip. nula se asume que los unieversos difieren solo por cuestion de tamanio. Esto da una estimacion optimista de
# H0 y deberia servir como cota minim a de nivel de significancia (en la cuenta, si los universos de partida
# difieren no solo en numero sino tambien en composicion
# el overlap esperado por azar seria mas bajo aun)
#NOTA: deberia hacer la cuenta bien: agarrar regiones al azar de cada metodo y comparar overlaps
mytest<-function(hits1,hits2,hits12,u1,u2,N=1000){
  res <- rep(0,N)
  for(i in 1:N){
   a <- sample(1:u1,hits1)  
   b <- sample(1:u2,hits2)
   res[i]<-length(intersect(a,b))
  }
  return(c(mu0=mean(res),sd0=sd(res),z=(hits12-mean(res))/sd(res)))
}


N<-1000
for(k in 1:2){
  if(k==1){
    m <- mOver.withinS
  }else{  
    m <- mOver.any
  }
  for(i in 1:8){
    for(j in (i+1):8){
      cat(colnames(m)[i],colnames(m)[j]," ")
      if(i==8 & j>7) next
      u1 <- universeNumbers.reg[[i]]
      u2 <- universeNumbers.reg[[j]]
      hits12 <- m[i,j]
      hits1  <- m[i,i]
      hits2  <- m[j,j]
      a<-mytest(hits1,hits2,hits12,u1,u2,N=N)
      if(i==1 & j==2){
        mres <- data.frame(method1=colnames(m)[i],method2=colnames(m)[j],mu0=a[1],sd0=a[2],hits12=hits12,z=a[3],pv=1-pnorm(hits12,a[1],a[2]),stringsAsFactors = FALSE)
      }else{
        mres <- rbind(mres,data.frame(method1=colnames(m)[i],method2=colnames(m)[j],mu0=a[1],sd0=a[2],hits12=hits12,z=a[3],pv=1-pnorm(hits12,a[1],a[2]),stringsAsFactors = FALSE))
      }   
      cat(a,"\n")
    }   
  }
  mres <- cbind(mres,qv=p.adjust(mres$pv,"bonferroni"))
  if(k==1){
    mres1k.ws  <- mres
  }else{
    mres1k.any <- mres
  }
}


##overlap-any
mm<-mOver.any
#Normalize each cell-ij by min(cell-ii, cell-jj) in order to estimate the overlap coefficient
for(i in 1:nrow(mm)){
  for(j in (i+1):ncol(mm)){
    if(j>ncol(mm)) next
    mm[i,j]<-mm[j,i]<-mm[i,j]/min(mm[i,i],mm[j,j])
  }
}
diag(mm)<-1
msignif <- mres1k.any
tm1 <- mOver.any
tm2 <- mm
tm<-apply(cbind(tm1,tm2),1,function(x){paste0(x[1:ncol(tm1)],"\n(",signif(x[ncol(tm2)+(1:ncol(tm1))],2),")")})
rownames(tm)<-colnames(tm)
for(i in 1:nrow(msignif)){
  m1 <- as.character(msignif$method1[i])
  m2 <- as.character(msignif$method2[i])
  saux <- tm[m1,m2]
  if(msignif$qv[i] < 0.01) tm[m1,m2] <- tm[m2,m1] <- sub("\n","*\n",saux) 
}
labeledHeatmap(log(mOver.any+0.1),xLabels=colnames(mOver.any),textMatrix = tm,plotLegend = FALSE,main="region overlap (any)",cex.text=1)

##overlap-within
mm<-mOver.withinS
#Normalize each cell-ij by min(cell-ii, cell-jj) in order to estimate the overlap coefficient
for(i in 1:nrow(mm)){
  for(j in (i+1):ncol(mm)){
    if(j>ncol(mm)) next
    mm[i,j]<-mm[i,j]/min(mm[i,i],mm[j,j])
    mm[j,i]<-mm[j,i]/min(mm[i,i],mm[j,j])
  }
}
diag(mm)<-1
msignif <- mres1k.ws
tm1     <- mOver.withinS
tm2     <- mm
tm      <-apply(cbind(tm1,tm2),1,function(x){paste0(x[1:ncol(tm1)],"\n(",signif(x[ncol(tm2)+(1:ncol(tm1))],2),")")})
rownames(tm)<-colnames(tm)
for(i in 1:nrow(msignif)){
  m1 <- as.character(msignif$method1[i])
  m2 <- as.character(msignif$method2[i])
  saux <- tm[m1,m2]
  if(msignif$qv[i] < 0.05) tm[m1,m2] <- tm[m2,m1] <- sub("\n","*\n",saux) 
}
labeledHeatmap(log(tm1+0.1),xLabels=colnames(tm1),textMatrix = tm,plotLegend = FALSE,main="region overlap (within)",cex.text=1.2)
}
 

# Reproducibility and ASpli data integration: 
# ASpli was used to consolidate datasets A and B considering 'experiment' as a fixed effect:
# $ \sim experiment + genotype$. 
{
  #load ASpli consolidated analysis
  saux <- "data/AB/reportes.Rdata"
  (load(saux))

  sr.AB <- reportes$prmt5_expEffect$sr                # splicing report object for the consolidated dataset
  is.AB <- signals(reportes[["prmt5_expEffect"]]$is)  #signal table
  as.AB <- unique(is.AB$locus)                        #genes harboring splicing alterations according the consolidated dataset 
}

# There were only 18 regions that display interaction effects
signals(reportes$prmt5_interaction$is)


# We want to enforce the exclusion of non-coherent signals, so we filtered out regions exhibiting
# mild interaction effect signals (fdr < 0.5) (Supplementary Figure S7)
{

 library(GenomicAlignments)
  
 ## Lets filter-out regions with bin signals with fdr>fdrb for the interaction term  (48 out of 4360 for fdrb=0.5)
 fdrb  <- 0.5
 a     <- binbased(reportes$prmt5_interaction$sr)
 bfdrb <- a[a$bin.fdr<fdrb,"bin"]
 bbinout <- bfdrb[!is.na(bfdrb)]
 #cuantos bines detectados en AB tienen efecto interaccion con fdr>0.5
 is.AB[,table(bin%in%bbinout)]
 binABok <- is.AB[!bin%in%bbinout,as.character(region)]
 
 
 ## Analyzing locale signals to exclude interaction events (28 out of 187)
 a      <- localebased(reportes$prmt5_interaction$sr)
 a      <- a[a$cluster.fdr < fdrb,]
 a      <- a[!duplicated(a$cluster.range),]
 locInt <- strsplit2(as.character(unique(a$cluster.range)),".",fixed=TRUE)
 regiones_int <- GRanges(locInt[,1], IRanges(start=as.numeric(locInt[,2]), end = as.numeric(locInt[,3])))

 aa <- is.AB[is.AB$jl==1,]
 locAvg <- strsplit2(as.character(unique(aa$J3)),".",fixed=TRUE)
 regiones_avg <- GRanges(locAvg[,1], IRanges(start=as.numeric(locAvg[,2]), end = as.numeric(locAvg[,3])))
 #regiones con posible senial de interaccion (cluster.fdr < fdrtb)
 regiones_overlap <- as.data.frame(findOverlaps(regiones_int,regiones_avg,type="within"))
 locABok <- aa[-regiones_overlap[,2],as.character(region)]

  
 ## Analyzing anchorage signal
 a  <- anchorbased(reportes$prmt5_interaction$sr)
 a  <- a[a$cluster.fdr < fdrb,]
 a  <- a[!duplicated(a$junction),"junction"] 
 #is.AB[,table(ja=is.AB$ja==1,nonInt=!(is.AB$J3%in%a))]
 ancABok <- is.AB[is.AB$ja==1 & !(is.AB$J3%in%a),as.character(region)]

 #these are all the regions passing the interaction filter
 ABok <- unique(c(binABok,locABok,ancABok))

 table(is.AB[,as.character(region)%in%ABok])

  # Let's VennDiagram the discoveries of experiment A, experiment B and consolidatedAB(filtering out mild interaction effects)
  x.AB <- is.AB[as.character(region)%in%ABok,]
  u <- unique(c(is.A[,as.character(region)],is.B[,as.character(region)],x.AB[,as.character(region)]))
  vennDiagram(cbind(A=u%in%is.A[,as.character(region)],
                    B=u%in%is.B[,as.character(region)],
                    AB=u%in%x.AB[,as.character(region)]),main=paste("fdr.interaction >",fdrb[i]),cex=0.8)
}


# PCR validation (Table 4)
{
  
  #regions with splicing alterations according to consolidated dataset AB
  as.AB<- is.AB[as.character(region)%in%ABok,unique(as.character(locus))] 
  
  #Sanchez
  nOK<-c()
  lRNAseq2<-c(list(as.AB=as.AB),lRNAseq)
  for(i in seq_along(lRNAseq2)){
    agreed <- pcr$locus%in%lRNAseq2[[i]] 
    nOK<-c(nOK,sum(agreed))
  }
  names(nOK)<-c("ASpli.AB","ASpli.A","ASpli.B","LeafCutter.A","LeafCutter.B","MAJIQ.A","MAJIQ.B","rMATS.A","rMATS.B")
  print(nOK)
  
  
  #Xiang paper PCR validation: 12 events
  unlist(lapply(lRNAseq2,function(x){
    length(intersect(x,xiang))
  }))
}
