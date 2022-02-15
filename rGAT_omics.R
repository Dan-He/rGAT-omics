extract_candidates <- function(args){
  gwas_file<-args[1]; flanking<-as.numeric(args[2])

  gwas<-read.delim(gwas_file,as.is=T)
  if( all(is.element(c("SNP","Chr","Pos_hg19"),colnames(gwas))) )                ### check column names
    gwas<-gwas[,c("SNP","Chr","Pos_hg19")]  else stop("Wrong column names!\n")
  if(substr(gwas$Chr[1],1,3)!="chr")  gwas$Chr<-paste("chr",gwas$Chr,sep="")     ### add prefix "chr" in column "Chr"

  if( flanking < 0 ) stop("Please verify the flank region number!\n")              ### check flanking region number

  gene_path<-"All_human_genes"
  gene<-read.delim(gene_path,as.is=T)

  output<-list()
  for(i in 1:nrow(gwas))
  {
    start <- max(gwas[i,]$Pos_hg19 - flanking,0)
    end <- gwas[i,]$Pos_hg19 + flanking
    temp <- gene[gene$chrom == gwas[i,]$Chr,]
    index1 <- temp$start_hg19 < end & temp$start_hg19 > start
    index2 <- temp$end_hg19 < end & temp$end_hg19 > start
    temp <- temp[index1|index2,]

    if(nrow(temp)>0)
    {
      temp$SNP <- gwas[i,]$SNP
      temp$SNP_chr <- gwas[i,]$Chr
      temp$SNP_pos_hg19 <- gwas[i,]$Pos_hg19
      output <- rbind(output,temp)
    }
  }
  output
}
extract_evi <- function(gene){
  gene$Name <- substr(gene$Name,1,15)
  # collect evidence in de novo mutations-----------------------
  path <- "SCZ_DNM/"
  dnm <- read.delim(paste(path,"2014_Feb_De_novo_mutation_in_SCZ_Nature_623_trios_s3.txt",sep=""),as.is=T)
  dnm <- dnm[,c("Child.phenotype","Study","Genes","Gene.annotations")]
  fromer <- read.delim(paste(path,"2014_Feb_De_novo_mutation_in_SCZ_Nature_623_trios_s2.txt",sep=""),as.is=T)
  fromer$Child.phenotype <- "SZ"; fromer$Study <- "Fromer"
  fromer <- fromer[,c("Child.phenotype","Study","Genes","Gene.annotations")]
  gir <- dnm[dnm$Study=="Girard",]
  xu <- dnm[dnm$Study=="Xu",]
  gul <- dnm[dnm$Study=="Gulsuner",]
  case <- rbind(gir,xu[xu[,1]=="SZ",],gul[gul[,1]=="SZ",],fromer)
  index2 <- is.element(case$Gene.annotations,c("esplice","frameshift","nonsense","missense","codon-deletion","code-insertion"))
  case <- case[index2,]
  p1 <- length(unique(case$Genes))/20000 # suppose there were 20000 unqiue genes in human
  gene$dnm <- p1
  gene$dnm[!is.element(gene$official_name,case$Genes)] <- 1-p1

  # collect evidence in gross deletion------------
  gros_data <- read.csv("HGMD-SCZMut/schizophrenia_grosdel.csv",header = T,stringsAsFactors = F)
  gros_gene <- unique(gros_data$gene)
  p2 <- length(gros_gene)/20000
  gene$gross <- p2
  gene$gross[!is.element(gene$official_name,gros_gene)] <- 1-p2

  # collect evidence in gene differentially expression analysis----------
  dexpr<-read.delim("SCZ_DE/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-differentialExpression-includeAncestry-DxSCZ-DE.tsv")
  gene$de <- dexpr[match(gene$Name,dexpr$genes),]$P.Value
  cat("End of collecting evidence in category 1\n")

  # collect distance between genes and index SNPs-------------
  index<-gene$strand=="+"
  tss_dist<-0
  tss_dist[index] <- abs(gene[index,]$start_hg19-gene[index,]$SNP_pos_hg19)
  tss_dist[!index] <- abs(gene[!index,]$end_hg19-gene[!index,]$SNP_pos_hg19)
  gene$dist <- tss_dist
  # collect HiC links from capture HiC------------
  cap4<-read.delim("capHiC/GM12878_DRE_number",as.is=T)
  gene$cap4_enhancer_no<-cap4[match(gene$Name,cap4$Name),]$cap4_enhancer_no
  # DRE promoters in the cortical and subcortical plate and the germinal zone
  cp<-read.delim("BrainHiC/S22_TSS_CP.txt",as.is=T)
  gz<-read.delim("BrainHiC/S23_TSS_GZ.txt",as.is=T)
  cp$enhancer<-paste(cp$chr,cp$interacting_bin_start,cp$interacting_bin_end,sep=":")
  gz$enhancer<-paste(gz$chr,gz$interacting_bin_start,gz$interacting_bin_end,sep=":")

  gene$brain_cp<-0; gene$brain_gz<-0
  for(i in 1:nrow(gene))
  {
    temp<-cp[cp$ENSGID_for_TSS==gene[i,]$Name,]
    if(nrow(temp)>0)
      gene[i,]$brain_cp<-length(unique(temp$enhancer))
  }

  for(i in 1:nrow(gene))
  {
    temp<-gz[gz$ENSGID_for_TSS==gene[i,]$Name,]
    if(nrow(temp)>0)
      gene[i,]$brain_gz<-length(unique(temp$enhancer))
  }
  # collect HiC links from FANTOM5-------------
  eep<-read.delim("Fantom5/enhancer_tss_associations.bed.txt",as.is=T)
  eep<-eep[unlist(lapply(strsplit(eep$name,";"),function(x) length(x)))==5,]
  pos<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[1]]))
  chr<-unlist(lapply(strsplit(pos,":"),function(x) x[[1]]))
  start<-as.numeric(lapply(strsplit(unlist(lapply(strsplit(pos,":"),function(x)x[[2]])),"-"),function(x) x[[1]]))
  end<-as.numeric(lapply(strsplit(unlist(lapply(strsplit(pos,":"),function(x)x[[2]])),"-"),function(x) x[[2]]))
  fdr<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[length(x)]]))
  fdr<-as.numeric(lapply(strsplit(fdr,":"),function(x) x[[2]]))
  genes<-unlist(lapply(strsplit(eep$name,";"),function(x) x[[length(x)-2]]))

  eep$chr<-chr; eep$start<-start; eep$end<-end; eep$fdr<-fdr; eep$gene<-genes
  eep <- eep[,13:17];  eep <- eep[eep$fdr<1,]
  eep$enhancer<-paste(eep$chr,eep$start,eep$end,sep=":")

  gene$fantom5_enhancer_no<-0
  for(i in 1:nrow(gene))
  {
    temp<-eep[eep$gene==gene[i,]$official_name,]

    if(nrow(temp)>0)
      gene[i,]$fantom5_enhancer_no<-length(unique(temp$enhancer))
  }
  # collect gene expression value from BrainSpan------------
  rpkm <- read.table("BrainSpan/RPKM_in_BrainSpan.txt",header = T,sep = "\t")
  map_rpkm <- rpkm[match(gene$Name,rownames(rpkm)),c(16:25)] # expression in adolescence and adulthood
  gene <- cbind(gene,map_rpkm)
  cat("End of collecting evidence in category 2\n")

  # Deal with missing values------------
  library("sampling")
  library("kknn")
  set.seed(2019)
  train <- gene[!is.na(gene$de),c("dist","de","cap4_enhancer_no","fantom5_enhancer_no","brain_cp","brain_gz")]
  test <- gene[is.na(gene$de),c("dist","de","cap4_enhancer_no","fantom5_enhancer_no","brain_cp","brain_gz")]
  knn.fit <- kknn(de~.,train = train,test = test,k = 5,kernel = "triangular") # use knn to handle NAs
  gene[is.na(gene$de),"de"] <- knn.fit$fitted.values
  set.seed(123)
  rpkm.ind <- which(is.na(gene$adolescence_parietal.lobe),arr.ind = T)
  col_ind <- which(colnames(gene)%in%c("dist","cap4_enhancer_no","fantom5_enhancer_no","brain_cp","brain_gz","de"),arr.ind = T)
  n <- which(colnames(gene)=="adolescence_parietal.lobe",arr.ind = T)
  for(i in n:(n+9)){
    data <- gene[,c(col_ind,i)]
    colnames(data)[7] <- "y"
    trainc <- data[-rpkm.ind,]
    testc <- data[rpkm.ind,]
    fit <- kknn(y~.,train = trainc,test = testc,kernel = "triangular")
    gene[rpkm.ind,i] <- fit$fitted.values
  }
  gene
}

rGAT_omics <- function(args) {
  load("GO_BioGrid_HIPP.RData")
  nodes <- colnames(Net)
  gene <- extract_candidates(args)
  gene <- gene[(!is.na(gene$official_name)),]
  gene <- gene[is.element(gene$official_name,nodes),]
  gene <- extract_evi(gene)
  Net <- Net[,(is.element(nodes,unique(gene$official_name)))]
  category1 <- gene[,c("de","dnm","gross")]
  n <- which(colnames(gene)=="adolescence_parietal.lobe",arr.ind = T)
  category2 <- gene[,c("dist","cap4_enhancer_no","fantom5_enhancer_no","brain_cp","brain_gz")]
  category2 <- cbind(category2,gene[,c(n:(n+9))])
  burnin_round <- 3000
  after_burnin_round <- 3000
  exclude_samegene <- T
  region <- split(gene$official_name,gene$SNP)
  # collect extra weight of candidate genes
  cat("Collecting extra weight of candidate gene...\n")
  library(MASS)
  library(moments)
  # hotelling transform
  data <- category2
  C <- cov(data)
  M <- colMeans(data)
  a.e <- eigen(C,symmetric = T)
  V <- a.e$vectors
  m <- round(0.6*ncol(category2))
  U <- V[1:m,]
  data_h <- U%*%(apply(data,1,function(x) x-M))
  data_h <- t(data_h)
  # box-cox transform
  for(i in 1:ncol(data_h)){
    lam <- -floor(min(c(data_h[,i],0)))
    y <- data_h[,i] + lam
    b <- boxcox(y~1)
    lambda <- b$x[which.max(b$y)]
    if(lambda!=0){
        data_h[,i] <- y^lambda
      }else{
        data_h[,i] <- log(y)
      }
  }
  data_h <- scale(data_h)
  data_p <- apply(data_h,2,function(x) unlist(sapply(x,pnorm)))
  extra_p <- cbind(data_p,category1)
  extra_weight<-(-log(apply(extra_p,1,function(x) { pchisq(-2*sum(log(x)),df=2*length(x),lower.tail=F)})))
  extra_weight <- data.frame(gene=gene$official_name,extra_weight=extra_weight)
  extra_weight <- extra_weight[match(unique(extra_weight$gene),extra_weight$gene),]
  ####================= end of loading data ===============####

  ####================= burn in step ======================####

  #t0<-proc.time()
  set.seed(123)
  thres <- 0.01; pickup <- 0;
  num_region <- length(region)
  circle <- 1; chosen <- NULL
  remaining <- unlist(lapply(region,function(x) sample(x,1)))
  num0 <- rep(0,sum(unlist(lapply(region,length))))

  dif <- thres + 1; dif_record <- NULL
  while(dif > thres && circle<(burnin_round + 1))
  {
    pickup <- pickup %% num_region + 1
    if(pickup == 1)
      if(!is.null(chosen))
      {
        ###================================= calculate frequency =========================###
        num1 <- NULL
        for(j in 1:length(region))
          num1 <- c(num1,unlist(lapply(region[[j]],function(x) sum(chosen[,j]==x))))
        num1 <- num1 + num0
        if(circle>1)
        {
          freq0 <- num0/(num_region * (circle-1))
          freq1 <- num1/(num_region * circle)
          dif <- (sum((freq0-freq1)^2))^0.5
          if( circle%%50==0 )
          {
            cat("Burnin sampling, sampling circle:",circle,"\n")
          }
          dif_record<-c(dif_record,dif)
        }
        num0<-num1; chosen<-NULL; circle<-circle+1
      }

    pickup_p<-Net[,is.element(colnames(Net),remaining[-pickup])]
    pickup_p<-pickup_p[is.element(rownames(pickup_p),unlist(region[pickup])),]

    if(exclude_samegene && !is.null(dim(pickup_p)))                              ### if there is overlap between candidate genes and
    {                                                                            ### conditional genes, exclude the same genes
      pickup_p<-pickup_p[,!is.element(colnames(pickup_p),rownames(pickup_p))]
      if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
    }
    if(is.null(dim(pickup_p))) { pickup_p<-1; names(pickup_p) <- region[[pickup]] } else    ### when there is only one candiate gene
    { pickup_p <- apply(pickup_p,1,sum)
    pickup_p <- extra_weight[match(names(pickup_p),extra_weight$gene),]$extra_weight * pickup_p }

    if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i] <- 1/length(pickup_p)

    remaining[pickup] <- sample(names(pickup_p),1,replace=T,prob=pickup_p)
    chosen<-rbind(chosen,remaining)
  }

  ###===================== end of burn in step ===============================###

  ###======================= post-burn in step ===================================###

  pickup <- 0; num_region <- length(region); circle<-1; chosen <- NULL
  num0<-rep(0,sum(unlist(lapply(region,length))))

  joi_dis<-matrix(0,nrow=nrow(gene),ncol=nrow(gene))
  temp <- NULL;
  for(j in 1:length(region))
    temp <- c(temp,paste(names(region[j]),region[[j]],sep="_"))
  colnames(joi_dis) <- temp; rownames(joi_dis) <- temp

  thres <- 0.01; dif <- thres + 1
  while(dif>thres && circle<(after_burnin_round + 1) )
  {
    pickup <- pickup %% num_region + 1
    if(pickup==1)
      if(!is.null(chosen))
      {
        ###================================= calculate frequency =========================###
        num1<-NULL
        for(j in 1:length(region))
          num1<-c(num1,unlist(lapply(region[[j]],function(x) sum(chosen[,j]==x))))
        num1 <- num1 + num0
        if(circle>1)
        {
          freq0 <- num0 / (num_region * (circle - 1))
          freq1 <- num1 / (num_region * circle)
          dif <- (sum((freq0 - freq1)^2))^0.5
          if( circle%%50 == 0 )
          {
            cat("Post-burnin sampling, sampling circle:",circle,"\n")
          }
        }
        num0<-num1; circle<-circle+1; chosen<-NULL
        ###============================= end of calculating frequency =======================###
      }

    pickup_p<-Net[,is.element(colnames(Net),remaining[-pickup])]
    pickup_p<-pickup_p[is.element(rownames(pickup_p),unlist(region[pickup])),]

    if(exclude_samegene && !is.null(dim(pickup_p)))                              ### if there is overlap between candidate genes and
    {                                                                            ### conditional genes, exclude the same genes
      pickup_p<-pickup_p[,!is.element(colnames(pickup_p),rownames(pickup_p))]
      if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
    }

    if(is.null(dim(pickup_p))) { pickup_p <- 1; names(pickup_p)<-region[[pickup]] } else    ### when there is only one candiate gene
    { pickup_p<-apply(pickup_p,1,sum); pickup_p <- extra_weight[match(names(pickup_p),extra_weight$gene),]$extra_weight * pickup_p }

    if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i] <- 1/length(pickup_p)     ### avoid probabilit=0

    remaining[pickup] <- sample(names(pickup_p),1,replace=T,prob=pickup_p)
    chosen<-rbind(chosen,remaining)

    ###=============================== calculating joint distribution  =====================###
    index_col <- match(paste(names(remaining[-pickup]),remaining[-pickup],sep="_"),colnames(joi_dis))
    index_row <- match(paste(names(remaining[pickup]),remaining[pickup],sep="_"),colnames(joi_dis))
    joi_dis[index_row,index_col] <- joi_dis[index_row,index_col] + 1
  }

  ###=====================  end of post-burnin step  =======================================###

  ###===================== summarize and record the results ================================###
  freq<-cbind(unlist(region),freq1)

  region_indicator<-NULL
  gene_num<-as.numeric(lapply(region,length))
  for(i in 1:length(gene_num))
    region_indicator<-c(region_indicator,rep(names(region[i]),gene_num[i]))

  freq <- cbind(freq,region_indicator)
  colnames(freq) <- c("gene","post_prob","region")
  freq <- as.data.frame(freq,stringsAsFactors=F)
  freq[,2] <- as.numeric(freq[,2])

  output <- NULL                                      #### sort according to posterior probability
  for(i in unique(freq$region))
  {
    temp<-freq[freq$region==i,]
    output<-rbind(output,temp[order(temp$post_prob,decreasing=T),])
  }
  freq<-output

  cat("Recording the results!\n")
  freq
}

select_genes <- function(freq){
  SNP <- unique(freq$region)
  temp <- NULL
  for(snp in SNP){
    snpi <- freq[freq$region==snp,2]
    ind <- which(freq$region==snp,arr.ind = T)
    maxind <- which(snpi==max(snpi),arr.ind = T)
    temp <- c(temp,ind[maxind])
  }
  HRG <- freq[temp,]
  HRG <- HRG[match(unique(HRG$gene),HRG$gene),]
  LBG <- freq[freq$post_prob<median(freq$post_prob),]
  LBG <- LBG[match(unique(LBG$gene),LBG$gene),]
  result <- list(HRG=HRG,LBG=LBG)
  result
}

args <- c("./SNP_file/SCZ_108_loci",1000000)
