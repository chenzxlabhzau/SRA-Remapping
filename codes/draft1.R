#sex bias analysis draft
library(dplyr)
library(sva)
metadata$group=paste(metadata$broad_stage,metadata$broad_tissue,sep = "_")

group<-metadata%>%group_by(broad_stage,broad_tissue,sex,group)%>%summarise(n=n())%>%
  na.omit()%>%filter(.,sex=="male"|sex=="female")

group2<-metadata%>%group_by(broad_stage,broad_tissue,sex,group)%>%summarise(n=n())%>%
  na.omit()%>%filter(.,sex=="male"|sex=="female")%>%
  group_by(broad_stage,broad_tissue,group)%>%summarise(sex_n=n())%>%filter(sex_n==2)

count<-fread("./data/GSE117217_gene_counts.tsv")%>%as.data.frame()
rownames(count)=count$FBgn
count<-count[,-1]
whole_body_info<-subset(group2,group2$broad_tissue=="whole body")
tissue_info<-subset(group2,group2$broad_tissue!="whole body")
library(DESeq2)
res <- purrr::map(tissue_info$group,function(i){
  coldata <- metadata[metadata$group==i,]%>%filter(.,sex=="male"|sex=="female")
  coldata$sex<-factor(coldata$sex,levels = c("male","female"))
  n_male=count(coldata$sex=="male")
  n_female=count(coldata$sex=="female")
  if (n_male<2|n_female<2) {
    print(paste(i,"have",count(coldata$sex=="male"),"males and",count(coldata$sex=="female"),"females"),sep = " ")
    next()
  }else{
    countData<-count[,coldata$curr_SRX]
    dds<- DESeqDataSetFromMatrix(
    countData =countData,
    colData=coldata,
    design = ~study+sex)
    dds<- DESeq(dds)
    normalized_counts<-counts(dds,normalized=TRUE)
    normalized_counts<-normalized_counts[rowSums(normalized_counts)>1,]
    conditions=factor(t(coldata$sex),levels = c("male","female"))
    if (n_male>=8&n_female>=8) {
      #remove batch effect
      mod = model.matrix(~as.factor(sex), data=coldata)
      combat_normalized_count<- ComBat(dat = as.matrix(log2(normalized_counts+1)), batch = coldata$study,mod = mod,par.prior = F)
      # Run the Wilcoxon rank-sum test for each gene
      pvalues<-sapply(1:nrow(combat_normalized_count),function(j){
        data<-cbind.data.frame(gene=as.numeric(t(combat_normalized_count[j,])),conditions)
        p=wilcox.test(gene~conditions, data)$p.value
        return(p)})
      fdr=p.adjust(pvalues,method = "fdr")
      # Calculate fold-change for each gene
      conditionsLevel<-levels(conditions)
      dataCon1=normalized_counts[,c(which(conditions==conditionsLevel[1]))]
      dataCon2=normalized_counts[,c(which(conditions==conditionsLevel[2]))]
      foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
      outRst<-data.frame(log2FoldChange=foldChanges, pvalue=pvalues, padj=fdr)
      rownames(outRst)=rownames(normalized_counts)
      outRst=na.omit(outRst)
      return(outRst)
    }else{
      outdds <- as.data.frame(results(dds,contrast = c("male","female")))
      return(outdds)
    }
  }
})
names(res)=paste(info_summary$broad_stage,info_summary$broad_tissue,sep = "_")
