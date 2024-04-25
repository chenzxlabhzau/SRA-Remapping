library(data.table)
library(dplyr)
library(sva)
library(DESeq2)
library(optparse)  #调用optpares包

option_list <- list(make_option(opt_str = c("-c", "--count"), type = "character", default = FALSE, help = "raw count文件"),
                    make_option(opt_str = c("-m", "--metadata"), type = "character", default = FALSE, help = "metainfo文件"),
                    make_option(opt_str = c("-g", "--group"), type = "character", default = FALSE, help = "分组"),
                    make_option(opt_str = c("-o", "--output"), type = "character", default = FALSE, help = "输出"))
opt_parser = OptionParser(option_list=option_list)  #解析参数
opt = parse_args(opt_parser)   #解析参数 
x1 = opt$count   #将输入的-d参数传递给x
y1 = opt$metadata #将输入的-c参数传递给y
group = opt$group
output = opt$output

count<-fread(x1)%>%as.data.frame()
rownames(count)=count$FBgn
count<-count[,-1]
metadata<-read.table(y1)


coldata <- metadata[metadata$group==group,]%>%filter(.,sex=="male"|sex=="female")
coldata$sex<-factor(coldata$sex,levels = c("male","female"))
n_male=count(coldata$sex=="male")
n_female=count(coldata$sex=="female")
if (n_male<2|n_female<2) {
  print(paste(group,"have",count(coldata$sex=="male"),"males and",count(coldata$sex=="female"),"females"),sep = " ")
}else{
  countData<-count[,coldata$curr_SRX]
  if(length(unique(coldata$study))>1){dds<- DESeqDataSetFromMatrix(
    countData =countData,
    colData=coldata,
    design = ~study+sex)
  }else{
    dds<- DESeqDataSetFromMatrix(
      countData =countData,
      colData=coldata,
      design = ~sex)}
  dds<- DESeq(dds)
  normalized_counts<-counts(dds,normalized=TRUE)
  normalized_counts<-normalized_counts[rowSums(normalized_counts)>1,]
  conditions=factor(t(coldata$sex),levels = c("male","female"))
  if (n_male>=8&n_female>=8){
    #remove batch effect
    if(length(unique(coldata$study))>1){
      mod = model.matrix(~as.factor(sex), data=coldata)
      normalized_counts<- ComBat(dat = as.matrix(log2(normalized_counts+1)), batch = coldata$study,mod = mod,par.prior = F)}
    # Run the Wilcoxon rank-sum test for each gene
    pvalues<-sapply(1:nrow(normalized_counts),function(j){
      data<-cbind.data.frame(gene=as.numeric(t(normalized_counts[j,])),conditions)
      p=wilcox.test(gene~conditions, data)$p.value
      return(p)})
    fdr=p.adjust(pvalues,method = "fdr")
    # Calculate fold-change for each gene
    conditionsLevel<-levels(conditions)
    dataCon1=normalized_counts[,c(which(conditions==conditionsLevel[1]))]
    dataCon2=normalized_counts[,c(which(conditions==conditionsLevel[2]))]
    foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))
    outRst<-data.frame(log2FoldChange=foldChanges, pvalue=pvalues, padj=fdr)
    rownames(outRst)=rownames(normalized_counts)
    res<-na.omit(outRst)
  }else{
    res<- as.data.frame(results(dds,contrast=c("sex","male","female")))
  }
}
res$gene_id=rownames(res)
group<-gsub(" ","_",group)
write.table(res,output,row.names = FALSE)


