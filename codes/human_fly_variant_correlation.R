#packages
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

metadata<-read.csv("./sra_fomal/data/metadata.csv")
debatched_fpkm<-fread("./sra_fomal/data/debatched_fpkm")%>%as.data.frame()
rownames(debatched_fpkm)<-debatched_fpkm$V1
debatched_fpkm<-debatched_fpkm[,-1]
adult_metadata<-subset(metadata,metadata$broad_stage=="adult")
adult_metadata<-adult_metadata[grep("brain|heart|gut|muscle|ovary|malpighian|testis|accessory|fat",adult_metadata$broad_tissue),]
rownames(adult_metadata)=adult_metadata$X
adult_metadata<-adult_metadata[,-c(1,9)]
debatched_fpkm<-debatched_fpkm[,adult_metadata$curr_SRX]
fly_expressed<-as.data.frame(apply(debatched_fpkm, 1, function(x){sum(x>1)}))
colnames(fly_expressed)="expressed_num"
fly_expressed<-subset(fly_expressed,fly_expressed$expressed_num>1)
debatched_fpkm<-subset(debatched_fpkm,rownames(debatched_fpkm)%in%rownames(fly_expressed))
pmbc<-CreateSeuratObject(counts=debatched_fpkm,project="rna_seq")
pmbc@meta.data$tissue<-adult_metadata$broad_tissue
pmbc<-FindVariableFeatures(pmbc,selection.method="vst",nfeatures = 2000)
#human-fly
fly_variant<-pmbc@assays$RNA@meta.features

human_all_fpkm<-read.table("/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq/Human.allgene.fpkm.txt")
tissue_info<-data.frame(human_sample=colnames(human_all_fpkm))
tissue_info$tissue<-unlist(lapply(strsplit(tissue_info$human_sample,"_"), "[[",1))
tissue_info<-tissue_info[grep("adipose|cerebral|small|heart|kidney|ovary|prostate|smooth|testis",tissue_info$tissue),]
human_all_fpkm<-human_all_fpkm[,tissue_info$human_sample]
human_expressed<-as.data.frame(apply(human_all_fpkm, 1, function(x){sum(x>1)}))
colnames(human_expressed)="expressed_num"
human_expressed<-subset(human_expressed,human_expressed$expressed_num>1)
human_all_fpkm<-subset(human_all_fpkm,rownames(human_all_fpkm)%in%rownames(human_expressed))
#seurat
pmbc<-CreateSeuratObject(counts=human_all_fpkm,project="rna_seq")
pmbc@meta.data$tissue<-tissue_info$tissue

#find high-variable gene
pmbc<-FindVariableFeatures(pmbc,selection.method="vst",nfeatures = 2000)

#human-fly
human_variant<-pmbc@assays$RNA@meta.features

fly_human_homo<-read.table("./sra_fomal/data/mart_export.txt",sep = ",",header = TRUE)
fly_human_one2one<-unique(subset(fly_human_homo,fly_human_homo$Drosophila.melanogaster.homology.type=="ortholog_one2one")[,c(1,6)])
fly_human_one2one$human_variant<-human_variant[fly_human_one2one$Gene.stable.ID,"vst.variance.standardized"]
fly_human_one2one$fly_variant<-fly_variant[fly_human_one2one$Drosophila.melanogaster.gene.stable.ID,"vst.variance.standardized"]

human_essential<-read.csv("./sra_fomal/data/Phenotypic.classes.csv",row.names = 2)
fly_human_one2one$human_essential<-human_essential[fly_human_one2one$Gene.stable.ID,"Pheno_class"]
fly_human_one2one_naomit<-na.omit(fly_human_one2one)

fly_human_one2one_naomit%>%group_by(human_essential)%>%summarise(cor=cor.test(human_variant,fly_variant,method = "spearman")$estimate,p=cor.test(human_variant,fly_variant,method = "spearman")$p.value)
cor1<-cor.test(lethal0$human_variant,lethal0$fly_variant,method = "spearman")
cor2<-cor.test(lethal_disease1$human_variant,lethal_disease1$fly_variant,method = "spearman")
cor3<-cor.test(Disease2$human_variant,Disease2$fly_variant,method = "spearman")
cor4<-cor.test(OtherGenes3$human_variant,OtherGenes3$fly_variant,method = "spearman")








