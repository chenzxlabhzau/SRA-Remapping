# code was written on 2022/04/14
# set dir
setwd("~/sra_fomal/")

#colorset:

#"egg" = "#C82423"
#"embryo"="#FF8884"
#"larva"="#F8AC8C"
#"pupa" = "#9AC9DB"
#"adult" = "#2878B5"

#tissue&cell_line (sort as sample size)
# "#b71c1c","#c62828","#ef5350","#ef9a9a","#880e4f","#c2185b","#ec407a","#f06292","#f48fb1","#f8bbd0",
# "#4a148c","#7b1fa2","#ab47bc","#ce93d8","#311b92","#673ab7","#7e57c2","#b39ddb","#1a237e","#303f9f",
# "#3f51b5","#c5cae9","#0d47a1","#1976d2","#2196f3","#90caf9","#0277bd","#039be5","#4fc3f7","#b3e5fc",
# "#006064","#0097a7","#00bcd4","#80deea","#004d40","#00897b","#4db6ac","#b2dfdb","#1b5e20","#43a047",
# "#81c784","#c8e6c9","#558b2f","#7cb342","#9ccc65","#c5e1a5","#827717","#afb42b","#d4e157","#e6ee9c",
# "#f57f17","#f9a825","#fdd835","#ffeb3b","#fff59d","#ff6f00","#ffa000","#ffc107","#ffd54f","#ffecb3",
# "#e65100","#ef6c00","#ff9800","#ffb74d","#ffe0b2","#bf360c","#e64a19","#ff7043","#ff8a65","#ffab91",
# "#ffccbc","#3e2723","#5d4037","#8d6e63","#bcaaa4","#424242","#757575","#bdbdbd","#eeeeee","#263238",
# "#455a64","#607d8b","#78909c","#b0bec5"

####################################
####### read count matrix   ########
####################################

#import packages
library(data.table)
library(dplyr)
#read count matrix
gene_count<-fread("~/sra_fomal/data/GSE117217_gene_counts.tsv") %>% as.data.frame()
rownames(gene_count)=gene_count[,1]
gene_count<-gene_count[,-1]

############################################
####### retained counts > 2000000   ########
############################################

#filter genes by read counts
gene_hist_table<-data.frame(count=apply(gene_count,2,sum), sample=colnames(gene_count))#calculate reads count for each sample
retain_by_count<-gene_count[,gene_hist_table$sample[gene_hist_table$count>2000000]]#remove samples with <2000000 reads

#clean
rm(gene_hist_table)

############################################
############# calculate fpkm  ##############
############################################

# prepare gene length infomation
# exon length from "dmel_r6-11.gtf", codes in linux is as follows:
# work dir:/home/wangdy/sra_fomal/data
# cat dmel_r6-11.gtf |awk '$3=="exon"{print $5-$4"\t"$10}'|sed "s/;//g"|sed "s/\"//g"|awk '{print $2"\t"$1}'>dmel_r6.11.exonlength

# aquire gene exon length
gene_exon_length<-read.table("~/sra_fomal/data/dmel_r6.11.exonlength",col.names = c("gene_id","length"))
gene_exon_length<-gene_exon_length%>%group_by(gene_id)%>%summarise(length=sum(length))%>%as.data.frame()
rownames(gene_exon_length)=gene_exon_length$gene_id
gene_id<-intersect(gene_exon_length$gene_id,rownames(retain_by_count)) #"FBgn0027554" not in annotation
retain_by_count<-retain_by_count[gene_id,]
gene_exon_length<-gene_exon_length[gene_id,]

# calculate fpkm
totalcounts<-colSums(retain_by_count)
fpkm <- t(do.call(rbind, lapply(1:length(totalcounts),function(i){10^9*retain_by_count[,i]/gene_exon_length$length/totalcounts[i]}))) %>% as.data.frame()
colnames(fpkm)=colnames(retain_by_count)
rownames(fpkm)=rownames(retain_by_count)
write.table(fpkm,"/home/wangdy/sra_fomal/data/fpkm")

#clean
rm(list=ls())

#################################################
############# remove batch effect  ##############
#################################################

#read in fpkm table
fpkm<-fread("~/sra_fomal/data/fpkm") %>% as.data.frame()
rownames(fpkm)=fpkm$V1
fpkm<-fpkm[,-1]

#read in batch infomation
batch<-read.csv("~/sra_fomal/data/metadata.csv",stringsAsFactors = FALSE,row.names = 1)

#Keep the samples of the two files consistent
fpkm<-fpkm[,batch$curr_SRX]

###Performing batch correction whthin tissues and stages

#import packages
library(sva)

#debatch
batch$group=paste(paste(batch$cell_line,batch$broad_stage,sep = ""),batch$broad_tissue,sep = "_")
debatched_fpkm<-data.frame(row.names = rownames(fpkm))
batch_summary<-batch%>%group_by(group)%>%summarise(batch_n=length(unique(study)))
debatch_groups<-as.data.frame(subset(batch_summary,batch_summary$batch_n!=1))
for (i in 1:nrow(debatch_groups)) {
  temp_batch<-subset(batch,batch$group==debatch_groups[i,"group"])
  temp_fpkm<-fpkm[,temp_batch$curr_SRX]
  unfilter_df = ComBat(temp_fpkm, batch = temp_batch$study)
  debatched_fpkm<-cbind(debatched_fpkm,unfilter_df)
}
debatched_fpkm<-cbind(debatched_fpkm,fpkm[,setdiff(batch$curr_SRX,colnames(debatched_fpkm))])
write.table(debatched_fpkm,"/home/wangdy/sra_fomal/data/debatched_fpkm")
#clean
rm(list = ls())

##########################################
############# data summary  ##############
##########################################

#read in metadata
metadata<-read.csv("~/sra_fomal/data//metadata.csv",stringsAsFactors = FALSE,row.names = 1)

#factor level
tissue_stage<-subset(metadata,is.na(metadata$cell_line)==T)
cell_line<-subset(metadata,is.na(metadata$cell_line)==F)
tissue_level<-tissue_stage%>%group_by(broad_tissue)%>%summarise(num = n())
tissue_level<-tissue_level[order(tissue_level$num,decreasing = T),]
cl_level<-cell_line%>%group_by(cell_line)%>%summarise(num = n())
cl_level<-cl_level[order(cl_level$num),]
colnames(tissue_level)=c("level","n")
colnames(cl_level)=c("level","n")
levels<-rbind(tissue_level,cl_level)


#prepare plot table
metadata$group=gsub("NA","",paste(metadata$broad_tissue,metadata$cell_line,sep = ""))
metadata[grep("RNA",metadata$cell_line),"group"]="Ras[V12]-wts[RNAi]"
plot_table<-metadata%>%group_by(group,broad_stage)%>%summarise(stage_n=n())
plot_table$broad_stage[which(is.na(plot_table$broad_stage)==T)]="cell line"
plot_table$broad_stage<-factor(plot_table$broad_stage,levels = c("adult","pupa","larva","embryo","egg","cell line"))
group_sum<-plot_table%>%group_by(group)%>%summarise(sum=sum(stage_n))
plot_table<-merge(plot_table,group_sum,by="group")
plot_table$group<-factor(plot_table$group,levels = levels$level)
levels$level<-factor(levels$level,levels=levels$level)

#plot
empty_bar <- 0
levels$id <-seq(1,nrow(levels))
number_of_bar <- nrow(levels)
levels$angle<-90-360*(levels$id-0.5)/number_of_bar    
levels$hjust<-ifelse((levels$angle<(-90)),1,0)
levels$angle<-ifelse(levels$angle<(-90),levels$angle+180,levels$angle)
metadata_plot<-ggplot() + 
  geom_bar(data = plot_table, aes(x=group, y=-stage_n/sum, fill=broad_stage),stat = "identity")+
  scale_fill_discrete(limits = c("egg", "embryo", "larva","pupa","adult"),guide=FALSE)+
  scale_fill_manual(values=c("egg" = "#fbefc4","embryo"="#52696f","larva"="#e57b7f","pupa" = "#9e3150", "adult" = "#88baa4","cell line"="white"))+
  # guides(fill=guide_legend(override.aes = list(size=0.01)))+
  # theme(legend.text = element_text(size = 2))+
  # theme(legend.key.size = unit(0.1,"inches"))+
  new_scale('fill')+
  geom_bar(data = levels, aes(x=level,y=log10(n+0.1),fill=level),stat = "identity")+
  scale_fill_discrete(limits=as.character(levels$level))+
  scale_fill_manual(values = c("#b71c1c","#c62828","#ef5350","#ef9a9a","#880e4f","#c2185b","#ec407a","#f06292","#f48fb1","#f8bbd0",
                               "#4a148c","#7b1fa2","#ab47bc","#ce93d8","#311b92","#673ab7","#7e57c2","#b39ddb","#1a237e","#303f9f",
                               "#3f51b5","#c5cae9","#0d47a1","#1976d2","#2196f3","#90caf9","#0277bd","#039be5","#4fc3f7","#b3e5fc",
                               "#006064","#0097a7","#00bcd4","#80deea","#004d40","#00897b","#4db6ac","#b2dfdb","#1b5e20","#43a047",
                               "#81c784","#c8e6c9","#558b2f","#7cb342","#9ccc65","#c5e1a5","#827717","#afb42b","#d4e157","#e6ee9c",
                               "#f57f17","#f9a825","#fdd835","#ffeb3b","#fff59d","#ff6f00","#ffa000","#ffc107","#ffd54f","#ffecb3",
                               "#e65100","#ef6c00","#ff9800","#ffb74d","#ffe0b2","#bf360c","#e64a19","#ff7043","#ff8a65","#ffab91",
                               "#ffccbc","#3e2723","#5d4037","#8d6e63","#bcaaa4","#424242","#757575","#bdbdbd","#eeeeee","#263238",
                               "#455a64","#607d8b","#78909c","#b0bec5"))+
  coord_polar()+
  theme_minimal()+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_text(data=levels,aes(x=level, y=log10(n+0.1)+0.5,label=paste(level," (",n,")",sep = ""),hjust=hjust,angle=angle),color="black",size=2,fontface="plain",inherit.aes=FALSE)+
  ylab("log10 ( n + 0.1 )")+
  theme(legend.title = element_blank(),legend.key.width=unit(0.1,"cm"),legend.key.height=unit(0.1,"cm"),legend.position = c(0.466, 0.5),legend.text = element_text(size = 5))
ggsave("./results/plot/Rplot.pdf",metadata_plot)

#########################################################
############# seurat - dimension reduction ##############
#########################################################

#packages
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)

#prepare data
metadata<-read.csv("~/sra_fomal/data//metadata.csv",stringsAsFactors = FALSE,row.names = 1)
cl_rows<-which(is.na(metadata$cell_line)==FALSE)
metadata[cl_rows,"sex"]="cell line"
metadata[cl_rows,"broad_tissue"]="cell line"
metadata[cl_rows,"broad_stage"]="cell line"

debatched_fpkm<-fread("./sra_fomal/data/debatched_fpkm") %>% as.data.frame()
rownames(debatched_fpkm)=debatched_fpkm[,1]
debatched_fpkm<-debatched_fpkm[,-1]
debatched_fpkm<-debatched_fpkm[,metadata$curr_SRX]

# #fraction
# #total gene
# type_length_table<-read.table("./data/type_length_table")
# debatched_fpkm<-debatched_fpkm[intersect(rownames(debatched_fpkm),type_length_table$gene_id),]
# type_length_table<-na.omit(type_length_table[rownames(debatched_fpkm),])
# table(type_length_table$gene_type)
# total_gene=2402+13806+318
# total_gene #16526
# 
# #only consider protein-coding/lnc/pseudo gene
# filtered_type_length_table<-subset(type_length_table,type_length_table$gene_type=="protein-coding"|type_length_table$gene_type=="lnc"|type_length_table$gene_type=="pseudogene")
# protein_lnc_pseudo_fpkm<-debatched_fpkm[rownames(debatched_fpkm)%in%filtered_type_length_table$gene_id,]
# expressed_sample<-data.frame(expressed_sample=apply(debatched_fpkm,1,function(x){sum(x>1)}))
# table(expressed_sample$expressed_sample>(total_gene*0.5))
# 6455/total_gene
# #39.06%
# table(expressed_sample$expressed_sample<(total_gene*0.1))
# 4559/total_gene
# #27.59%


debatched_fpkm<-fread("./sra_fomal/data/debatched_fpkm") %>% as.data.frame()
rownames(debatched_fpkm)=debatched_fpkm[,1]
debatched_fpkm<-debatched_fpkm[,-1]

#seurat
pmbc<-CreateSeuratObject(counts=debatched_fpkm,project="rna_seq")
pmbc@meta.data$sex<-metadata$sex
pmbc@meta.data$broad_tissue<-metadata$broad_tissue
pmbc@meta.data$broad_stage<-metadata$broad_stage
pmbc@meta.data$study<-metadata$study


#find high-variable gene
pmbc<-FindVariableFeatures(pmbc,selection.method="vst",nfeatures = 3)
#high variable gene plot
top10 <- head(VariableFeatures(pmbc), 3)
#plot variable features with and without labels
plot1 <- VariableFeaturePlot(pmbc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("./sra_fomal/results/plot/003_dynamic_gene.pdf",width = 6,height = 4.5)
plot2
dev.off()

#Scaling the data
all.genes<-rownames(pmbc)
pmbc<-ScaleData(pmbc,features=all.genes)
#Perform linear dimensional reduction
pmbc<-RunPCA(pmbc,features=VariableFeatures(object=pmbc))
pmbc <- FindNeighbors(pmbc, dims = 1:50)
pmbc <- FindClusters(pmbc, resolution = 0.1)

#UMAP
pmbc <- RunUMAP(pmbc, label = TRUE, dims = 1:50)

FeaturePlot(object = pmbc, 
            reduction = "umap",
            features = "FBgn0033450",
            pt.size=0.6)+
  # scale_colour_gradientn(colours=c(c('azure4','yellow','orange',"red","#b20000")),
  #                        limits = c(0, 200),
  #                        breaks = c(0, 50,100,150,200),
  #                        labels = c(0, 50,100,150,200)) +  ## 停在0.175处
  scale_colour_gradient2(low = c('azure4','azure3','azure2',"azure"), mid = "azure",
                         high = c("yellow","coral3","red","#b20000"), midpoint = 20,  # midpoint =0.09在0.13处停  =0.1在0.145处停   在0.095在0.14处停
                         limits = c(0, 200),
                         breaks = c(0, 50,100,150,200),
                         labels = c(0, 50,100,150,200))+
  ggtitle(label = "cluster1 blueprint score")+
  theme(axis.title=element_text(size=40,face="bold"),    ## x,y坐标轴标题字体
        axis.text=element_text(vjust=1,size=40,face = "bold"),    ## x,y坐标轴刻度线字体
        # plot.title = element_text(size=40),
        plot.title=element_blank(),
        legend.text=element_text(size=35),    ## 图例字体
        legend.title = element_blank(),
        legend.key.size = unit(1.4,"line"))

pdf("./sra_fomal/results/plot/003_dy_FBgn0033450.pdf",width = 5,height = 4)
FeaturePlot(pmbc,"FBgn0033450",min.cutoff = 0, max.cutoff = 200,pt.size=0.4)#,min.cutoff = 0, max.cutoff = 2000
dev.off()
pdf("./sra_fomal/results/plot/003_cs_FBgn0000559.pdf",width = 5,height = 4)
FeaturePlot(pmbc,"FBgn0000559",min.cutoff = 0, max.cutoff =1000,pt.size=0.2)
dev.off()

#cluster by batch
pmbc@meta.data$seurat_clusters<-as.character(pmbc@meta.data$study)
pmbc@meta.data$seurat_clusters<-as.factor(pmbc@meta.data$study)
pmbc@active.ident<-pmbc@meta.data$seurat_clusters
debatch_umap<-DimPlot(pmbc, reduction = "umap",pt.size=0.1,label = FALSE)
pdf("./results/plot/002_debatch_umap.pdf",width = 50,height = 10)
debatch_umap
dev.off()

#clustered by sex
pmbc@meta.data$seurat_clusters<-as.character(pmbc@meta.data$sex)
pmbc@meta.data$seurat_clusters<-as.factor(pmbc@meta.data$sex)
pmbc@active.ident<-pmbc@meta.data$seurat_clusters
sex_umap<-DimPlot(pmbc, reduction = "umap",label = FALSE)
pdf("./results/plot/002_sex_umap.pdf",width = 7,height = 5)
sex_umap
dev.off()

#clustered by stage
pmbc@meta.data$seurat_clusters<-as.character(pmbc@meta.data$broad_stage)
pmbc@meta.data$seurat_clusters<-factor(pmbc@meta.data$broad_stage,levels=c("egg","embryo","larva","pupa","adult","cell_line"))
pmbc@active.ident<-pmbc@meta.data$seurat_clusters
stage_umap<-DimPlot(pmbc, reduction = "umap", label = FALSE)
pdf("./results/plot/002_stage_umap.pdf",width = 7,height = 5)
stage_umap
dev.off()

#clustered by tissue
pmbc@meta.data$seurat_clusters<-as.character(pmbc@meta.data$broad_tissue)
pmbc@meta.data$seurat_clusters<-as.factor(pmbc@meta.data$broad_tissue)
pmbc@active.ident<-pmbc@meta.data$seurat_clusters
tissue_umap<-DimPlot(pmbc, reduction = "umap", label = TRUE)
pdf("./results/plot/002_tissue_umap.pdf",width = 14,height = 10)
tissue_umap
dev.off()

#Extract the list of highly variable genes
meta_features<-pmbc[["RNA"]]@meta.features
write.table(meta_features,"/home/wangdy/sra_fomal/results/result_table/seurat_variant_table")

#########################################################
############# expression pattern - variant ##############
#########################################################

library(rstatix)
library(ggpubr)
library(tidyverse)
#dN/dS
variant_table<-read.table("./sra_fomal//results/result_table/seurat_variant_table")
variant_table$gene_id<-rownames(variant_table)
evo_info<-read.table("./sra_fomal/data/evo_info",header = TRUE)
variant_table<-merge(variant_table,evo_info,by="gene_id")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
variant_table<-subset(variant_table,variant_table$gene_id%in%rownames(expressed_detect))

# cor<-cor.test(variant_table$omega,variant_table$vst.variance.standardized,method = "spearman")
# ggplot(variant_table,mapping = aes(x = log2(omega), y = log2(vst.variance.standardized)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
#   theme_bw()+
#   ylab("Standardized expression variation (log2)")+
#   xlab("dn/ds (log2)")+
#   theme(axis.title.x=element_text(size=18), 
#         axis.title.y=element_text(size=16),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16))+
#   geom_smooth(method = "lm",colour="black")+
#   annotate("text", x=3.5, y=6, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)
# ggsave(filename = "./sra_fomal/results/plot/004_std_var_log2dnds.pdf",device = "pdf",width = 5,height = 4.5)
# #ggsave("~/sra/plot/Fig2_stage_age.pdf",width = 12,height = 6)
variant_table<-variant_table[order(variant_table$vst.variance.standardized,decreasing = TRUE),]
variant_table$group<-factor(c(rep(c(6,5,4,3,2,1),each=1731),1),levels = c(1,2,3,4,5,6))
# my_comparisons = list(c("1","2"),c("2","3"),c("3","4"),c("4","5"),c("5","6"))

df_p_val <- variant_table%>%
  wilcox_test(formula = omega ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position()%>%
  slice(1,6,10,13,15)%>%
  mutate(y.position=c(0.3,0.33,0.4,0.47,0.5))
ggplot(variant_table,aes(x=group,y=omega))+geom_boxplot(aes(fill=group),outlier.alpha = 0, notch=TRUE)+
  stat_pvalue_manual(data = df_p_val,label = '{p.adj}',label.size=5)+
  scale_fill_brewer(palette = "YlGnBu")+
  ylim(0,0.5)+
  ylab("dN/dS")+
  xlab("")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = 'none')
  
  #stat_compare_means(comparisons = my_comparisons,label.y = 0.5)+
  # geom_signif(comparisons = my_comparisons,
  #             map_signif_level=T,
  #             y_position = 0.3,
  #             textsize=6,test=wilcox.test,step_increase=0)
ggsave(filename = "./sra_fomal/results/plot/004_m1_dnds_variance.pdf",device = "pdf",width = 5,height = 4.5)


#promoter conservation
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
variant_table$gene_id<-rownames(variant_table)
promoter_conservation<-read.table("./sra_fomal/data/Dmel_BDGP_promoter_phylop_result",col.names = c("name","size","covered","sum","mean0","gene_mean_phylop"))[,c(1,6)]
variant_table<-merge(variant_table,promoter_conservation,by.x="gene_id",by.y="name")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
variant_table<-subset(variant_table,variant_table$gene_id%in%rownames(expressed_detect))
variant_table<-variant_table[order(variant_table$vst.variance.standardized,decreasing = TRUE),]
variant_table$group<-factor(c(rep(c(6,5,4,3,2),each=2756),rep(1,2755)),levels = c(1,2,3,4,5,6))

df_p_val <- variant_table%>%
  wilcox_test(formula = gene_mean_phylop ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position()%>%
  slice(1,6,10,13,15)%>%
  mutate(y.position=c(2.95,2.7,2.5,2.35,2.2))
ggplot(variant_table,aes(x=group,y=gene_mean_phylop))+geom_boxplot(aes(fill=group),outlier.alpha = 0, notch=TRUE)+
  stat_pvalue_manual(data = df_p_val,label = '{p.adj}',label.size=4.5,tip.length = 0)+
  scale_fill_brewer(palette = "YlGnBu")+
  ylim(0,3)+
  ylab("Promoter conservation (phyloP)")+
  xlab("")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = 'none')
ggsave(filename = "./sra_fomal/results/plot/004_m1_phylop_variance.pdf",device = "pdf",width = 5,height = 4.5)

# cor<-cor.test(variant_table$gene_mean_phylop,variant_table$vst.variance.standardized,method = "spearman")
# ggplot(variant_table,mapping = aes(x = log2(gene_mean_phylop), y = log2(vst.variance.standardized)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
#   xlim(c(-2.5,1.5))+
#   theme_bw()+
#   ylab("Standardized expression variation (log2)")+
#   xlab("Promoter conservation (phyloP, log2)")+
#   theme(axis.title.x=element_text(size=18), 
#         axis.title.y=element_text(size=16),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16))+
#   annotate("text", x=0.8, y=6.5, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)+
#   geom_smooth(method = "lm",colour="black")

#phastcons
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
variant_table$gene_id<-rownames(variant_table)
promoter_conservation<-read.table("/home/wangdy/sra/conservation/result/Dmel_BDGP_promoter_phastcons_result",col.names = c("name","size","covered","sum","mean0","gene_mean_phastcons"))[,c(1,6)]
variant_table<-merge(variant_table,promoter_conservation,by.x="gene_id",by.y="name")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
variant_table<-subset(variant_table,variant_table$gene_id%in%rownames(expressed_detect))
variant_table<-variant_table[order(variant_table$vst.variance.standardized,decreasing = TRUE),]
variant_table$group<-factor(c(rep(c(6,5,4,3,2),each=2756),rep(1,2755)),levels = c(1,2,3,4,5,6))

# cor<-cor.test(variant_table$gene_mean_phastcons,variant_table$vst.variance.standardized,method = "spearman")
# ggplot(variant_table,mapping = aes(x = gene_mean_phastcons, y = log2(vst.variance.standardized)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
#   theme_bw()+
#   ylab("Standardized expression variation (log2)")+
#   xlab("Promoter conservation (phastcons)")+
#   theme(axis.title.x=element_text(size=18), 
#         axis.title.y=element_text(size=16),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16))+
#   annotate("text", x=0.2, y=6.5, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)+
#   geom_smooth(method = "lm",colour="black")

df_p_val <- variant_table%>%
  wilcox_test(formula = gene_mean_phastcons ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position()%>%
  slice(1,6,10,13,15)%>%
  mutate(y.position=c(1,1,1,1,1))
ggplot(variant_table,aes(x=group,y=gene_mean_phastcons))+geom_boxplot(aes(fill=group),outlier.alpha = 0, notch=TRUE)+
  stat_pvalue_manual(data = df_p_val,label = '{p.adj}',label.size=4,tip.length = 0,bracket.shorten = 0.2)+
  scale_fill_brewer(palette = "YlGnBu")+
  #ylim(0,3)+
  ylab("Promoter conservation (phastcons)")+
  xlab("")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = 'none')

ggsave(filename = "./sra_fomal/results/plot/004_m1_phastcons.pdf",device = "pdf",width = 5,height = 4.5)

#TF counts diversity
tf_bd_gene_id<-read.table("./sra_fomal/data/tf_bd_gene_id",col.names = c("gene_id","tf_bd"))
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
expressed_variant<-variant_table[rownames(expressed_detect),]
expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$gene_id=rownames(expressed_variant)
gene_tf_count<-tf_bd_gene_id%>%group_by(gene_id)%>%summarise(tf_n=length(tf_bd))
gene_tf_type<-tf_bd_gene_id%>%group_by(gene_id)%>%summarise(type=length(unique(tf_bd)))%>%na.omit()
expressed_variant<-merge(expressed_variant,gene_tf_count,by="gene_id")
expressed_variant<-merge(expressed_variant,gene_tf_type,by="gene_id")
expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$group<-factor(c(rep(c(6,5,4,3,2,1),each=2726),rep(1,2)),levels = c(1,2,3,4,5,6))


# cor<-cor.test(expressed_variant$vst.variance.standardized,expressed_variant$tf_n,method = "spearman")
# ggplot(expressed_variant,mapping = aes(x = log2(tf_n), y = log2(vst.variance.standardized)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
#   theme_bw()+
#   ylab("Standardized expression variation (log2)")+
#   xlab("TF diversity (log2)")+
#   theme(axis.title.x=element_text(size=18), 
#         axis.title.y=element_text(size=16),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16))+
#   annotate("text", x=7, y=6.5, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)+
#   geom_smooth(method = "lm",colour="black")

df_p_val <- expressed_variant%>%
  wilcox_test(formula = tf_n ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position()%>%
  slice(1,6,10,13,15)%>%
  mutate(y.position=c(550,545,540,400,380))
ggplot(expressed_variant,aes(x=group,y=tf_n))+geom_boxplot(aes(fill=group),outlier.alpha = 0, notch=TRUE)+
  stat_pvalue_manual(data = df_p_val,label = '{p.adj}',label.size=4,tip.length = 0,bracket.shorten = 0.1)+
  scale_fill_brewer(palette = "YlGnBu")+
  #ylim(0,3)+
  ylab("TF counts")+
  xlab("")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = 'none')

ggsave(filename = "./sra_fomal/results/plot/004_m1_tfbd_variance.pdf",device = "pdf",width = 5,height = 4.5)

#RNA tf bd dyversity
tf_bd_gene_id<-read.table("./sra_fomal/results/result_table/RNAbd_geneid",col.names = c("gene_id","tf_bd"))
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
expressed_variant<-variant_table[rownames(expressed_detect),]
expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$gene_id=rownames(expressed_variant)
gene_tf_count<-tf_bd_gene_id%>%group_by(gene_id)%>%summarise(tf_n=length(tf_bd))
gene_tf_type<-tf_bd_gene_id%>%group_by(gene_id)%>%summarise(type=length(unique(tf_bd)))%>%na.omit()
expressed_variant<-merge(expressed_variant,gene_tf_count,by="gene_id")
expressed_variant<-merge(expressed_variant,gene_tf_type,by="gene_id")
expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$group<-factor(c(rep(c(6,5,4,3,2,1),each=2711),1),levels = c(1,2,3,4,5,6))

df_p_val <- expressed_variant%>%
  wilcox_test(formula = tf_n ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position()%>%
  slice(1,6,10,13,15)%>%
  mutate(y.position=c(260,260,260,255,250))
ggplot(expressed_variant,aes(x=group,y=tf_n))+geom_boxplot(aes(fill=group),outlier.alpha = 0, notch=TRUE)+
  stat_pvalue_manual(data = df_p_val,label = '{p.adj}',label.size=4,tip.length = 0,bracket.shorten = 0.1,step.increase = FALSE)+
  scale_fill_brewer(palette = "YlGnBu")+
  #ylim(0,3)+
  ylab("Promoter conservation (TF counts)")+
  xlab("")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = 'none')

# cor<-cor.test(expressed_variant$vst.variance.standardized,expressed_variant$tf_n,method = "spearman")
# ggplot(expressed_variant,mapping = aes(x = tf_n, y = log2(vst.variance.standardized)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
#   theme_bw()+
#   ylab("Standardized expression variation (log2)")+
#   xlab("RBP diversity")+
#   theme(axis.title.x=element_text(size=18), 
#         axis.title.y=element_text(size=16),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16))+
#   annotate("text", x=180, y=6.5, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)+
#   geom_smooth(method = "lm",colour="black")
ggsave(filename = "./sra_fomal/results/plot/004_std_var_RNA_tfbd.pdf",device = "pdf",width = 5,height = 4.5)

#TF counts
tf_bd_gene_id<-read.table("./sra_fomal/results/result_table/tf_bd_gene_id_forcount",col.names = c("gene_id","tf_bd"))
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
expressed_variant<-variant_table[rownames(expressed_detect),]
expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$gene_id=rownames(expressed_variant)
gene_tf_type<-tf_bd_gene_id%>%group_by(gene_id)%>%summarise(tf_n=length(tf_bd))%>%na.omit()
expressed_variant<-merge(expressed_variant,gene_tf_type,by="gene_id")
expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$group<-factor(c(rep(c(6,5,4,3,2,1),each=2726),rep(1,2)),levels = c(1,2,3,4,5,6))

df_p_val <- expressed_variant%>%
  wilcox_test(formula = tf_n ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position()%>%
  slice(1,6,10,13,15)%>%
  mutate(y.position=c(1140,1100,1000,700,500))
ggplot(expressed_variant,aes(x=group,y=tf_n))+geom_boxplot(aes(fill=group),outlier.alpha = 0, notch=TRUE)+
  stat_pvalue_manual(data = df_p_val,label = '{p.adj}',label.size=4,tip.length = 0,bracket.shorten = 0.1)+
  scale_fill_brewer(palette = "YlGnBu")+
  #ylim(0,3)+
  ylab("TF counts")+
  xlab("")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = 'none')

# cor<-cor.test(expressed_variant$vst.variance.standardized,expressed_variant$tf_n,method = "spearman")
# ggplot(expressed_variant,mapping = aes(x = log2(tf_n), y = log2(vst.variance.standardized)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
#   theme_bw()+
#   ylab("Standardized expression variation (log2)")+
#   xlab("TF counts (log2)")+
#   theme(axis.title.x=element_text(size=18), 
#         axis.title.y=element_text(size=16),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16))+
#   annotate("text", x=7, y=6.5, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)+
#   geom_smooth(method = "lm",colour="black")
ggsave(filename = "./sra_fomal/results/plot/004_m1_tfbd_variance.pdf",device = "pdf",width = 5,height = 4.5)

#RNA TF counts
tf_bd_gene_id<-read.table("./sra_fomal/results/result_table/RNAbd_geneid_for_count",col.names = c("gene_id","tf_bd"))
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
expressed_variant<-variant_table[rownames(expressed_detect),]
expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$gene_id=rownames(expressed_variant)
gene_tf_type<-tf_bd_gene_id%>%group_by(gene_id)%>%summarise(tf_n=length(tf_bd))%>%na.omit()
expressed_variant<-merge(expressed_variant,gene_tf_type,by="gene_id")

expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$group<-factor(c(rep(c(6,5,4,3,2,1),each=2708),rep(1,2)),levels = c(1,2,3,4,5,6))

df_p_val <- expressed_variant%>%
  wilcox_test(formula = tf_n ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position()%>%
  slice(1,6,10,13,15)%>%
  mutate(y.position=c(2000,2000,1800,1600,1200))
ggplot(expressed_variant,aes(x=group,y=tf_n))+geom_boxplot(aes(fill=group),outlier.alpha = 0, notch=TRUE)+
  stat_pvalue_manual(data = df_p_val,label = '{p.adj}',label.size=4,tip.length = 0,bracket.shorten = 0.1)+
  scale_fill_brewer(palette = "YlGnBu")+
  ylim(0,2000)+
  ylab("Promoter conservation (RBP counts)")+
  xlab("")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = 'none')

# cor<-cor.test(expressed_variant$vst.variance.standardized,expressed_variant$tf_n,method = "spearman")
# ggplot(expressed_variant,mapping = aes(x = log2(tf_n), y = log2(vst.variance.standardized)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
#   theme_bw()+
#   ylab("Standardized expression variation (log2)")+
#   xlab("RBP counts (log2)")+
#   theme(axis.title.x=element_text(size=18), 
#         axis.title.y=element_text(size=16),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16))+
#   annotate("text", x=10, y=6.5, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)+
#   geom_smooth(method = "lm",colour="black")

ggsave(filename = "./sra_fomal/results/plot/004_m1_RBP_counts_variance.pdf",device = "pdf",width = 5,height = 4.5)

#gene_age
# variant_table<-read.table("./sra_fomal//results/result_table/seurat_variant_table")
# variant_table$gene_id=rownames(variant_table)
# gene_age<-read.csv("./sra_fomal/data/gene_age.csv")
# gene_age$gene_age<-as.numeric(gsub(">","",gene_age$gene_age))
# variant_table<-merge(variant_table,gene_age,by.x="gene_id",by.y = "ensembl_id")
# expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
# variant_table<-subset(variant_table,variant_table$gene_id%in%rownames(expressed_detect))
# 
# cor<-cor.test(variant_table$vst.variance.standardized,variant_table$gene_age,method = "spearman")
# ggplot(variant_table,mapping = aes(x = log2(gene_age), y = log2(vst.variance.standardized)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
#   theme_bw()+
#   ylab("Standardized expression variation (log2)")+
#   xlab("Gene age (Myr, log2)")+
#   theme(axis.title.x=element_text(size=18),
#         axis.title.y=element_text(size=16),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16))+
#   annotate("text", x=10.5, y=6, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)+
#   geom_smooth(method = "lm",colour="black")+
#   scale_x_continuous(breaks = seq(0,12,2))
# ggsave(filename = "./sra_fomal/results/plot/004_std_var_age.pdf",device = "pdf",width = 5,height = 4.5)

library(tidyverse)
variant_table<-read.table("./sra_fomal//results/result_table/seurat_variant_table")
variant_table$gene_id=rownames(variant_table)
gentree_age<-readxl::read_xlsx("./sra_fomal/data/dm6_ver78_age.xlsx")
variant_table<-merge(variant_table,gentree_age,by.x = "gene_id",by.y = "gene")
variant_table[variant_table$branch==0,"age_class"]="Older"
variant_table[variant_table$branch!=0,"age_class"]="Younger"
variant_table$age_class<-factor(variant_table$age_class,levels = c("Younger","Older"))

df_p_val <- variant_table%>%
  wilcox_test(formula = vst.variance.standardized ~ age_class) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position()%>%
  mutate(y.position=12)

ggplot(variant_table,aes(x=age_class,y=vst.variance.standardized))+
  geom_boxplot(aes(fill=age_class),outlier.alpha = 0,notch=TRUE)+
  stat_pvalue_manual(data = df_p_val,label = '{p}',label.size=4,tip.length = 0,bracket.shorten = 0.1)+
  ylim(0,15)+
  scale_fill_manual(values=c("#8EE5EE","#53868B"))+
  ylab("Standardized expression variation")+
  xlab("")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position = 'none')


ggsave(filename = "./sra_fomal/results/plot/004_m1_genorinage_variance.pdf",device = "pdf",width = 5,height = 4.5)


#gene type
library(ggpubr)
library(ggsci)
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
variant_table$gene_id<-rownames(variant_table)
type_length_table<-read.table("./sra_fomal/data/type_length_table")
variance_type<-merge(variant_table,type_length_table,by="gene_id")
plot_table<-subset(variance_type,variance_type$gene_type=="protein-coding"|variance_type$gene_type=="lncRNA"|variance_type$gene_type=="pseudogene")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
plot_table<-subset(plot_table,plot_table$gene_id%in%rownames(expressed_detect))

my_comparisons = list(c("protein-coding","lncRNA"),c("pseudogene","lncRNA"),c("protein-coding","pseudogene"))
plot_table$gene_type<-factor(plot_table$gene_type,levels = c("protein-coding","lncRNA","pseudogene"))
ggplot(plot_table,aes(gene_type,log2(vst.variance.standardized+0.01),fill=gene_type))+
  geom_boxplot()+
  theme_classic()+ylab("Standardized expression variation (log2)")+
  #guides(fill=guide_legend(title = "gene type"))+
  xlab("")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  scale_fill_aaas()+
  scale_x_discrete(labels=c("Coding","LncRNA","Pseudogene"))+
  stat_compare_means(comparisons = my_comparisons,label.y = c(9,10,11),size=5)+
  theme(legend.position = 'none')
ggsave(filename = "./sra_fomal/results/plot/004_std_var_gene_type.pdf",device = "pdf",width = 5,height = 4.5)

#gene essential
essential_table<-read.csv("./sra/elife-53865-supp2-v2.csv")
plot_table<-merge(variance_type,essential_table,by.x="gene_id",by.y="FBgn")
colnames(plot_table)[11]="phenotype"
my_comparisons = list(c("Lethal","Semi-lethal"),c("Viable","Semi-lethal"),c("Lethal","Viable"))
plot_table$phenotype<-factor(plot_table$phenotype,levels = c("Lethal","Semi-lethal","Viable"))
ggplot(plot_table,aes(phenotype,log2(vst.variance.standardized+0.01),fill=phenotype))+
  geom_boxplot()+
  theme_classic()+
  ylab("Standardized expression variation")+
  #guides(fill=guide_legend(title = "gene type"))+
  xlab("")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=18),
        axis.text.x = element_text(size=16,angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=16))+
  # scale_fill_manual(values = c("#EE1F26","#862461","#61439A","#4F69B5"))+
  # scale_fill_manual(values = c("#C1C2E2","#8CB78D","#525252"))+
  # scale_fill_manual(values =c("#d95f0d","#fc9272","#9ecae1")) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(9,10,11),size=4.0)+
  theme(legend.position = 'none')
ggsave(filename = "./sra_fomal/results/plot/005_coding_var_essential.pdf",device = "pdf",width = 4,height = 5)

gene_age_essential<-merge(variant_table,essential_table,by.x="gene_id",by.y="FBgn")
colnames(gene_age_essential)[11]="phenotype"
ggplot(gene_age_essential,aes(phenotype,gene_age,fill=phenotype))+geom_boxplot()+theme_classic()+ylab("Gene age")+
  #guides(fill=guide_legend(title = "gene type"))+
  xlab("")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=18),
        axis.text.x = element_text(size=16,angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=16))+
  # scale_fill_manual(values = c("#EE1F26","#862461","#61439A","#4F69B5"))+
  # scale_fill_manual(values = c("#C1C2E2","#8CB78D","#525252"))+
  # scale_fill_manual(values =c("#d95f0d","#fc9272","#9ecae1")) +
  stat_compare_means(comparisons = my_comparisons,label.y = c(2500,3000,3500),size=4.0)+
  theme(legend.position = 'none')

#variation&development
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
fpkm<-fread("~/sra_fomal/data/fpkm") %>% as.data.frame()
rownames(fpkm)=fpkm$V1
fpkm<-fpkm[,-1]
metadata<-read.csv("./sra_fomal/data/metadata.csv",row.names = 1)
metadata<-metadata[colnames(fpkm),]
###mean
mean_vairance_by_sample<-data.frame(sample_id=metadata$curr_SRX,simpled_dev=metadata$simpled_dev,
                               broad_stage=metadata$broad_stage,
                               mean_variace=apply(fpkm, 2, function(x){mean(na.omit(variant_table[rownames(fpkm)[which(x>1)],"vst.variance.standardized"]))}))
plot_table <- mean_vairance_by_sample%>%na.omit()%>%group_by(simpled_dev)%>%summarise(mean = mean(mean_variace))%>%as.data.frame()
#                                                        quantile50=quantile(mean_age)[3],
#                                                         quantile75=quantile(mean_age)[4])
plot_table<-merge(plot_table,unique(metadata[,c(5,6)]),by="simpled_dev")
simpled_dev_levels<-c("unfertilized egg","activated egg","e1","e2","e3","e4","e5","e6","e7","e8","e11","e12","e13","e14","e15","e16","e17","l1","l2","l3","pu6","pu7","pu8","pu9","pu10","a0","a1","a2","a3","a4","a5","a6","a7","a8","a10","a11","a14","a20","a24","a25","a28","a42")
plot_table<-subset(plot_table,plot_table$simpled_dev%in%simpled_dev_levels)
rownames(plot_table)=plot_table$simpled_dev
plot_table<-plot_table[simpled_dev_levels,]
plot_table$simpled_dev<-factor(plot_table$simpled_dev,levels = simpled_dev_levels)

cor<-cor.test(1:nrow(plot_table),plot_table$mean,method = "spearman")
ggplot(plot_table,mapping = aes(x = simpled_dev, y = mean, group=1, color=broad_stage))+geom_line(size=1)+
  theme_classic2()+
  ylab("Mean expressed variance")+
  labs(color = "Stage")+
  scale_color_discrete(breaks = c('egg','embryo','larva','pupae',"adult"))+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=15,angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=15),
        axis.title.y=element_text(size=17),
        legend.text=element_text(size=15),
        legend.title=element_text(size=17))+
  geom_smooth(method = "lm",linetype=3,colour="black")+
  annotate(geom="text", x=35, y=1, size=6,
           label=paste0("rho=",round(cor$estimate,2),",\n P=",format(cor$p.value,2,scientific = TRUE,digits=3)))
###median
median_vairance_by_sample<-data.frame(sample_id=metadata$curr_SRX,simpled_dev=metadata$simpled_dev,
                                    broad_stage=metadata$broad_stage,
                                    median_variace=apply(fpkm, 2, function(x){median(na.omit(variant_table[rownames(fpkm)[which(x>1)],"vst.variance.standardized"]))}))
plot_table <- median_vairance_by_sample%>%na.omit()%>%group_by(simpled_dev)%>%summarise(median = median(median_variace))%>%as.data.frame()
#                                                        quantile50=quantile(mean_age)[3],
#                                                         quantile75=quantile(mean_age)[4])
plot_table<-merge(plot_table,unique(metadata[,c(5,6)]),by="simpled_dev")
simpled_dev_levels<-c("unfertilized egg","activated egg","e1","e2","e3","e4","e5","e6","e7","e8","e11","e12","e13","e14","e15","e16","e17","l1","l2","l3","pu6","pu7","pu8","pu9","pu10","a0","a1","a2","a3","a4","a5","a6","a7","a8","a10","a11","a14","a20","a24","a25","a28","a42")
plot_table<-subset(plot_table,plot_table$simpled_dev%in%simpled_dev_levels)
rownames(plot_table)=plot_table$simpled_dev
plot_table<-plot_table[simpled_dev_levels,]
plot_table$simpled_dev[1]<-"egg0"
plot_table$simpled_dev[2]<-"egg1"
plot_table$simpled_dev<-factor(plot_table$simpled_dev,levels = plot_table$simpled_dev)
cor<-cor.test(1:nrow(plot_table),plot_table$median,method = "spearman")
ggplot(plot_table,mapping = aes(x = simpled_dev, y = median, group=1, color=broad_stage))+geom_line(size=1)+
  theme_classic2()+
  ylab("Standardized expression variation")+
  labs(color = "Stage")+
  scale_color_discrete(breaks = c('egg','embryo','larva','pupa',"adult"))+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=15,angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=15),
        axis.title.y=element_text(size=17),
        legend.text=element_text(size=15),
        legend.title=element_text(size=17))+
  geom_smooth(method = "lm",colour="black")+
  annotate(geom="text", x=8, y=0.45, size=6,
           label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)))
ggsave(filename = "./sra_fomal/results/plot/004_std_var_develop.pdf",device = "pdf",width = 10,height = 6)

#tissue&variation
variant_table<-read.table("./results/result_table/seurat_variant_table")
fpkm<-fread("~/sra_fomal/data/fpkm") %>% as.data.frame()
rownames(fpkm)=fpkm$V1
fpkm<-fpkm[,-1]
metadata<-read.csv("./data/metadata.csv",row.names = 1)
metadata<-metadata[colnames(fpkm),]

#modified metadata
#metadata[grep("carcass",metadata$broad_tissue),"broad_tissue"]=metadata[grep("carcass",metadata$broad_tissue),"simplified_tissue"]
#metadata[grep("multiple tissue",metadata$broad_tissue),"broad_tissue"]=metadata[grep("multiple tissue",metadata$broad_tissue),"simplified_tissue"]
#metadata[grep("adult abdomen minus digestive system and reproductive system",metadata$broad_tissue),"broad_tissue"]="abdomen minus digestive system and reproductive system"
# metadata[is.na(metadata$broad_tissue),"broad_tissue"]=metadata[is.na(metadata$broad_tissue),"cell_line"]
# metadata[is.na(metadata$broad_stage),"broad_stage"]=metadata[is.na(metadata$broad_stage),"cell_line"]

#calculate expressed gene variance
median_vairance_by_sample<-data.frame(sample_id=metadata$curr_SRX,broad_tissue=metadata$broad_tissue,
                                    broad_stage=metadata$broad_stage,
                                    mean_variace=apply(fpkm, 2, function(x){median(na.omit(variant_table[rownames(fpkm)[which(x>1)],"vst.variance.standardized"]))}))
plot_table<-median_vairance_by_sample%>%na.omit()%>%group_by(broad_stage,broad_tissue)%>%summarise(quantile25=quantile(mean_variace)[2],
                                                                                                   quantile50=quantile(mean_variace)[3],
                                                                                                   quantile75=quantile(mean_variace)[4])%>%as.data.frame()
adult_tissue_plot_table<-subset(plot_table,plot_table$broad_stage=="adult")
adult_tissue_plot_table<-adult_tissue_plot_table[order(adult_tissue_plot_table$quantile50),]
adult_tissue_plot_table<-adult_tissue_plot_table[grep("carcass|multiple|abdomen|whole body",adult_tissue_plot_table$broad_tissue,invert = T),]
adult_tissue_plot_table$broad_tissue<-factor(adult_tissue_plot_table$broad_tissue,levels = adult_tissue_plot_table$broad_tissue)

ggplot(adult_tissue_plot_table, aes(x=broad_tissue, y=quantile50)) +
  geom_segment(aes(y = quantile25,
                   x = broad_tissue,
                   yend = quantile75,
                   xend = broad_tissue),
               size=1,color = "gray") +
  geom_point(aes(size=1),color="black",fill="black",stat='identity') +
  xlab("")+
  ylab("Standardized expression variation")+
  theme_bw()+
  theme(axis.text = element_text(size = 16),
        text = element_text(size=18),
        legend.text = element_text(size = 15),
        legend.position = 'none',
        axis.text.x = element_text(size = 16,angle = 60,hjust = 1))
ggsave(filename = "./sra_fomal/results/plot/004_std_var_expressed.pdf",device = "pdf",width = 8,height = 6)

##############################################
############# sex bias analysis ##############
##############################################

#run in linux (dywang@211.69.141.147 /home/dywang/project/sra/deg_results/deseq2)

#read in metadata
metadata<-read.csv("./data/metadata.csv",row.names = 1)
#read in deseq2 results
files<-list.files("./results/deseq2/")
data.list<-list()
for (i in files){
  data.list[[i]]=read.table(paste("./results/deseq2/",i,sep = ""),header = TRUE)
}

# for (i in names(data.list)) {
#   data.list[[i]]$sexbias=ifelse(data.list[[i]]$padj<0.05&abs(data.list[[i]]$log2FoldChange)>2,ifelse(data.list[[i]]$log2FoldChange>2,"male","female"),"unbias")
# }

sex_bias<-lapply(lapply(data.list,function(x){data.list$x$sexbias=ifelse(x$padj<0.05&abs(x$log2FoldChange)>2,ifelse(x$log2FoldChange>2,"male","female"),"unbias")}), table)
sex_bias_summary<-data.frame(tissue=names(sex_bias),female=unlist(lapply(sex_bias,"[[",1)),male=unlist(lapply(sex_bias,"[[",2)))
#adult sexbias gene count in each tissue
adult_summary<-sex_bias_summary[grep("adult",sex_bias_summary$tissue),]
adult_summary$tissue<-gsub("adult_","",adult_summary$tissue)
adult_summary$sum<-apply(adult_summary[c(2,3)], 1, sum)
adult_summary<-adult_summary[order(adult_summary$sum),]
adult_summary<-adult_summary[grep("abdomen|thorax|organism|minus|and|whole",adult_summary$tissue,invert = TRUE),]
colnames(adult_summary)[2:3]=c("Female","Male")
adult_plot_table<-melt(adult_summary,id.vars = c("sum","tissue"))
adult_plot_table$tissue<-gsub("_"," ",adult_plot_table$tissue)
adult_plot_table$tissue<-factor(adult_plot_table$tissue,levels = c("neuron","head","brain","gut","genitalia","malpighian tubule","fat body","gonad"))
ggplot(adult_plot_table,aes(x=tissue,y=value,fill=variable))+geom_bar(stat="identity",position = position_dodge())+
  theme_bw()+
  scale_fill_manual(values = c("#E0848F","#889DB8"))+
  ylab("Sex-bias gene count")+
  labs("")+
  geom_text(aes(label=value), color="black", size=3.5, vjust=-0.5,position = position_dodge(0.9))+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=18),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 16,angle = 60,hjust = 1))

#chr distribution
gene_loc<-read.table("./data/gene_loc",col.names = c("loc","gene_id"))
for (i in names(data.list)){
  data.list[[i]]$sex_bias=ifelse(data.list[[i]]$padj<0.05&abs(data.list[[i]]$log2FoldChange)>2,ifelse(data.list[[i]]$log2FoldChange>2,"male","female"),"unbias")
  data.list[[i]]$gene_loc=gene_loc[match(data.list[[i]]$gene_id,gene_loc$gene_id),"loc"]
  data.list[[i]]$gene_loc[grep("chrX",data.list[[i]]$gene_loc)]="chrX"
  data.list[[i]]$gene_loc[grep("chrY",data.list[[i]]$gene_loc)]="chrY"
  data.list[[i]]<-na.omit(data.list[[i]])
}
sex_bias_summary<-data.frame()
for (i in names(data.list)) {
  table<-data.list[[i]]
  table<-table%>%group_by(gene_loc)%>%summarise(male=count(sex_bias=="male"),female=sum(sex_bias=="female"),unbias=sum(sex_bias=="unbias"))%>%
    filter(.,gene_loc%in%c("chr2L","chr2R","chr3L","chr3R","chrX"))%>%mutate(tissue=i)
  sex_bias_summary<-rbind(sex_bias_summary,table)
}

sex_bias_summary<-sex_bias_summary[grep("adult",sex_bias_summary$tissue),]
sex_bias_summary$tissue<-gsub("adult_","",sex_bias_summary$tissue)
sex_bias_summary<-sex_bias_summary[grep("abdomen|thorax|organism|minus|and|whole",sex_bias_summary$tissue,invert = TRUE),]
sex_bias_summary$sum=apply(sex_bias_summary[,c("male","female","unbias")],1,sum)
fisher_summary<-sex_bias_summary
fisher_sex_bias_chr<-melt(fisher_summary,id.var=(c("gene_loc","tissue","sum","unbias")))
sex_bias_summary$male=sex_bias_summary$male/sex_bias_summary$sum
sex_bias_summary$female=sex_bias_summary$female/sex_bias_summary$sum
adult_sex_bias_chr<-melt(sex_bias_summary,id.var=(c("gene_loc","tissue","sum","unbias")))
adult_sex_bias_chr$variable=factor(adult_sex_bias_chr$variable,levels = c("female","male"))
value<-tapply(fisher_sex_bias_chr$value, paste(fisher_sex_bias_chr$tissue,fisher_sex_bias_chr$variable,fisher_sex_bias_chr$gene_loc=="chrX",sep = "_"), sum)
unbias<-tapply(fisher_sex_bias_chr$unbias, paste(fisher_sex_bias_chr$tissue,fisher_sex_bias_chr$variable,fisher_sex_bias_chr$gene_loc=="chrX",sep = "_"), sum)
fisher_matrix<-data.frame(rbind(matrix(value,nrow = 2),matrix(unbias,nrow = 2)))
colnames(fisher_matrix)=unique(gsub("_TRUE","",gsub("_FALSE","",names(value))))
for (i in 1:ncol(fisher_matrix)) {
 a<-fisher.test(matrix(fisher_matrix[1:4,i],nrow = 2)) 
 fisher_matrix[5,i]=a$p.value
 if(a$p.value<0.001){
   fisher_matrix[6,i]="***"}else if(a$p.value<0.01&a$p.value>0.001){
     fisher_matrix[6,i]="**"}else if(a$p.value<0.05&a$p.value>0.01){
       fisher_matrix[6,i]="*"
     }else{fisher_matrix[6,i]="n.s."}
 fisher_matrix[7,i]=a$estimate
}
sig_table<-data.frame(group=colnames(fisher_matrix),sig=as.character(fisher_matrix[6,]))
adult_sex_bias_chr$group=paste(adult_sex_bias_chr$tissue,adult_sex_bias_chr$variable,sep = "_")
adult_sex_bias_chr<-merge(adult_sex_bias_chr,sig_table,by="group")
adult_sex_bias_chr[grep("chrX",adult_sex_bias_chr$gene_loc,invert = T),"sig"]=""
ggplot(adult_sex_bias_chr,aes(x=gene_loc,y=value,fill=variable))+
  geom_bar(stat = "identity",position = "stack")+
  scale_fill_manual(values = c("#E0848F","#889DB8"))+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=18),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 16,angle = 60,hjust = 1))+
  facet_wrap(.~tissue,ncol = 2,scales='free_y')+
  geom_text(aes(label = sig),
            position = position_stack(vjust = .5))

#whole body

whole_body_summary<-sex_bias_summary[grep("whole_body",sex_bias_summary$tissue),]
colnames(whole_body_summary)[2:3]=c("Female","Male")
whole_body_summary$sum<-apply(whole_body_summary[c(2,3)], 1, sum)
whole_body_summary<-whole_body_summary[order(whole_body_summary$sum),]
whole_body_summary$tissue<-gsub("_whole_body","",whole_body_summary$tissue)
whole_body_plot_table<-melt(whole_body_summary,id.vars = c("sum","tissue"))
whole_body_plot_table$tissue<-factor(whole_body_plot_table$tissue,levels = whole_body_summary$tissue)
ggplot(whole_body_plot_table,aes(x=tissue,y=value,fill=variable))+geom_bar(stat="identity",position = position_dodge())+
  theme_bw()+
  scale_fill_manual(values = c("#E0848F","#889DB8"))+
  ylab("Sex-bias gene count")+
  labs("")+
  geom_text(aes(label=value), color="black", size=3.5, vjust=-0.5,position = position_dodge(0.9))+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=18),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 16,angle = 60,hjust = 1))

whole_body_summary$ratio<-whole_body_summary$female/whole_body_summary$male
whole_body_summary$tissue<-factor(whole_body_summary$tissue,levels = c("embryo","larva","pupa","adult"))
ggplot(whole_body_summary,aes(x=tissue,y=ratio))+geom_bar(stat="identity")+
  theme_bw()+
  ylab("female/male ratio")+
  labs("")+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=18),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 16,angle = 60,hjust = 1))

whole_body_summary<-sex_bias_summary[grep("whole_body",sex_bias_summary$tissue),]
colnames(whole_body_summary)[2:3]=c("Female","Male")
whole_body_summary$tissue<-gsub("_whole_body","",whole_body_summary$tissue)
whole_body_summary$sum=apply(whole_body_summary[,c("Male","Female","unbias")],1,sum)
fisher_summary<-whole_body_summary
fisher_sex_bias_chr<-melt(fisher_summary,id.var=(c("gene_loc","tissue","sum","unbias")))


whole_body_summary$Male=whole_body_summary$Male/whole_body_summary$sum
whole_body_summary$Female=whole_body_summary$Female/whole_body_summary$sum
whole_sex_bias_chr<-melt(whole_body_summary,id.var=(c("gene_loc","tissue","sum","unbias")))
whole_sex_bias_chr$variable=factor(whole_sex_bias_chr$variable,levels = c("Male","Female"))
whole_sex_bias_chr$tissue=factor(whole_sex_bias_chr$tissue,levels = c("embryo","larva","pupa","adult"))

value<-tapply(fisher_sex_bias_chr$value, paste(fisher_sex_bias_chr$tissue,fisher_sex_bias_chr$variable,fisher_sex_bias_chr$gene_loc=="chrX",sep = "_"), sum)
unbias<-tapply(fisher_sex_bias_chr$unbias, paste(fisher_sex_bias_chr$tissue,fisher_sex_bias_chr$variable,fisher_sex_bias_chr$gene_loc=="chrX",sep = "_"), sum)
fisher_matrix<-data.frame(rbind(matrix(value,nrow = 2),matrix(unbias,nrow = 2)))
colnames(fisher_matrix)=unique(gsub("_TRUE","",gsub("_FALSE","",names(value))))
for (i in 1:ncol(fisher_matrix)) {
  a<-fisher.test(matrix(fisher_matrix[1:4,i],nrow = 2)) 
  fisher_matrix[5,i]=a$p.value
  if(a$p.value<0.001){
    fisher_matrix[6,i]="***"}else if(a$p.value<0.01&a$p.value>0.001){
      fisher_matrix[6,i]="**"}else if(a$p.value<0.05&a$p.value>0.01){
        fisher_matrix[6,i]="*"
      }else{fisher_matrix[6,i]="n.s."}
  fisher_matrix[7,i]=a$estimate
}
sig_table<-data.frame(group=colnames(fisher_matrix),sig=as.character(fisher_matrix[6,]))
colnames(whole_sex_bias_chr)[5]="sex"
whole_sex_bias_chr$group=paste(whole_sex_bias_chr$tissue,whole_sex_bias_chr$sex,sep = "_")
whole_sex_bias_chr<-merge(whole_sex_bias_chr,sig_table,by="group")
whole_sex_bias_chr[grep("chrX",whole_sex_bias_chr$gene_loc,invert = T),"sig"]=""
ggplot(whole_sex_bias_chr,aes(x=gene_loc,y=value,fill=sex))+
  geom_bar(stat = "identity",position = "stack")+
  scale_fill_manual(values = c("#E0848F","#889DB8"))+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=18),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 16,angle = 60,hjust = 1),
        legend.position = "top")+
  facet_grid(.~tissue,scales='free_y')+
  geom_text(aes(label = sig),
            position = position_stack(vjust = .5))

####bias gene exist tissues
for (i in names(data.list)) {
  data.list[[i]]$sexbias=ifelse(data.list[[i]]$padj<0.05&abs(data.list[[i]]$log2FoldChange)>2,ifelse(data.list[[i]]$log2FoldChange>2,"male","female"),"unbias")
  data.list[[i]]$tissue=i
}
all_sex_bias<-na.omit(do.call(rbind,data.list))
adult_sex_bias<-all_sex_bias[grep("adult",all_sex_bias$tissue),]
adult_sex_bias<-subset(adult_sex_bias,adult_sex_bias$sexbias=="male"|adult_sex_bias$sexbias=="female")
rownames(adult_sex_bias)=1:nrow(adult_sex_bias)
adult_sex_bias$tissue<-gsub("adult_","",adult_sex_bias$tissue)
adult_sex_bias<-adult_sex_bias[grep("abdomen|thorax|organism|minus|and|whole",adult_sex_bias$tissue,invert = TRUE),]
switched_genes<-adult_sex_bias%>%group_by(gene_id)%>%summarise(class=length(table(sexbias)))%>%filter(.,class==2)
pure_genes<-adult_sex_bias%>%group_by(gene_id)%>%summarise(class=length(table(sexbias)))%>%filter(.,class==1)
switched_genes<-subset(adult_sex_bias,adult_sex_bias$gene_id%in%switched_genes$gene_id)
pure_genes<-subset(adult_sex_bias,adult_sex_bias$gene_id%in%pure_genes$gene_id)
pure_genes<-pure_genes%>%group_by(gene_id)%>%
  summarise(male=count(sexbias=="male"),female=count(sexbias=="female"))%>%
  mutate(n=male+female)%>%group_by(n)%>%
  summarise(male=count(male!=0),female=count(female!=0))%>%
  mutate(sum=male+female)
colnames(pure_genes)[c(2,3)]=c("Male","Female")
sex_bias_conservasion<-melt(pure_genes[,1:3],id.vars = c("n"))
sex_bias_conservasion<-subset(sex_bias_conservasion,sex_bias_conservasion$value!=0)
sex_bias_conservasion$n=factor(sex_bias_conservasion$n,levels = 8:1)
ggplot(sex_bias_conservasion,aes(x=n,y=value,fill=variable))+geom_bar(stat = "identity",width = 0.9)+coord_flip()+
  theme_bw()+
  scale_fill_manual(values = c("#889DB8","#E0848F"))+
  xlab("Shared tissues")+
  labs(title = "Number of Sex-Biased Genes")+
  geom_text(aes(label=value),hjust=-0.1,color="black",size=5)+
  ylim(0,7000)+
  scale_y_continuous(expand = c(0.01,0),limits = c(0,7500))+
  theme(axis.text.y = element_text(size=18),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=18),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 0.2,  # 顶部边缘距离
                             r = 0.2,  # 右边边缘距离
                             b = 0.2,  # 底部边缘距离
                             l = 0.2,  # 左边边缘距离
                             unit = "cm"))+
  
  facet_grid(.~variable)

###whole_body
whole_sex_bias<-all_sex_bias[grep("whole_body",all_sex_bias$tissue),]
whole_sex_bias<-subset(whole_sex_bias,whole_sex_bias$sexbias=="male"|whole_sex_bias$sexbias=="female")
rownames(whole_sex_bias)=1:nrow(whole_sex_bias)
whole_sex_bias$tissue<-gsub("_whole_body","",whole_sex_bias$tissue)
switched_genes<-whole_sex_bias%>%group_by(gene_id)%>%summarise(class=length(table(sexbias)))%>%filter(.,class==2)
pure_genes<-whole_sex_bias%>%group_by(gene_id)%>%summarise(class=length(table(sexbias)))%>%filter(.,class==1)
switched_genes<-subset(whole_sex_bias,whole_sex_bias$gene_id%in%switched_genes$gene_id)
pure_genes<-subset(whole_sex_bias,whole_sex_bias$gene_id%in%pure_genes$gene_id)
pure_genes<-pure_genes%>%group_by(gene_id)%>%
  summarise(male=count(sexbias=="male"),female=count(sexbias=="female"))%>%
  mutate(n=male+female)%>%group_by(n)%>%
  summarise(male=count(male!=0),female=count(female!=0))%>%
  mutate(sum=male+female)
colnames(pure_genes)[c(2,3)]=c("Male","Female")
sex_bias_conservasion<-melt(pure_genes[,1:3],id.vars = c("n"))
sex_bias_conservasion<-subset(sex_bias_conservasion,sex_bias_conservasion$value!=0)
sex_bias_conservasion$n=factor(sex_bias_conservasion$n,levels = 8:1)
ggplot(sex_bias_conservasion,aes(x=n,y=value,fill=variable))+geom_bar(stat = "identity",width = 0.9)+coord_flip()+
  theme_bw()+
  scale_fill_manual(values = c("#889DB8","#E0848F"))+
  xlab("Shared stages")+
  labs(title = "Number of Sex-Biased Genes")+
  geom_text(aes(label=value),hjust=-0.1,color="black",size=5)+
  ylim(0,7000)+
  scale_y_continuous(expand = c(0.01,0),limits = c(0,7500))+
  theme(axis.text.y = element_text(size=18),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=18),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 0.2,  # 顶部边缘距离
                             r = 0.2,  # 右边边缘距离
                             b = 0.2,  # 底部边缘距离
                             l = 0.2,  # 左边边缘距离
                             unit = "cm"))+
  facet_grid(.~variable)


melt_table<-reshape2::melt(pure_genes,id.vars=c("n","sum"))
melt_table$n<-factor(melt_table$n,levels = c("switched",8:1))
melt_table$variable<-factor(melt_table$variable,levels = c("female","male"))
ggplot(melt_table,aes(n,value,fill=variable))+
  geom_bar(stat = "identity",position = "fill",width=0.9)+
  geom_hline(aes(yintercept= 0.5),linetype="dashed")+
  xlab("#Tissue")+
  ylab("")+
  scale_fill_manual(values = c("#889DB8","#E0848F"),breaks=c("male","female"))+
  coord_flip()+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        text = element_text(size = 18),
        legend.position = 'top',
        plot.margin = margin(t = 0.3,  # 顶部边缘距离
                             r = 0.3,  # 右边边缘距离
                             b = 0.3,  # 底部边缘距离
                             l = 0.3,  # 左边边缘距离
                             unit = "cm"))
###X-link&sexbias

gene_loc<-read.table("./data/gene_loc",col.names = c("loc","gene_id"))
pure_genes<-adult_sex_bias%>%group_by(gene_id)%>%summarise(class=length(table(sexbias)))%>%filter(.,class==1)
pure_genes<-subset(adult_sex_bias,adult_sex_bias$gene_id%in%pure_genes$gene_id)
gene_loc[grep("chrX",gene_loc$loc),"loc"]="chrX"
gene_loc[grep("chrY",gene_loc$loc),"loc"]="chrY"
tissue_freq<-as.data.frame(table(pure_loc$gene_id))
colnames(tissue_freq)=c("gene_id","n")
n_loc<-merge(tissue_freq,gene_loc,by="gene_id",all.x=TRUE)
n_loc[grep("chrX",n_loc$loc,invert = T),"loc"]="Autosomal"
n_loc[grep("chrX",n_loc$loc),"loc"]="X-link"
n_loc_plot<-as.data.frame(table(n_loc$n,n_loc$loc))
ggplot(n_loc_plot,aes(Var1,Freq,fill=Var2))+
  geom_bar(stat = "identity",position = "fill",width=0.9)+
  geom_hline(aes(yintercept= 0.5),linetype="dashed")+
  xlab("#Tissue")+
  ylab("")+
  scale_fill_manual(values = c("#686868","#BEBEBE"),breaks=c("X-link","Autosomal"))+
  coord_flip()+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.position = 'top',
        plot.margin = margin(t = 0.3,  # 顶部边缘距离
                             r = 0.3,  # 右边边缘距离
                             b = 0.3,  # 底部边缘距离
                             l = 0.3,  # 左边边缘距离
                             unit = "cm"))
###
evo_info<-read.table("./data/evo_info",header = TRUE)
n_evo<-rbind(tissue_freq,data.frame(gene_id=unique(switched_genes$gene_id),n="switched"))
n_evo<-merge(n_evo,evo_info[1:2],by="gene_id")
n_evo<-merge(n_evo,gene_loc,by="gene_id",all.x=TRUE)
n_evo[grep("chrX",n_evo$loc,invert = T),"loc"]="Autosomal"
n_evo[grep("chrX",n_evo$loc),"loc"]="X-link"
gene_id_sex<-unique(pure_genes[,c(7,8)])
n_evo<-merge(n_evo,gene_id_sex,by="gene_id",all.x=TRUE)
n_evo[is.na(n_evo$sexbias),"sexbias"]<-"switched"
ggplot(n_evo,aes(x=n,y=log2(omega),group=n))+geom_boxplot()+facet_grid(sexbias~loc)

###gene_age
files<-list.files("./results/deseq2/")
data.list<-list()
for (i in files){
  data.list[[i]]=read.table(paste("./results/deseq2/",i,sep = ""),header = TRUE)
}
all_sex_bias<-na.omit(do.call(rbind,data.list))
all_sex_bias<-subset(all_sex_bias,all_sex_bias$baseMean>0.1)
expressed_id<-unique(all_sex_bias$gene_id)
all_sex_bias$sexbias<-ifelse(all_sex_bias$padj<0.05&abs(all_sex_bias$log2FoldChange)>2,ifelse(all_sex_bias$log2FoldChange>2,"male","female"),"unbias")
all_sex_bias$group<-rownames(all_sex_bias)
all_sex_bias<-all_sex_bias[grep("abdomen|thorax|organism|minus|and",all_sex_bias$group,invert = TRUE),]
all_sex_bias$group<-unlist(lapply(strsplit(all_sex_bias$group,"[.]"),"[[",1))
all_sex_bias<-subset(all_sex_bias,all_sex_bias$sexbias!="unbias")
rownames(all_sex_bias)=1:nrow(all_sex_bias)
pure_genes<-subset(all_sex_bias,all_sex_bias$group=="adult_gonad")
#pure_genes<-all_sex_bias%>%group_by(gene_id)%>%summarise(class=length(table(sexbias)))%>%filter(.,class==1)
#pure_genes<-subset(all_sex_bias,all_sex_bias$gene_id%in%pure_genes$gene_id)
gene_age<-read.csv("./data/gene_age.csv",stringsAsFactors = FALSE)
gene_age$gene_age<-as.numeric(gsub(">","",gene_age$gene_age))
gene_age<-gene_age[gene_age$ensembl_id%in%expressed_id,]
gene_age_n<-gene_age%>%group_by(gene_age)%>%summarise(n=n())
gene_age<-merge(gene_age_n,gene_age[1:2],by="gene_age")
pure_genes<-merge(pure_genes,gene_age,by.x="gene_id",by.y="ensembl_id")
pure_genes<-unique(pure_genes[,c(1,8,10,11)])
gene_loc<-read.table("./data/gene_loc",col.names = c("loc","gene_id"))
gene_loc[grep("chrX",gene_loc$loc),"loc"]="X-link"
gene_loc[grep("2L|2R|3L|3R|4",gene_loc$loc),"loc"]="Autosomal"
gene_loc<-gene_loc[grep("chr",gene_loc$loc,invert = T),]
pure_genes<-merge(pure_genes,gene_loc,by="gene_id")
plot_table<-pure_genes%>%group_by(gene_age,loc,sexbias,n)%>%summarise(count=length(sexbias))%>%mutate(fraction=count/n)
ggplot(plot_table,aes(x=gene_age,y=fraction,color=loc))+geom_point()+facet_wrap(sexbias~.)
# gene_type<-read.table("./data/type_length_table")[,1:2]
# a<-merge(pure_genes,gene_type,by="gene_id")
male_x<-subset(plot_table,plot_table$sexbias=="male"&plot_table$loc=="X-link")
male_a<-subset(plot_table,plot_table$sexbias=="male"&plot_table$loc=="Autosomal")
female_x<-subset(plot_table,plot_table$sexbias=="female"&plot_table$loc=="X-link")
female_a<-subset(plot_table,plot_table$sexbias=="female"&plot_table$loc=="Autosomal")
cor.test(male_x$fraction,male_x$gene_age,method = "spearman")
cor.test(male_a$fraction,male_a$gene_age,method = "spearman")
cor.test(female_x$fraction,female_x$gene_age,method = "spearman")
cor.test(female_a$fraction,female_a$gene_age,method = "spearman")
male_plot<-subset(plot_table,plot_table$sexbias=="male")
female_plot<-subset(plot_table,plot_table$sexbias=="female")
ggplot(male_plot,aes(log10(gene_age),fraction))+geom_point(aes(color=loc))+geom_smooth(aes(colour = loc), method="lm", se = F)+theme_bw()+xlab("Age (myr)")+ylab("Proportion of male-biased genes")+
  scale_x_continuous(breaks = c(1,2,3,log10(5000)),labels = c(10,100,1000,5000))+
  theme(axis.title.x= element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.title.y=element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14))+
  annotation_logticks(sides = "b")+scale_colour_discrete("Chr")+
  annotate("text", x=3, 0.68, label="R(X)=-0.82, P=1.994e-06",size=5)+
  annotate("text", x=3, 0.61, label="R(A)=-0.87, P=1.434e-06",size=5)
ggplot(female_plot,aes(log10(gene_age),fraction))+geom_point(aes(color=loc))+geom_smooth(aes(colour = loc), method="lm", se = F)+theme_bw()+xlab("Age (myr)")+ylab("Proportion of female-biased genes")+
  scale_x_continuous(breaks = c(1,2,3,log10(5000)),labels = c(10,100,1000,5000))+
  theme(axis.title.x= element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.title.y=element_text(size=16),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14))+
  annotation_logticks(sides = "b")+scale_colour_discrete("Chr")+
  annotate("text", x=3, 0.24, label="R(X)=0.51, P=0.0085",size=5)+
  annotate("text", x=3, 0.21, label="R(A)=0.61, P=0.0008",size=5)

whole_body_sex_bias<-all_sex_bias[grep("whole_body",all_sex_bias$tissue),]
whole_body_sex_bias<-subset(whole_body_sex_bias,whole_body_sex_bias$sexbias=="male"|whole_body_sex_bias$sexbias=="female")
rownames(whole_body_sex_bias)=1:nrow(whole_body_sex_bias)
whole_body_sex_bias$tissue<-gsub("_whole_body","",whole_body_sex_bias$tissue)
whole_body_sex_bias_freq<-as.data.frame(table(table(whole_body_sex_bias$gene_id)))
whole_body_sex_bias_freq$Var1<-factor(whole_body_sex_bias_freq$Var1,levels = 4:1)
ggplot(whole_body_sex_bias_freq,aes(x=Var1,y=Freq))+geom_bar(stat = "identity",color="#96CDCD",fill="#96CDCD")+coord_flip()+
  theme_void()+
  xlab("Shared tissues")+
  labs("")+
  geom_text(aes(label=Freq),hjust=-0.1,color="black",size=5)+
  ylim(0,3000)+
  scale_y_continuous(expand = c(0.01,0),limits = c(0,7500))+
  theme(axis.text = element_text(size = 16,hjust =-80),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=18),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 0.2,  # 顶部边缘距离
                             r = 0.2,  # 右边边缘距离
                             b = 0.2,  # 底部边缘距离
                             l = 0.2,  # 左边边缘距离
                             unit = "cm"))
#(250,500)


#in linux
#scp debatched_fpkm_with_intergenic dywang@211.69.141.147:/home/dywang/project/sra/data
#scp tissue_stage_batch dywang@211.69.141.147:/home/dywang/project/sra/data
#sed -i "s/ /\t/g" tissue_stage_batch
#echo -e "source activate R4.12\nRscript /home/dywang/project/sra/codes/anova.R -i /home/dywang/project/sra/data/debatched_fpkm_with_intergenic -o /home/dywang/project/sra/debatched_withintergenic_anova_res -l -p 0.0001 -r -m /home/dywang/project/sra/data/tissue_stage_batch --merge_mdata_on curr_SRX -F broad_stage+broad_tissue+sex\nconda deactivate">debatched_withintergenic_anova_r.pbs
#qsub -q batch -V -l nodes=1:ppn=10 debatched_withintergenic_anova_r.pbs





















