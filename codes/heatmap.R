#debatched fpkm
debatched_fpkm<-fread("~/sra_fomal/data/debatched_fpkm") %>% as.data.frame()
rownames(debatched_fpkm)=debatched_fpkm$V1
debatched_fpkm<-debatched_fpkm[,-1]
sample<-as.data.frame(apply(debatched_fpkm, 2,function(x){sum(x>1)}))
colnames(sample)="expressed_gene"
sample<-subset(sample,sample$expressed_gene>5000)
debatched_fpkm<-debatched_fpkm[,rownames(sample)]
detect<-as.data.frame(apply(debatched_fpkm, 1,function(x){sum(x>1)}))
colnames(detect)="count"
expressed_detect<-subset(detect,detect$count>2)
expressed_detect$detect_ratio=expressed_detect$count/12026
ggplot(expressed_detect, aes(x=detect_ratio)) +
  geom_histogram(binwidth = 0.005,color="indianred",fill="indianred", alpha=0.7, position = 'identity')+
  theme_classic()+
  ylab("Count")+
  xlab("Gene detection rate")+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_blank())

##Correlation 
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
expressed_detect$variance.sd<-variant_table[rownames(expressed_detect),"vst.variance.standardized"]
cor<-cor.test(expressed_detect$detect_ratio,expressed_detect$variance.sd,method = "spearman")
ggplot(expressed_detect,mapping = aes(x = detect_ratio, y = log10(variance.sd+0.1)))+
         geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
  ylim(c(-1,2))+
  theme_bw()+
  ylab("log10 ( Standard variance )")+
  xlab("Gene detection rate")+
  theme(axis.text.x = element_text(size=15, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=15),
        axis.title.x=element_text(size=17),
        axis.title.y=element_text(size=17),
        legend.text=element_text(size=15),
        legend.title=element_text(size=17),
        title = element_text(size = 17))+
  geom_smooth(method = "lm",linetype=1,colour="blue")+
  annotate("text", x=0.8, y=1.8, label="Spearman correlation\n Rho=-0.57, P<2.2e-16",size=5)

##
expressed_detect<-expressed_detect[order(expressed_detect$variance.sd,decreasing = T),]
dynamic_500<-head(rownames(expressed_detect),round(nrow(expressed_detect)*0.05))
constraint_500<-tail(rownames(expressed_detect),round(nrow(expressed_detect)*0.05))
dynamic_500_fpkm<-as.data.frame(t(debatched_fpkm[dynamic_500,]))
constraint_500_fpkm<-as.data.frame(t(debatched_fpkm[constraint_500,]))
metadata<-read.csv("./sra_fomal/data/metadata.csv",row.names = 1)
metadata$group=paste(metadata$broad_tissue,metadata$broad_stage,metadata$sex,sep = "_")
cell_line<-which(is.na(metadata[,"broad_stage"]==TRUE))
metadata<-metadata[-cell_line,]
metadata<-metadata[grep("unknown|carcass|multiple|intersex",metadata$group,invert = TRUE),]



dynamic_500_fpkm<-dynamic_500_fpkm[intersect(rownames(dynamic_500_fpkm),rownames(metadata)),]
metadata<-metadata[intersect(rownames(dynamic_500_fpkm),rownames(metadata)),]
group_median<-apply(dynamic_500_fpkm,2,function(x){tapply(x,metadata$group,median)})
group_median<-group_median[grep("mix",rownames(group_median),invert = T),]
group_median<-as.data.frame(t(group_median[,which(apply(group_median, 2, max)>0)]))
pdf("./sra_fomal/results/plot/003_dy_heatmap.pdf",width = 10,height = 6)
pheatmap::pheatmap(group_median,scale = "row",show_rownames = FALSE,angle_col = 90)
dev.off()

constraint_500_fpkm<-constraint_500_fpkm[intersect(rownames(constraint_500_fpkm),rownames(metadata)),]
metadata<-metadata[intersect(rownames(constraint_500_fpkm),rownames(metadata)),]
group_median<-apply(constraint_500_fpkm,2,function(x){tapply(x,metadata$group,median)})
group_median<-group_median[grep("mix",rownames(group_median),invert = T),]
group_median<-as.data.frame(t(group_median[,which(apply(group_median, 2, max)>0)]))
pheatmap::pheatmap(group_median,scale = "row",show_rownames = FALSE,angle_col = 90)
pdf("./sra_fomal/results/plot/003_cs_heatmap.pdf",width = 10,height = 6)
pheatmap::pheatmap(group_median,scale = "row",show_rownames = FALSE,angle_col = 90)
dev.off()



#write.table(dynamic_500_fpkm,"./results/result_table/dynamic_500_fpkm")
#write.table(constraint_500_fpkm,"./results/result_table/constraint_500_fpkm")















