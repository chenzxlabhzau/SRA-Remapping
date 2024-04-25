library(data.table)
debatched_fpkm<-fread("./sra_fomal/data/debatched_fpkm") %>% as.data.frame()
rownames(debatched_fpkm)=debatched_fpkm[,1]
debatched_fpkm<-debatched_fpkm[,-1]
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
debatched_fpkm<-debatched_fpkm[rownames(expressed_detect),]
metadata<-read.csv("./sra_fomal/data/metadata.csv",row.names = 1)
metadata$group=paste(metadata$broad_stage,metadata$broad_tissue,sep = "_")
metadata<-subset(metadata,metadata$broad_tissue!="tissue carcass"&metadata$broad_tissue!="multiple tissue")
debatched_fpkm<-debatched_fpkm[,metadata$curr_SRX]
median_fpkm=as.data.frame(t(apply(debatched_fpkm,1,function(a){tapply(a,metadata$broad_tissue,median)})))
tau <- function(x){
  ncol=ncol(x)
  x <- x[apply(x, 1, function(y)sum(y >0.1 )>=1),]
  x$max=apply(x[,1:ncol],1,max)
  x$tissue=colnames(x[,1:ncol])[apply(x[,1:ncol],1,which.max)]
  x$tau=apply(1-x[,1:ncol]/x$max,1,sum)/(ncol-1)
  return(as.data.frame(x))
}
median_tau <- tau(median_fpkm)
all_tau<-data.frame(gene_id=rownames(median_tau),tissue=median_tau$tissue,tau=median_tau$tau)

variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
all_tau$standard.variant=variant_table[all_tau$gene_id,"vst.variance.standardized"]

cor<-cor.test(all_tau$tau,all_tau$standard.variant,method = "spearman")
ggplot(all_tau,mapping = aes(x = log2(tau), y = log2(standard.variant)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
  theme_bw()+
  ylab("Standardized expression variation (log2)")+
  xlab("Tau (log2)")+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  geom_smooth(method = "lm",colour="black")+
  annotate("text", x=-1.4, y=5, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)
ggsave(filename = "./sra_fomal/results/plot/004_std_var_log2tau.pdf",device = "pdf",width = 5,height = 4.5)


#######no median

debatched_tau <- tau(debatched_fpkm)
total_tau<-data.frame(gene_id=rownames(debatched_tau),tissue=debatched_tau$tissue,tau=debatched_tau$tau)
total_tau$broad_tissue=metadata[total_tau$tissue,"group"]
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
total_tau$standard.variant=variant_table[total_tau$gene_id,"vst.variance.standardized"]

cor<-cor.test(total_tau$tau,total_tau$standard.variant,method = "spearman")
ggplot(total_tau,mapping = aes(x = tau, y = log2(standard.variant)))+geom_point(shape=21,size=1,color="#2F4F4F",fill="cadetblue1",alpha=0.2)+
  theme_bw()+
  ylab("Standardized expression variation (log2)")+
  xlab("Expression specificity (tau value)")+
  ylim(-6,8)+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  geom_smooth(method = "lm",colour="black")+
  annotate("text", x=0.7, y=5, label=paste0("rho = ",round(cor$estimate,2),",\n P = ",format(cor$p.value,2,scientific = TRUE,digits=3)),size=6)
ggsave(filename = "./sra_fomal/results/plot/004_std_var_log2tau.pdf",device = "pdf",width = 5,height = 4.5)






