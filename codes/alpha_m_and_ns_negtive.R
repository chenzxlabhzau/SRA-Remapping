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


dmel_alpha_m_and_NS_negative<-read.csv("./sra_fomal/data/dmel_alpha_m_and_NS_negative.csv",row.names = 1)
seurat_variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
dmel_alpha_m_and_NS_negative$standardized_variance<-seurat_variant_table[rownames(dmel_alpha_m_and_NS_negative),"vst.variance.standardized"]
dmel_alpha_m_and_NS_negative<-na.omit(dmel_alpha_m_and_NS_negative)
dmel_alpha_m_and_NS_negative<-dmel_alpha_m_and_NS_negative[intersect(rownames(dmel_alpha_m_and_NS_negative),rownames(expressed_detect)),]
