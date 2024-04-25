#TF counts
tf_bd_gene_id<-read.table("./sra_fomal/data/tf_bd_gene_id",col.names = c("gene_id","tf_bd"))
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
expressed_variant<-variant_table[rownames(expressed_detect),]
expressed_variant<-expressed_variant[order(expressed_variant$vst.variance.standardized,decreasing = TRUE),]
expressed_variant$group=c(rep(1,2780),rep(2,2780),rep(3,2780),rep(4,2780),rep(5,2780),rep(6,2779))
gene_tfbd_count<-tf_bd_gene_id%>%group_by(gene_id)%>%summarise(type=length(unique(tf_bd)))%>%na.omit()
gene_tfbd_count$group<-expressed_variant[gene_tfbd_count$gene_id,"group"]
gene_tfbd_count<-na.omit(gene_tfbd_count)
gene_tfbd_count$group<-factor(gene_tfbd_count$group)
ggplot(gene_tfbd_count,mapping = aes(x = group, y = type, fill = group))+geom_boxplot()+
  theme_bw()+
  ylab("TF types on promoter")+
  xlab("Group")+
  theme(axis.text.x = element_text(size=15, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=15),
        axis.title.x=element_text(size=17),
        axis.title.y=element_text(size=17),
        legend.text=element_text(size=15),
        legend.title=element_text(size=17),
        title = element_text(size = 17))
