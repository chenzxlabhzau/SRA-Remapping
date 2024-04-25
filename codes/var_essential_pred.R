essential<-read.csv("./sra_fomal/data/1-s2.0-S2001037019305628-mmc3.csv",header = TRUE)[-1,c(1,5)]
colnames(essential)=c("gene_id","class")
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
variant_table$gene_id=rownames(variant_table)
variant_table[essential$gene_id,"class"]="Essential"
variant_table[grep("E",variant_table$class,invert = TRUE),"class"]="Non-Essential"
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
variant_table<-subset(variant_table,variant_table$gene_id%in%rownames(expressed_detect))
my_comparisons = list(c("Essential","Non-Essential"))
ggplot(variant_table,aes(x=class,y=log2(vst.variance.standardized),fill=class))+geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,size=5)+theme_bw()+
  ylab("Standardized expression variation (log2)")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))+
  theme(legend.position = 'none')
ggsave(filename = "./sra_fomal/results/plot/var_essential_pred.pdf",device = "pdf",width = 4,height = 5)
