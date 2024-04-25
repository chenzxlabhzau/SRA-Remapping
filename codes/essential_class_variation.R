#####read in files
human_fly_homolog<-read.csv("./sra_fomal/data/mart_export.txt",sep = ",",header = TRUE)[,c(1,6,5)]
human_fly_one2one<-unique(subset(human_fly_homolog,human_fly_homolog$Drosophila.melanogaster.homology.type=="ortholog_one2one"))
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
variant_table<-variant_table[intersect(rownames(expressed_detect),rownames(expressed_detect)),]
variant_table<-variant_table[intersect(rownames(variant_table),human_fly_one2one$Drosophila.melanogaster.gene.stable.ID),]
variance_type<-data.frame(fly_id=rownames(variant_table),standard_variance=variant_table$vst.variance.standardized)

pheno_class<-read.csv("./sra_fomal/data/Phenotypic.classes.csv",header = TRUE,row.names = 2)[,-1]
variance_type<-merge(variance_type,human_fly_one2one,by.x="fly_id",by.y="Drosophila.melanogaster.gene.stable.ID")
variance_type$pheno_type<-pheno_class[variance_type$Gene.stable.ID,"Pheno_class"]
variance_type<-na.omit(variance_type)
ggplot(variance_type,aes(x=pheno_type,y=log2(standard_variance+0.001),fill=pheno_type))+geom_violin()+geom_boxplot()
