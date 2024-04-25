#####read in files
metadata<-read.csv("./sra_fomal/data/metadata.csv")
debatched_fpkm<-fread("./sra_fomal/data/debatched_fpkm") %>% as.data.frame()
rownames(debatched_fpkm)=debatched_fpkm[,1]
debatched_fpkm<-debatched_fpkm[,-1]
table(metadata$broad_stage,metadata$broad_tissue)
human_fly_homolog<-read.csv("./sra_fomal/data/mart_export.txt",sep = ",",header = TRUE)[,c(1,6,5)]
human_fly_one2one<-unique(subset(human_fly_homolog,human_fly_homolog$Drosophila.melanogaster.homology.type=="ortholog_one2one"))

#####human
human_all_fpkm<-read.table("/home/qians/MamDC/Data/Seqdata/WangKuster2019MSBRNAseq/Human.allgene.fpkm.txt")
tissue_info<-data.frame(human_sample=colnames(human_all_fpkm))
tissue_info$human_tissue<-gsub("_rep.","",tissue_info$human_sample)
human_median_fpkm<-t(apply(human_all_fpkm, 1, function(x){tapply(x, tissue_info$human_tissue, median)}))
human_fly_tissue<-read.csv("./sra_fomal/data/human_fly_tissue.csv")
human_co_fpkm<-as.data.frame(human_median_fpkm[,human_fly_tissue$human])

#####fly
metadata$broad_tissue<-gsub("_"," ",metadata$broad_tissue)
fly_all_meta<-subset(metadata,metadata$broad_tissue%in%human_fly_tissue$fly&metadata$broad_stage=="adult")
fly_all_fpkm<-debatched_fpkm[,fly_all_meta$curr_SRX]
fly_median_fpkm<-t(apply(fly_all_fpkm, 1, function(x){tapply(x, fly_all_meta$broad_tissue, median)}))
fly_co_fpkm<-as.data.frame(fly_median_fpkm[,unique(human_fly_tissue$fly)])

#####correlation
human_fly_one2one<-subset(human_fly_one2one,human_fly_one2one$Gene.stable.ID%in%rownames(human_co_fpkm)&human_fly_one2one$Drosophila.melanogaster.gene.stable.ID%in%rownames(fly_co_fpkm))
human_co_fpkm<-human_co_fpkm[human_fly_one2one$Gene.stable.ID,]
fly_co_fpkm<-fly_co_fpkm[human_fly_one2one$Drosophila.melanogaster.gene.stable.ID,]
for (i in 1:nrow(human_fly_tissue)) {
  human_expr<-human_co_fpkm[,human_fly_tissue[i,"human"]]
  fly_expr<-fly_co_fpkm[,human_fly_tissue[i,"fly"]]
  pearson<-cor.test(human_expr,fly_expr,method = "pearson")
  spearman<-cor.test(human_expr,fly_expr,method = "spearman")
  human_fly_tissue[i,"pearson-r"]=pearson$estimate
  human_fly_tissue[i,"pearson-p"]=pearson$p.value
  human_fly_tissue[i,"spearman-r"]=spearman$estimate
  human_fly_tissue[i,"spearman-p"]=spearman$p.value
}

all_cor<-data.frame()
for (i in unique(human_fly_tissue$human)) {
  for(j in unique(human_fly_tissue$fly)){
    all_cor[j,i]=cor.test(human_co_fpkm[,i],fly_co_fpkm[,j],method = "spearman")$estimate
  }
}
pheatmap::pheatmap(all_cor,display_numbers = TRUE)

all_cor_pearson<-data.frame()
for (i in unique(human_fly_tissue$human)) {
  for(j in unique(human_fly_tissue$fly)){
    all_cor_pearson[j,i]=cor.test(human_co_fpkm[,i],fly_co_fpkm[,j],method = "pearson")$estimate
  }
}
pheatmap::pheatmap(all_cor_pearson,display_numbers = TRUE)

#####classified by nessesary
essential_gene_info<-read.csv("./sra_fomal/data/Phenotypic.classes.csv",header = TRUE,row.names = 1)
essential_gene_info<-merge(human_fly_one2one,essential_gene_info,by.x="Gene.stable.ID",by.y="Human_ID",all.x=TRUE)

#0lethal
lethal_gene<-subset(essential_gene_info,essential_gene_info$Pheno_class=="0Lethal")
lethal_human<-human_co_fpkm[lethal_gene$Gene.stable.ID,]
lethal_fly<-fly_co_fpkm[lethal_gene$Drosophila.melanogaster.gene.stable.ID,]
lethal_spearman<-data.frame()
for (i in unique(human_fly_tissue$human)) {
  for(j in unique(human_fly_tissue$fly)){
    lethal_spearman[j,i]=cor.test(lethal_human[,i],lethal_fly[,j],method = "spearman")$estimate
  }
}
pheatmap::pheatmap(lethal_spearman,display_numbers = TRUE,main="0lethal")
#1lethal_disease
semi_gene<-subset(essential_gene_info,essential_gene_info$Pheno_class=="1Lethal_disease")
semi_human<-human_co_fpkm[semi_gene$Gene.stable.ID,]
semi_fly<-fly_co_fpkm[semi_gene$Drosophila.melanogaster.gene.stable.ID,]
semi_spearman<-data.frame()
for (i in unique(human_fly_tissue$human)){
  for(j in unique(human_fly_tissue$fly)){
    semi_spearman[j,i]=cor.test(semi_human[,i],semi_fly[,j],method = "spearman")$estimate
  }
}
pheatmap::pheatmap(semi_spearman,display_numbers = TRUE,main="1Lethal_disease")
#2disease
disease_gene<-subset(essential_gene_info,essential_gene_info$Pheno_class=="2Disease")
disease_human<-human_co_fpkm[disease_gene$Gene.stable.ID,]
disease_fly<-fly_co_fpkm[disease_gene$Drosophila.melanogaster.gene.stable.ID,]
disease_spearman<-data.frame()
for (i in unique(human_fly_tissue$human)){
  for(j in unique(human_fly_tissue$fly)){
    disease_spearman[j,i]=cor.test(disease_human[,i],disease_fly[,j],method = "spearman")$estimate
  }
}
pheatmap::pheatmap(disease_spearman,display_numbers = TRUE,main="2Disease")
#3othergene
other_gene<-subset(essential_gene_info,essential_gene_info$Pheno_class=="3OtherGenes")
other_human<-human_co_fpkm[other_gene$Gene.stable.ID,]
other_fly<-fly_co_fpkm[other_gene$Drosophila.melanogaster.gene.stable.ID,]
other_spearman<-data.frame()
for (i in unique(human_fly_tissue$human)){
  for(j in unique(human_fly_tissue$fly)){
    other_spearman[j,i]=cor.test(other_human[,i],other_fly[,j],method = "spearman")$estimate
  }
}
pheatmap::pheatmap(other_spearman,display_numbers = TRUE,main="3OtherGenes")

####
phenotypic_class<-read.csv("./sra_fomal/data/Phenotypic.classes.csv",row.names = 1,header = TRUE)
for (i in 1:nrow(human_fly_tissue)) {
  table<-data.frame(human=human_co_fpkm[,human_fly_tissue[i,"human"]],fly=fly_co_fpkm[,human_fly_tissue[i,"fly"]],human_id=rownames(human_co_fpkm))
  table<-merge(table,phenotypic_class,by.x="human_id",by.y="Human_ID")
  cor1<-cor.test(table[which(table$Pheno_class=="0Lethal"),"human"],table[which(table$Pheno_class=="0Lethal"),"fly"],method = "pearson")
  cor2<-cor.test(table[which(table$Pheno_class=="1Lethal_disease"),"human"],table[which(table$Pheno_class=="1Lethal_disease"),"fly"],method = "pearson")
  cor3<-cor.test(table[which(table$Pheno_class=="2Disease"),"human"],table[which(table$Pheno_class=="2Disease"),"fly"],method = "pearson")
  cor4<-cor.test(table[which(table$Pheno_class=="3OtherGenes"),"human"],table[which(table$Pheno_class=="3OtherGenes"),"fly"],method = "pearson")
  table[grep("0",table$Pheno_class),"Pheno_class"]="Vital genes"
  table[grep("1",table$Pheno_class),"Pheno_class"]="Vital & disease-suppressing genes"
  table[grep("2",table$Pheno_class),"Pheno_class"]="Disease-suppressing genes"
  table[grep("3",table$Pheno_class),"Pheno_class"]="Other genes"
  table$Pheno_class=factor(table$Pheno_class,levels = c("Vital genes","Vital & disease-suppressing genes","Disease-suppressing genes","Other genes"))
  p1<-ggplot(table,mapping = aes(x = log2(human), y = log2(fly)))+geom_point(aes(color=Pheno_class),size=1,alpha=0.5)+
    theme_bw()+
    ylab(paste("human ",gsub("_"," ",human_fly_tissue[i,"human"]),sep = "- "))+
    xlab(paste("fly ",gsub("_"," ",human_fly_tissue[i,"fly"]),sep = "- "))+
    labs(color="Phenotypic class")+
    theme(axis.text.x = element_text(size=15, hjust = 1, vjust = 1),
          axis.text.y = element_text(size=15),
          axis.title.x=element_text(size=17),
          axis.title.y=element_text(size=17),
          legend.text=element_text(size=15),
          legend.title=element_text(size=17),
          title = element_text(size = 17))+
    geom_smooth(aes(color=Pheno_class),method = "lm")+
    annotate(geom="text", x=-4, y=7.5, size=3,
             label=paste0("rho1=",round(cor1$estimate,2),", P=",format(cor1$p.value,2,scientific = TRUE,digits=3),"\nrho2=",round(cor2$estimate,2),", P=",format(cor2$p.value,2,scientific = TRUE,digits=3),
             "\nrho3=",round(cor3$estimate,2),", P=",format(cor3$p.value,2,scientific = TRUE,digits=3),"\nrho4=",round(cor4$estimate,2),", P=",format(cor4$p.value,2,scientific = TRUE,digits=3)))
    #labs(title = c("Pearson correlation"),subtitle = paste0("rho=",round(cor$estimate,2),", P=",format(cor$p.value,2,scientific = TRUE,digits=3)))
  ggsave(gsub(" ","_",paste("/home/wangdy/sra_fomal/results/plot/expr_cor/expr_cor_",human_fly_tissue[i,"human"],"-",human_fly_tissue[i,"fly"],".pdf",sep = "")),p1,device = "pdf",width = 8,height = 3.5)
}

pearson_cor_p<-data.frame()
pearson_cor_r<-data.frame()
for (i in 1:nrow(human_fly_tissue)) {
  for (j in 1:nrow(human_fly_tissue)) {
    table<-data.frame(human=human_co_fpkm[,human_fly_tissue[i,"human"]],fly=fly_co_fpkm[,human_fly_tissue[j,"fly"]],human_id=rownames(human_co_fpkm))
    table<-merge(table,phenotypic_class,by.x="human_id",by.y="Human_ID")
    table<-subset(table,table$Pheno_class=="0Lethal")
    cor1<-cor.test(table[which(table$Pheno_class=="0Lethal"),"human"],table[which(table$Pheno_class=="0Lethal"),"fly"],method = "pearson")
    pearson_cor_p[human_fly_tissue[i,"human"],human_fly_tissue[j,"fly"]]=cor1$p.value
    pearson_cor_r[human_fly_tissue[i,"human"],human_fly_tissue[j,"fly"]]=cor1$estimate
  }
}



table<-data.frame(human=human_co_fpkm[,human_fly_tissue[3,"human"]],fly=fly_co_fpkm[,human_fly_tissue[8,"fly"]],human_id=rownames(human_co_fpkm))
table<-merge(table,phenotypic_class,by.x="human_id",by.y="Human_ID")
cor1<-cor.test(table[which(table$Pheno_class=="0Lethal"),"human"],table[which(table$Pheno_class=="0Lethal"),"fly"],method = "pearson")
cor2<-cor.test(table[which(table$Pheno_class=="1Lethal_disease"),"human"],table[which(table$Pheno_class=="1Lethal_disease"),"fly"],method = "pearson")
cor3<-cor.test(table[which(table$Pheno_class=="2Disease"),"human"],table[which(table$Pheno_class=="2Disease"),"fly"],method = "pearson")
cor4<-cor.test(table[which(table$Pheno_class=="3OtherGenes"),"human"],table[which(table$Pheno_class=="3OtherGenes"),"fly"],method = "pearson")
p1<-ggplot(table,mapping = aes(x = log2(human), y = log2(fly)))+geom_point(aes(color=Pheno_class),size=1,alpha=0.5)+
  theme_bw()+
  ylab(paste("human ",gsub("_"," ",human_fly_tissue[3,"human"]),sep = "- "))+
  xlab(paste("fly ",gsub("_"," ",human_fly_tissue[8,"fly"]),sep = "- "))+
  theme(axis.text.x = element_text(size=15, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=15),
        axis.title.x=element_text(size=17),
        axis.title.y=element_text(size=17),
        legend.text=element_text(size=15),
        legend.title=element_text(size=17),
        title = element_text(size = 17))+
  geom_smooth(aes(color=Pheno_class),method = "lm")+
  annotate(geom="text", x=-4, y=7.5, size=3,
           label=paste0("rho1=",round(cor1$estimate,2),", P=",format(cor1$p.value,2,scientific = TRUE,digits=3),"\nrho2=",round(cor2$estimate,2),", P=",format(cor2$p.value,2,scientific = TRUE,digits=3),
                        "\nrho3=",round(cor3$estimate,2),", P=",format(cor3$p.value,2,scientific = TRUE,digits=3),"\nrho4=",round(cor4$estimate,2),", P=",format(cor4$p.value,2,scientific = TRUE,digits=3)))
p1

###################################################
cor_max_min<-data.frame()
for (i in 1:nrow(human_fly_tissue)) {
  table<-data.frame(human=human_co_fpkm[,human_fly_tissue[i,"human"]],fly=fly_co_fpkm[,human_fly_tissue[i,"fly"]],human_id=rownames(human_co_fpkm))
  table<-merge(table,phenotypic_class,by.x="human_id",by.y="Human_ID")
  cor1<-cor.test(table[which(table$Pheno_class=="0Lethal"),"human"],table[which(table$Pheno_class=="0Lethal"),"fly"],method = "pearson")
  cor2<-cor.test(table[which(table$Pheno_class=="1Lethal_disease"),"human"],table[which(table$Pheno_class=="1Lethal_disease"),"fly"],method = "pearson")
  cor3<-cor.test(table[which(table$Pheno_class=="2Disease"),"human"],table[which(table$Pheno_class=="2Disease"),"fly"],method = "pearson")
  cor4<-cor.test(table[which(table$Pheno_class=="3OtherGenes"),"human"],table[which(table$Pheno_class=="3OtherGenes"),"fly"],method = "pearson")
  max_cor<-max(cor1$estimate,cor2$estimate,cor3$estimate,cor4$estimate)
  min_cor<-min(cor1$estimate,cor2$estimate,cor3$estimate,cor4$estimate)
  cor_max_min[i,"fly"]=human_fly_tissue[i,"fly"]
  cor_max_min[i,"human"]=human_fly_tissue[i,"human"]
  cor_max_min[i,"max"]=max_cor
  cor_max_min[i,"min"]=min_cor
}


























