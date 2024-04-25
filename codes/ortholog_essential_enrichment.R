human_fly_homolog<-read.csv("./sra_fomal/data/mart_export.txt",sep = ",",header = TRUE)
essential_gene_info<-read.csv("./sra_fomal/data/Phenotypic.classes.csv",header = TRUE,row.names = 1)
human_fly_homolog[human_fly_homolog==""]<-NA
human_fly_homolog<-human_fly_homolog[which(is.na(human_fly_homolog$Drosophila.melanogaster.gene.stable.ID)==FALSE),]
human_fly_essential<-merge(human_fly_homolog,essential_gene_info,by.x="Gene.stable.ID",by.y="Human_ID",all.y=TRUE)
human_essential<-unique(human_fly_essential[,c(1,5,8)])
human_essential[which(is.na(human_essential$Drosophila.melanogaster.homology.type)),"Drosophila.melanogaster.homology.type"]="no_homolog"

plot_table<-human_essential%>%group_by(Pheno_class,Drosophila.melanogaster.homology.type)%>%
  summarise(n=length(Gene.stable.ID))
sum<-plot_table%>%group_by(Pheno_class)%>%summarise(sum=sum(n))
plot_table<-merge(plot_table,sum,by="Pheno_class")
plot_table$percent=plot_table$n/plot_table$sum
plot_table$diff=plot_table$sum-plot_table$n

for (i in 1:nrow(plot_table)) {
  pheno<-plot_table[i,"Pheno_class"]
  homo_type<-plot_table[i,"Drosophila.melanogaster.homology.type"]
  temp<-plot_table%>%filter(.,Pheno_class!=pheno)%>%filter(.,Drosophila.melanogaster.homology.type==homo_type)
  plot_table[i,"n"]
  test<-fisher.test(matrix(c(plot_table[i,"n"],plot_table[i,"diff"],sum(temp$n),sum(temp$diff)),nrow = 2))
  plot_table[i,"pvalue"]=ifelse(test$p.value<0.05,ifelse(test$p.value<0.01,ifelse(test$p.value<0.001,"***","** "),"*  "),"   ")
  plot_table[i,"or"]=test$estimate
}
plot_table$label<-paste(percent(plot_table$percent,accuracy = 0.01),plot_table$pvalue,sep = " ")
plot_table[grep("0Lethal",plot_table$Pheno_class),"Pheno_class"]="Vital genes"
plot_table[grep("1Lethal_disease",plot_table$Pheno_class),"Pheno_class"]="vital & disease-\nsuppressing genes"
plot_table[grep("2Disease",plot_table$Pheno_class),"Pheno_class"]="Disease-suppressing\n genes"
plot_table[grep("3OtherGenes",plot_table$Pheno_class),"Pheno_class"]="Other genes"
plot_table[grep("no_homolog",plot_table$Drosophila.melanogaster.homology.type),"Drosophila.melanogaster.homology.type"]="No homolog"
plot_table[grep("ortholog_many2many",plot_table$Drosophila.melanogaster.homology.type),"Drosophila.melanogaster.homology.type"]="Many2many"
plot_table[grep("ortholog_one2many",plot_table$Drosophila.melanogaster.homology.type),"Drosophila.melanogaster.homology.type"]="One2many"
plot_table[grep("ortholog_one2one",plot_table$Drosophila.melanogaster.homology.type),"Drosophila.melanogaster.homology.type"]="One2one"
plot_table$Drosophila.melanogaster.homology.type=factor(plot_table$Drosophila.melanogaster.homology.type,levels = c("No homolog","Many2many","One2many","One2one"))
plot_table$Pheno_class<-factor(plot_table$Pheno_class,levels = c("Vital genes","vital & disease-\nsuppressing genes","Disease-suppressing\n genes","Other genes"))
ggplot(plot_table,aes(x=Pheno_class,y=percent,fill=Drosophila.melanogaster.homology.type))+
  geom_bar(stat = "identity")+geom_text(label=plot_table$label,position=position_stack(1), vjust=1.2,size=4)+
  theme_bw()+
  xlab("")+
  ylab("Proportion")+
  labs(fill="Homology type")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16,angle=45,hjust = 1),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))+
  scale_fill_manual(values = c("#A9A9A9","#A9BAE5","#BAE5A9","#e5a9ba"))
ggsave(filename = "./sra_fomal/results/plot/006_homo_type_essential.pdf",device = "pdf",width = 7,height = 5)











