#debatched_fpkm
library(data.table)
debatched_fpkm<-fread("./sra_fomal/data/debatched_fpkm") %>% as.data.frame()
rownames(debatched_fpkm)=debatched_fpkm[,1]
debatched_fpkm<-debatched_fpkm[,-1]

#fpkm
fpkm<-fread("~/sra_fomal/data/fpkm") %>% as.data.frame()
rownames(fpkm)=fpkm$V1
fpkm<-fpkm[,-1]

#metadata
metadata<-read.csv("~/sra_fomal/data/metadata.csv",stringsAsFactors = FALSE,row.names = 1)

#whole_body
whole_body<-subset(metadata,metadata$broad_tissue=="whole body")
whole_body_fpkm<-fpkm[,whole_body$curr_SRX]
whole_body_debatched_fpkm<-debatched_fpkm[,whole_body$curr_SRX]

#UMAP

#fpkm
library(umap)
fpkm_umap<-umap(t(whole_body_fpkm))
fpkm_umap_plot=as.data.frame(fpkm_umap$layout)
colnames(fpkm_umap_plot)=c("UMAP_1","UMAP_2")
fpkm_umap_plot$broad_stage=whole_body$broad_stage
fpkm_umap_plot$broad_stage=factor(fpkm_umap_plot$broad_stage,levels = c("egg","embryo","larva","pupa","adult"),labels=c("Egg","Embyro","Larva","Pupa","Adult"))
ggplot(fpkm_umap_plot,aes(UMAP_1,UMAP_2,color=broad_stage))+geom_point()+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))+
  guides(color=guide_legend(title="Stage"))
ggsave("./sra_fomal/results/plot/suppliment_01_fpkm.pdf",device = "pdf",width=6,height = 4.5)

#debatched_fpkm
debatched_fpkm_umap<-umap(t(whole_body_debatched_fpkm))
debatched_fpkm_umap_plot=as.data.frame(debatched_fpkm_umap$layout)
colnames(debatched_fpkm_umap_plot)=c("UMAP_1","UMAP_2")
debatched_fpkm_umap_plot$broad_stage=whole_body$broad_stage
debatched_fpkm_umap_plot$broad_stage=factor(debatched_fpkm_umap_plot$broad_stage,levels = c("egg","embryo","larva","pupa","adult"),labels=c("Egg","Embyro","Larva","Pupa","Adult"))
ggplot(debatched_fpkm_umap_plot,aes(UMAP_1,UMAP_2,color=broad_stage))+geom_point()+
  theme_classic()+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))+
  guides(color=guide_legend(title="Stage"))
ggsave("./sra_fomal/results/plot/suppliment_01_debatch_fpkm.pdf",device = "pdf",width=6,height = 4.5)

