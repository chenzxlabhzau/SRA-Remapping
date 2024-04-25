files<-list.files("./data/drosophila_sexbias/")
species.list<-list()
adult_whole_body<-read.table("../sra_fomal/results/deseq2/adult_whole_body",header = TRUE)
for (i in files){
  temp=read.csv(paste("./data/drosophila_sexbias/",i,sep = ""),header = TRUE)%>%na.omit()
  temp$genus=unlist(strsplit(i,"[.]"))[1]
  temp$sexbias=ifelse(temp$padj<0.05&abs(temp$log2FoldChange)>2,ifelse(temp$log2FoldChange>2,"male","female"),"unbias")
}
orthlog<-fread("./data/drosophila_sexbias/mart_export.txt")%>%as.data.frame()
orthlog<-orthlog[,c(1,3,6,9,12,15,18)]
orthlog[orthlog==""]<-NA
orthlog<-na.omit(orthlog)
dmel<-as.data.frame(table(orthlog$`Gene stable ID`))


###############################################################################
variant_table<-read.table("./sra_fomal//results/result_table/seurat_variant_table")
variant_table<-variant_table[order(variant_table$vst.variance.standardized,decreasing = TRUE),]
fpkm<-fread("~/sra_fomal/data/fpkm")%>%as.data.frame()
rownames(fpkm)=fpkm$V1
fpkm<-fpkm[,-1]
a<-as.data.frame(apply(fpkm, 1, function(x){sum(x>1)}))
expressed_gene_id<-rownames(a)[a>2]
variant_rank<-subset(variant_table,rownames(variant_table)%in%expressed_gene_id)
varia_gene_id<-rownames(variant_rank)[1:(nrow(variant_rank)*0.05)]
stable_gene_id<-tail(rownames(variant_rank),500)[1:(nrow(variant_rank)*0.05)]

library(clusterProfiler)
library(org.Dm.eg.db)
gene.df <- bitr(stable_gene_id,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Dm.eg.db)                   
gene <- gene.df$ENTREZID
ego_ALL <- enrichGO(gene = gene,
                    OrgDb=org.Dm.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.05,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)
ego_ALL <- as.data.frame(ego_ALL)

BP<-ego_ALL[rownames(ego_ALL)[ego_ALL$ONTOLOGY=="BP"],]
# MF<-ego_ALL[rownames(ego_ALL)[ego_ALL$ONTOLOGY=="MF"][1:20],]
# CC<-ego_ALL[rownames(ego_ALL)[ego_ALL$ONTOLOGY=="CC"][1:20],]
# go_enrich_df <- data.frame(
#   ID=c(BP$ID,CC$ID,MF$ID),                         
#   Description=c(BP$Description,CC$Description,MF$Description),
#   GeneNumber=c(BP$Count,CC$Count,MF$Count),  (count)
#   type=factor(c(rep("biological process",20),  (Ontology)
#                 rep("cellular component",20),
#                 rep("molecular function",20)), 
#               levels=c("biological process", "cellular component","molecular function" )))
# go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
# COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色
# go_enrich_df<-subset(go_enrich_df,go_enrich_df$type=="biological process")
# ggdotchart(go_enrich_df, x = "Description", y = "GeneNumber",
#            color = "indianred",                       
#            add = "segments",                             # 添加棒子
#            ggtheme = theme_pubr(),                        # 改变主题
#            xlab=""
# )
# write.table(BP,"./results/result_table/stable_go_BP")
# library(simplifyEnrichment)
# mat<-GO_similarity(BP$ID)
# df = simplifyGO(mat)

BP<-BP[order(BP$p.adjust)[1:10],]
BP<-BP[order(BP$Count),]
BP$Description<-factor(BP$Description,levels = BP$Description)
ggplot(data = BP, # 绘图使用的数据
       aes(x = Description ,y = Count))+ #横纵坐标及排序
  geom_point(aes(size = Count,color = p.adjust))+ # 气泡大小及颜色设置
  theme_bw()+ # 去除背景色
  coord_flip()+
  scale_colour_gradient(low = "red",high = "purple")+ # 设置气泡渐变色
  labs(x = "", y = "")+ # 设置图例颜色及大小
  theme(axis.title = element_text(size = 13),axis.text = element_text(size = 11),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
        legend.title = element_text(size = 13),legend.text = element_text(size = 11),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
ggsave(filename = "./sra_fomal/results/plot/005_cscoding_go.pdf",device = "pdf",width = 8,height = 4.5)



















ggplot(data=ego_ALL, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.6) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  xlab("") + 
  ylab("Gene Number") + 
  # scale_y_continuous(limits=c(0, 350),expand = c(0,0))+
  # ylim(0,350)+
  labs(title = "The Most Enriched GO Terms In House-keeping Genes")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=20),
        axis.text.x = element_text(size=18,angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=18),
        panel.grid=element_blank(),
        panel.background=element_blank())+
  theme(legend.position = 'none')+
  theme(plot.title = element_text(size = 18))+
  coord_flip()