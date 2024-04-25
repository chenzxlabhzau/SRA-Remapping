library(Hmisc)
library(dplyr)
library(data.table)
rm(list = ls())
gc()
exp = fread("./sra_fomal/data/debatched_fpkm") %>% as.data.frame()
row.names(exp) = exp$V1
exp =exp[,-1]
exp = exp[apply(exp, 1, function(x)sum(x>1)>= 2),]
exp = exp[apply(exp, 2, function(x)sum(x>1)>= 5000),]
#exp %<>% dplyr::select(.,!starts_with("Testis")) # exclude testis sample
gene_type<-read.table("./sra_fomal/data/type_length_table")
exp$type <- gene_type[match(row.names(exp),gene_type$gene_id),"gene_type"]
coding = dplyr::filter(exp,grepl("protein-coding",type,ignore.case = T)) %>% dplyr::select(.,-type)
coding1 = coding[apply(coding, 1, function(x)sum(x > 0.1) >=5),]

variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
variant_table<-variant_table%>%filter(.,rownames(variant_table)%in%rownames(exp))
variant_table$type<-gene_type[match(rownames(variant_table),gene_type$gene_id),"gene_type"]
lnc_variant_table<-variant_table%>%filter(.,grepl("lncRNA",type,ignore.case = T))%>%arrange(.,desc(vst.variance.standardized))
dynamic_lnc<-head(lnc_variant_table,round(nrow(lnc_variant_table)*0.05))
constraint_lnc<-tail(lnc_variant_table,round(nrow(lnc_variant_table)*0.05))
dynamic_lnc_fpkm<-exp[rownames(dynamic_lnc),-12050]
constraint_lnc_fpkm<-exp[rownames(constraint_lnc),-12050]

dynamic_lnc_fpkm = dynamic_lnc_fpkm[apply(dynamic_lnc_fpkm, 1, function(x)sum(x > 0.1) >=5),]
constraint_lnc_fpkm = constraint_lnc_fpkm[apply(constraint_lnc_fpkm, 1, function(x)sum(x > 0.1) >=5),]
merge_dat_dy <- cbind(t(coding1), t(dynamic_lnc_fpkm))
merge_dat_cs <- cbind(t(coding1), t(constraint_lnc_fpkm))

#dy
cc_dy<- rcorr(merge_dat_dy, type="pearson")
trans <- function(object){
  pvalue = object[,!(colnames(object) %in% rownames(dynamic_lnc_fpkm))]
  pvalue = pvalue[rownames(object) %in% rownames(dynamic_lnc_fpkm),]
  result = data.frame(lnc = rep(rownames(pvalue), ncol(pvalue)),
                      mRNA = rep(colnames(pvalue), each = nrow(pvalue)),
                      values = matrix(pvalue, ncol = 1)[,1])
  return(result)
}
PVALUE <- trans(cc_dy$P)
Correlation <- trans(cc_dy$r)
PVALUE$ID<- paste(PVALUE$lnc,PVALUE$mRNA,sep = "_")
Correlation$ID<- paste(Correlation$lnc,Correlation$mRNA,sep = "_")
if (all(Correlation$ID == PVALUE$ID)) {
  Pair<- cbind(Correlation,PVALUE)
}
Pair = Pair[,c(1,2,3,7)]
colnames(Pair) = c("lnc","coding","R","Pvalue")
fwrite(Pair,file = "./sra_fomal/results/result_table/dylnc_coding_pair",nThread = 5)

#cs
cc_cs<- rcorr(merge_dat_cs, type="pearson")
trans <- function(object){
  pvalue = object[,!(colnames(object) %in% rownames(constraint_lnc_fpkm))]
  pvalue = pvalue[rownames(object) %in% rownames(constraint_lnc_fpkm),]
  result = data.frame(lnc = rep(rownames(pvalue), ncol(pvalue)),
                      mRNA = rep(colnames(pvalue), each = nrow(pvalue)),
                      values = matrix(pvalue, ncol = 1)[,1])
  return(result)
}
PVALUE <- trans(cc_cs$P)
Correlation <- trans(cc_cs$r)
PVALUE$ID<- paste(PVALUE$lnc,PVALUE$mRNA,sep = "_")
Correlation$ID<- paste(Correlation$lnc,Correlation$mRNA,sep = "_")
if (all(Correlation$ID == PVALUE$ID)) {
  Pair<- cbind(Correlation,PVALUE)
}
Pair = Pair[,c(1,2,3,7)]
colnames(Pair) = c("lnc","coding","R","Pvalue")
fwrite(Pair,file = "./sra_fomal/results/result_table/dylnc_coding_pair",nThread = 5)

dylnc_coding_pair = fread("./sra_fomal/results/result_table/dylnc_coding_pair") %>% as.data.frame()
cslnc_coding_pair = fread("./sra_fomal/results/result_table/cs_lnc_coding_pair") %>% as.data.frame()
dylnc_coding_pair<-dylnc_coding_pair%>%filter(.,Pvalue<0.01)%>%filter(.,R>0.65)
cslnc_coding_pair<-cslnc_coding_pair%>%filter(.,Pvalue<0.01)%>%filter(.,R>0.65)

library(clusterProfiler)
library(org.Dm.eg.db)
gene.df <- bitr(unique(dylnc_coding_pair$coding),fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Dm.eg.db)                   
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
dy_BP<-ego_ALL[rownames(ego_ALL)[ego_ALL$ONTOLOGY=="BP"],]

dy_BP<-dy_BP[order(dy_BP$p.adjust)[1:10],]
dy_BP<-dy_BP[order(dy_BP$Count),]
dy_BP$Description<-factor(dy_BP$Description,levels = dy_BP$Description)
ggplot(data = dy_BP, # 绘图使用的数据
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
ggsave(filename = "./sra_fomal/results/plot/005_dylnc_go.pdf",device = "pdf",width = 8,height = 4.5)
##############

gene.df <- bitr(unique(cslnc_coding_pair$coding),fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Dm.eg.db)                   
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
cs_BP<-ego_ALL[rownames(ego_ALL)[ego_ALL$ONTOLOGY=="BP"],]

cs_BP<-cs_BP[order(cs_BP$p.adjust)[1:10],]
cs_BP<-cs_BP[order(cs_BP$Count),]
cs_BP$Description<-factor(cs_BP$Description,levels = cs_BP$Description)
ggplot(data = cs_BP, # 绘图使用的数据
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
ggsave(filename = "./sra_fomal/results/plot/005_cslnc_go.pdf",device = "pdf",width = 8,height = 4.5)













