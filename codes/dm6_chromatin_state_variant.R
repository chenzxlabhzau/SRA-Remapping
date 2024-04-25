dm6_chr_state<-read.table("./sra_fomal/results/result_table/dm6.chromatin_9state_S2.gene")
colnames(dm6_chr_state)=c("peak_chr","peak_start","peak_end","state","gene_chr","gene_start","gene_end","gene_id")
variance_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
chr_stat_gene<-subset(dm6_chr_state,dm6_chr_state$gene_id%in%rownames(expressed_detect))
gene_state<-chr_stat_gene%>%group_by(gene_id,state)%>%summarise(n=n())

gene_variance<-data.frame(gene_id=rownames(variance_table),standardized_variance=variance_table$vst.variance.standardized)
gene_variance<-subset(gene_variance,gene_variance$gene_id%in%rownames(expressed_detect))
gene_variance<-gene_variance[order(gene_variance$standardized_variance,decreasing = TRUE),]
gene_variance$group=c(rep(1,2780),rep(2,2780),rep(3,2780),rep(4,2780),rep(5,2780),rep(6,2779))
gene_state<-merge(gene_state,gene_variance,by="gene_id")
gene_state_mean<-gene_state%>%group_by(group,state)%>%summarise(mean=mean(n))
sum<-gene_state_mean%>%group_by(group)%>%summarise(sum=sum(mean))
gene_state_mean<-merge(gene_state_mean,sum,by="group")
gene_state_mean$perc<-gene_state_mean$mean/gene_state_mean$sum
gene_state_mean[grep(1,gene_state_mean$state),"state"]="1_TssA"
gene_state_mean[grep(2,gene_state_mean$state),"state"]="2_TssAFlnk"
gene_state_mean[grep(3,gene_state_mean$state),"state"]="3_TxFlnk"
gene_state_mean[grep(4,gene_state_mean$state),"state"]="4_Tx"
gene_state_mean[grep(5,gene_state_mean$state),"state"]="5_TxWk"
gene_state_mean[grep(6,gene_state_mean$state),"state"]="6_EnhG"
gene_state_mean[grep(7,gene_state_mean$state),"state"]="7_Enh"
gene_state_mean[grep(8,gene_state_mean$state),"state"]="8_ZNF/Rpts"
gene_state_mean[grep(9,gene_state_mean$state),"state"]="9_Het states"
gene_state_mean$group<-factor(gene_state_mean$group,levels = c(1,2,3,4,5,6))
ggplot(gene_state_mean,aes(x=group,y=perc,fill=state))+
  geom_bar(stat = "identity",position = "stack")+
  geom_text(aes(label=percent(perc,accuracy = 0.01)),position = position_stack(vjust = 0.5), size = 3)+
  theme_bw()
