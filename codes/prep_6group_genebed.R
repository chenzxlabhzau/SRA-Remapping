tss1000<-read.table("./sra_fomal/data/tss1000.bed")
tss1000[tss1000<0]<-0
variant_table<-read.table("./sra_fomal/results/result_table/seurat_variant_table")
expressed_detect<-read.table("./sra_fomal/results/result_table/expressed_detect")
tss1000$standard_variance<-variant_table[tss1000$V4,"vst.variance.standardized"]
tss1000<-tss1000[order(tss1000$standard_variance,decreasing = TRUE),]
tss1000<-unique(subset(tss1000,tss1000$V4%in%rownames(expressed_detect)))
tss1000$group=c(rep(1,2780),rep(2,2780),rep(3,2780),rep(4,2780),rep(5,2780),rep(6,2778))
group1<-subset(tss1000,tss1000$group==1)[,1:4]
group2<-subset(tss1000,tss1000$group==2)[,1:4]
group3<-subset(tss1000,tss1000$group==3)[,1:4]
group4<-subset(tss1000,tss1000$group==4)[,1:4]
group5<-subset(tss1000,tss1000$group==5)[,1:4]
group6<-subset(tss1000,tss1000$group==6)[,1:4]
write.table(group1,"./sra_fomal/results/result_table/group1.bed",sep = "\t",quote = FALSE, col.names = FALSE,row.names = FALSE)
write.table(group2,"./sra_fomal/results/result_table/group2.bed",sep = "\t",quote = FALSE, col.names = FALSE,row.names = FALSE)
write.table(group3,"./sra_fomal/results/result_table/group3.bed",sep = "\t",quote = FALSE, col.names = FALSE,row.names = FALSE)
write.table(group4,"./sra_fomal/results/result_table/group4.bed",sep = "\t",quote = FALSE, col.names = FALSE,row.names = FALSE)
write.table(group5,"./sra_fomal/results/result_table/group5.bed",sep = "\t",quote = FALSE, col.names = FALSE,row.names = FALSE)
write.table(group6,"./sra_fomal/results/result_table/group6.bed",sep = "\t",quote = FALSE, col.names = FALSE,row.names = FALSE)