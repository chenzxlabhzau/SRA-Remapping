library(readxl)
library(dplyr)

GSEtime = readr::read_table("/home/shimw/project/MassiveQC/GES60314_time.txt",col_names = c("srr","time"))
GSEtime = left_join(GSEtime,srrs[,c("Run","Bytes")], by = c("srr"="Run"))


ll = names(table(GSEtime$srr)[table(GSEtime$srr)>1])
nn = GSEtime[GSEtime$srr%in%ll,]

GSEtime = GSEtime %>% group_by(srr)  %>% slice( 1 )
GSEtime$`Size (Mb)` = GSEtime$Bytes/(1024*1024)
readr::write_csv(GSEtime,"/home/shimw/project/MassiveQC/GES60314_run_time.csv")
GSEtime = readr::read_csv("/home/shimw/project/MassiveQC/GES60314_run_time.csv")

GSEtime$time
library("scales")
library("lubridate")
library("ggplot2")

GSEtime$time = as.numeric(GSEtime$time)/3600
median(GSEtime$`Size (Mb)`)

density(GSEtime$time)
MaxY1_index <- which.max(density(GSEtime$time)$y)

ggplot(as.data.frame(GSEtime), aes(time)) + 
  geom_density(color= "#d88c9a")+
  geom_vline(xintercept = density(GSEtime$time)$x[MaxY1_index],linetype="dashed")+
  theme_classic()+
  xlab("Run time (minutes)")+ylab("Density")+
  theme(axis.text.x = element_text(size = 13 , color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        # 加粗坐标轴
        panel.grid = element_blank() # 删除网格线
  )

ggsave("/home/shimw/project/MassiveQC/A.pdf",height = 3.7,width = 4.1)


MaxY2_index <- which.max(density(GSEtime$`Size (Mb)`)$y)
ggplot(as.data.frame(GSEtime), aes(`Size (Mb)`)) + 
  geom_density(color= "#99c1b9")+
  geom_vline(xintercept = density(GSEtime$`Size (Mb)`)$x[MaxY2_index],linetype="dashed")+
  theme_classic()+
  ylab("Density")+
  theme(axis.text.x = element_text(size = 13 , color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        # 加粗坐标轴
        panel.grid = element_blank() # 删除网格线
  )+
  scale_x_continuous(breaks=c(0,339,1000,2000,3000))
ggsave("/home/shimw/project/MassiveQC/B.pdf",height = 3.7,width = 4.1)



stat.Accuracy = readr::read_csv("/home/shimw/project/MassiveQC/stat.csv")


stat.Accuracy = stat.Accuracy %>% group_by(Num_contamination)%>% summarise("mean_accuracy" = mean(Accuracy))

ggplot(stat.Accuracy)

ggplot(stat.Accuracy, aes(x=Num_contamination, y=mean_accuracy)) +
  geom_line(color="#db7c26") + geom_point(color="#db7c26")+
  theme_classic()+
  ylab("Accuracy")+xlab("Number of contamination")+
  theme(axis.text.x = element_text(size = 13 , color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        # 加粗坐标轴
        panel.grid = element_blank() # 删除网格线
  )+
  scale_x_continuous(breaks=c(5,10,15,20,25,30))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),limits = c(0,1.1))
ggsave("/home/shimw/project/MassiveQC/C.pdf",height = 3.7,width = 4.1)


sample.Accuracy = readr::read_csv("/home/shimw/project/MassiveQC/sample_stat.csv")
aa=sample.Accuracy[sample.Accuracy$Num_sample==400,]
sample.Accuracy = sample.Accuracy %>% 
  group_by(Num_sample)%>%
  summarise("prob" = mean(Accuracy),sd = sd(Accuracy))


ggplot(sample.Accuracy, aes(x=Num_sample, y=prob)) +
  geom_line(color="#669bbc") + geom_point(color="#669bbc")+
  geom_errorbar(aes(ymin=prob-sd, ymax=prob+sd), width=10,color="#669bbc")+
  theme_classic()+
  ylab("Accuracy")+xlab("Number of samples")+
  theme(axis.text.x = element_text(size = 13 , color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        # 加粗坐标轴
        panel.grid = element_blank() # 删除网格线
  )+
  scale_x_continuous(breaks=c(0,100,200,300,400,500,600,700,800))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),limits = c(0,1.1))
ggsave("/home/shimw/project/MassiveQC/D.pdf",height = 3.7,width = 4.1)

library(dplyr)
library(SEGtool)
data(EbiHumanExpression)

debatched_fpkm = data.table::fread("/home/wangdy/sra_fomal/data/debatched_fpkm",data.table = F)
row.names(debatched_fpkm) = debatched_fpkm$V1
debatched_fpkm = debatched_fpkm[,-1]
meta.info = read.csv("/home/wangdy/sra_fomal/data/metadata.csv",row.names = 1)

meta.info = meta.info[meta.info$broad_stage=="adult"&meta.info$broad_tissue!="whole body",]%>%
  tidyr::drop_na(broad_tissue)%>%
  dplyr::arrange(broad_tissue)
debatched_fpkm=debatched_fpkm[,meta.info$curr_SRX]

k = as.data.frame(t(debatched_fpkm["FBgn0033450",]))
m = cbind(k, meta.info[row.names(k),])%>%tidyr::drop_na(broad_stage)
m%>%dplyr::filter(broad_tissue=="ovary")%>%
  group_by(broad_stage)%>%summarise(mean_fpkm = mean(FBgn0033450))


exp_all = SEGtool::replicates_value_integration(debatched_fpkm,factor(meta.info$broad_tissue))

lllll = purrr::map(levels(factor(meta.info$broad_tissue)),function(x){
  aa = which(meta.info$broad_tissue == x)
  if (length(aa)==1) {
    bb = data.frame(row.names = row.names(debatched_fpkm),
                    "mm"=debatched_fpkm[,aa])
    names(bb) = x
    return(bb)
  }
  bb = apply(debatched_fpkm[,aa], 1, mean)%>%as.data.frame()
  names(bb) = x
  return(bb)
})%>%bind_cols()

data(EbiHumanExpression)
EbiHumanExpression = EbiHumanExpression[1:10000,]
SEGtool_result <- SEGtool(EbiHumanExpression, exp_cutoff = 3,multi_cpu = 1, detect_mod=2, result_outdir='SEGtool_result', draw_heatmap=T, draw_pca=T, draw_plot=FALSE, html_report=FALSE)









