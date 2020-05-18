library(dplyr)
library(tidyr)
library(umap)
library(ggplot2)
library(Rtsne)
library(randomcoloR)
library(DDRTree)
g_1000<-read.delim("recoded_1000G.raw",header = T,sep = " ")
rownames(g_1000) <- g_1000[,2]
unprocessed_single_cell_data_total_transposed <-g_1000[,7:39188]
population<- read.delim("igsr_samples.tsv",header = T)
population$Sample.name
write.table(population$Sample.name,file="population_EAS.txt",quote = F,sep = "\t",row.names = F)

population<-as.data.frame(population[which(population$Superpopulation.code=="EAS"),])

area_unprocessed_single_cell_data_total_transposed<-intersect(population$Sample.name,rownames(unprocessed_single_cell_data_total_transposed))

g_1000_area<-unprocessed_single_cell_data_total_transposed[area_unprocessed_single_cell_data_total_transposed,]


datadata<-as.data.frame(as.character(colnames(g_1000_area)) )
mutation_genotype<-datadata %>% separate("as.character(colnames(g_1000_area))", into=c("ID", "genotype"), "_")
mutation_genotype <-mutation_genotype[7:length(mutation_genotype[,1]),]



ncol=dim(g_1000_area)[2]
nrow=dim(g_1000_area)[1]
one_data_max<-as.matrix(g_1000_area, ncol=ncol, nrow=nrow)# no need this step
one_data_max_num <- apply(one_data_max,2,as.numeric)
one_data_max_num[is.na(one_data_max_num)] <- 0

tsne_one_data_umap <- umap(one_data_max_num)


tsne_one_data_umap <- Rtsne(one_data_max_num, dims=2, Perplexity=30, max_iter = 5000,theta=0.0,check_duplicates = F)
#tsne_one_data_umap <-prcomp(t(one_data_max_num),rank.=2)
#PC_table <- data.frame(tsne_one_data_umap$rotation)

bind_population_umap<-read.delim("bind_population_umap_DDR_EAS_.txt",header = T)


population_position <- match(rownames(g_1000_area),population$Sample.name)
#bind_population_umap <- cbind(PC_table$PC1,PC_table$PC2,population[population_position,])
bind_population_umap <- cbind(tsne_one_data_umap$Y,population[population_position,])
colnames(bind_population_umap) <- c("X1","X2",c(colnames(population)))


umap_wait<- cbind(as.numeric(bind_population_umap$X1),as.numeric(bind_population_umap$X2)) 

DDRTree_umap_wait<-DDRTree(t(umap_wait))

DDRTree_umap_wait_Z<-DDRTree_umap_wait$Z

bind_population_umap_DDR<-cbind(Z1=DDRTree_umap_wait_Z[1,],Z2=DDRTree_umap_wait_Z[2,],bind_population_umap)


library(randomcoloR)
n <- length(unique(bind_population_umap$Population.code ))
palette <- distinctColorPalette(n)
ggplot(bind_population_umap,mapping = aes(x= X1, y= X2, color= bind_population_umap$Population.code) )+
  geom_point(size = 0.1)+
  scale_color_manual(values =  palette)



bind_population_umap<-read.delim(file="recoded_1000Gtsne.txt",header = T)

means_group <- bind_population_umap %>% 
  group_by(Population.name)%>%
  summarise(mean_t1 = mean(X1),
            mean_t2 = mean(X2))%>% print(n=1)

require("ggrepel")
pdf("EAS_tsne_plot.pdf")
ggplot(bind_population_umap, aes(x=X1, y=X2,colour=bind_population_umap$Population.code))+
  geom_point()+
  scale_color_manual(values =  palette)+
  geom_label_repel(data = means_group, aes(x = mean_t1,y = mean_t2),label =means_group$Population.name, colour=1,size = 5,alpha = .8)+ theme(plot.title=element_text(size = 30,hjust=0.5),legend.title=element_blank())


dev.off()

bind_population_umap_DDR <-bind_population_umap

means_group_Z <- bind_population_umap_DDR %>% 
  group_by(Population.name)%>%
  summarise(mean_t1 = mean(Z1),
            mean_t2 = mean(Z2))%>% print(n=1)
n <- length(unique(bind_population_umap_DDR$Population.code ))
palette <- distinctColorPalette(n)
pdf("EAS_tsne_ddrtree_plot.pdf")
ggplot(bind_population_umap_DDR, aes(x=Z1, y=Z2,colour=bind_population_umap$Population.code))+
  geom_point()+
  scale_color_manual(values =  palette)+
  geom_label_repel(data = means_group_Z, aes(x = mean_t1,y = mean_t2),label =means_group_Z$Population.name, colour=1,size = 5,alpha = .8)+ theme(plot.title=element_text(size = 30,hjust=0.5),legend.title=element_blank())


dev.off()




E11_stable_all<-read.table("e11.alleles_e11_E11_E11_E11_stable_all.txt",head=T)


plotbind_population_umap_DDR<-bind_population_umap_DDR[match(rownames(E11_stable_all),bind_population_umap_DDR$Sample.name),]

plotbind_population_umap_DDR_dim<- cbind(plotbind_population_umap_DDR,E11_stable_all)

plotbind_population_umap_DDR_dim<-as.data.frame(plotbind_population_umap_DDR_dim)



for (area_e11 in colnames(E11_stable_all))
  
{
  
  pdf(paste("EAS ddrtree plot ofplot_e11_each_",area_e11,".pdf",sep=""))
  plot<-ggplot(plotbind_population_umap_DDR_dim, aes(x=Z1, y=Z2,colour=plotbind_population_umap_DDR_dim[,area_e11]))+
    geom_point(size=3)+scale_colour_gradientn(colours = c("#253494","#ffffcc", "#dd1c77"),values = scales::rescale(c(0,0.5)))+labs(title=paste0("E11 ",area_e11))+ theme(plot.title=element_text(size = 30,hjust=0.5),legend.title=element_blank())
  
  
  print(plot)
  dev.off()
}


highlight <- means_group$Population.name =="Esan"

ggplot(bind_population_umap, aes(x=X1, y=X2,colour=bind_population_umap$Population.code))+
  geom_point(size = 1)+
  scale_color_manual(values =  palette)+
  theme(legend.position = "none") +
  geom_label_repel(data = means_group, aes(x = mean_t1,y = mean_t2,fill=highlight,label =means_group$Population.name), colour=1,size = 5,alpha = .6)+
  scale_fill_manual(values = c("Grey","Red"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

)



write.table(bind_population_umap,file="recoded_1000Gtsne.txt",quote = F,sep = "\t")
