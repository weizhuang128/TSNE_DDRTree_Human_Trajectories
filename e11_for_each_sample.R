
####################### do e11 for each sample

load('/Users/weiz/Desktop/R_code/e11/e11.alleles.RData')
load('/Users/weiz/Desktop/R_code/e11/e11.11.F.RData')

Pseudotime_data_table_EAS <- read.table("/Users/weiz/Desktop/R_code/e11/Pseudotime_data_table_EAS.txt",header = T)


setwd("/Users/weiz/Desktop/R_code/e11/23andme/")

human_name<-data.frame(file_name=list.files("/Users/weiz/Desktop/R_code/e11/23andme/"))

human_name_inter<-intersect(Pseudotime_data_table_EAS$sample_name,substring(human_name$file_name, 1, 7))
rownames(Pseudotime_data_table_EAS) <-Pseudotime_data_table_EAS$sample_name
Pseudotime_data_table<-Pseudotime_data_table_EAS[human_name_inter,]

######
E11_states_table <- data.frame(matrix(NA, ncol = 11, nrow = 0)) 
colnames(E11_states_table) <-c("African","European","India","Malay","SouthChineseDai","SouthwestChineseYi","EastChinese","Japanese","NorthChineseOroqen","Yakut","American")

for (states_num in unique(Pseudotime_data_table$State))
  
{
  
  E11_table <-  data.frame(matrix(NA, ncol = 11, nrow = 0)) 
  colnames(E11_table) <-c("African","European","India","Malay","SouthChineseDai","SouthwestChineseYi","EastChinese","Japanese","NorthChineseOroqen","Yakut","American")
  
  #states_num <- "5"
  grouped_person<-Pseudotime_data_table[Pseudotime_data_table$State==states_num,1]
  print(states_num)
  for(person_name in grouped_person)
  {
    print(person_name)
    #person_name<-"HG02046"
    genotype <- read.delim(file = paste(person_name,"_my_23andme.1.txt",sep = ""),col.names = c("rsid", "chr", "site","genotype"),row.names = NULL,header = F)
    
    # Use E11
    res <- tfrdpub(genotype, 11, e11.alleles, e11.11.F)
    # Use E11
    ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "E11")
    E11_table[person_name,]<-ances$q
  }
  write.table(E11_table,paste("E11_table_states",states_num,sep = "",".txt"),quote = F,sep = "\t")
  E11_table_sum <- apply(E11_table,2,sum)
  E11_states_table[states_num,]<-E11_table_sum
  
}


######


E11_states_table_reshaped<- melt(t(na.omit(E11_states_table)))

colnames(E11_states_table_reshaped) <- c("Group","State","Percentage")

n <- length(unique(Pseudotime_data_table$State))
palette <- distinctColorPalette(11)
pdf("E11_states_table_reshaped.pdf")
ggplot(E11_states_table_reshaped, aes(fill=Group, y=(Percentage)^2, x=as.character(State) )) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values =  palette)
dev.off()


#######################
library(fpc)
library(dbscan)
library(RColorBrewer)
library(randomcoloR)

Pseudotime_data_table_EAS <- read.table("/Users/weiz/Desktop/R_code/e11/Pseudotime_data_table_EAS.txt",header = T)
dbscan_result <- fpc::dbscan(Pseudotime_data_table_EAS[,c(2,3)], eps =0.5, MinPts = 1)

Pseudotime_data_table_EAS <-cbind(Pseudotime_data_table_EAS,dbscan=dbscan_result$cluster)
n <- length(unique(dbscan_result$cluster ))
palette <- distinctColorPalette(n)

ggplot(Pseudotime_data_table_EAS,mapping = aes(x= dim1, y= dim2, color= as.character(dbscan_result$cluster ) ))+
  geom_point(size = 1)+
  scale_color_manual(values =  palette)


means_group <- Pseudotime_data_table_EAS %>% 
  group_by(dbscan)%>%
  summarise(mean_t1 = mean(dim1),
            mean_t2 = mean(dim2))%>% print(n=1)
setwd("/Users/weiz/Desktop/R_code/e11/new/")
require("ggrepel")
pdf("dbscan_Pseudotime_data_table_EASeps_0.5_MinPts_1plot_1.pdf")
ggplot(Pseudotime_data_table_EAS, aes(x=dim1, y=dim2,colour=as.character(Pseudotime_data_table_EAS$dbscan)))+
  geom_point(size=2)+
  scale_color_manual(values =  palette)+
  geom_label_repel(data = means_group, aes(x = mean_t1,y = mean_t2),label =means_group$dbscan, colour=1,size = 10,alpha = .5)+ theme(legend.title=element_blank())

dev.off()




######
E11_dbscan_table <- data.frame(matrix(NA, ncol = 11, nrow = 0)) 
colnames(E11_dbscan_table) <-c("African","European","India","Malay","SouthChineseDai","SouthwestChineseYi","EastChinese","Japanese","NorthChineseOroqen","Yakut","American")

E11_stable_all<-  data.frame(matrix(NA, ncol = 11, nrow = 0)) 
colnames(E11_stable_all) <-c("African","European","India","Malay","SouthChineseDai","SouthwestChineseYi","EastChinese","Japanese","NorthChineseOroqen","Yakut","American")

for (dbscan_num in unique(Pseudotime_data_table$dbscan))
  
{
  #dbscan_num <- "1"
  
  E11_table <-  data.frame(matrix(NA, ncol = 11, nrow = 0)) 
  colnames(E11_table) <-c("African","European","India","Malay","SouthChineseDai","SouthwestChineseYi","EastChinese","Japanese","NorthChineseOroqen","Yakut","American")
  
  grouped_person<-Pseudotime_data_table[which(Pseudotime_data_table$dbscan==dbscan_num),1]
  print(dbscan_num)
  for(person_name in grouped_person)
  {
    print(person_name)
    #person_name<-"HG02046"
    genotype <- read.delim(file = paste(person_name,"_my_23andme.1.txt",sep = ""),col.names = c("rsid", "chr", "site","genotype"),row.names = NULL,header = F)
    
    # Use E11
    res <- tfrdpub(genotype, 11, e11.alleles, e11.11.F)
    # Use E11
    ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "E11")
    E11_table[person_name,]<-ances$q
  }
  write.table(E11_table,paste("E11_table_dbscan",dbscan_num,sep = "",".txt"),quote = F,sep = "\t")
  E11_table_sum <- apply(E11_table,2,sum)
  E11_stable_all<-rbind(E11_stable_all,E11_table)
  E11_dbscan_table[dbscan_num,]<-E11_table_sum
  
}


write.table(E11_stable_all,"/Users/weiz/Desktop/R_code/e11/e11.alleles_e11_E11_E11_E11_stable_all.txt",quote = F,sep = "\t")


E11_dbscan_table_reshaped<- melt(t(na.omit(E11_dbscan_table)))

colnames(E11_dbscan_table_reshaped) <- c("Group","dbscan","Percentage")

n <- length(unique(Pseudotime_data_table$dbscan))
palette <- distinctColorPalette(11)
pdf("E11_dbscan_table_reshaped.pdf")
ggplot(E11_dbscan_table_reshaped, aes(fill=Group, y=(Percentage)^2, x=as.character(dbscan) )) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values =  palette)
dev.off()


######



E11_stable_all<-read.table("/Users/weiz/Desktop/R_code/e11/e11.alleles_e11_E11_E11_E11_stable_all.txt",head=T)


plotPseudotime_data_table_EAS<-Pseudotime_data_table_EAS[match(rownames(E11_stable_all),Pseudotime_data_table_EAS$sample_name),]

plotPseudotime_data_table_EAS_dim<- cbind(plotPseudotime_data_table_EAS,E11_stable_all)

plotPseudotime_data_table_EAS_dim<-as.data.frame(plotPseudotime_data_table_EAS_dim)



for (area_e11 in colnames(E11_stable_all))
  
{
  
  pdf(paste("Human trajectory plot ofplot_e11_each_",area_e11,".pdf",sep=""))
  plot<-ggplot(plotPseudotime_data_table_EAS_dim, aes(x=dim1, y=dim2,colour=plotPseudotime_data_table_EAS_dim[,area_e11]))+
    geom_point(size=3)+scale_colour_gradientn(colours = c("#253494","#ffffcc", "#dd1c77"),values = scales::rescale(c(0,0.5)))+labs(title=paste0("E11 ",area_e11))+ theme(plot.title=element_text(size = 30,hjust=0.5),legend.title=element_blank())
  
  
  print(plot)
  dev.off()
}






write.table(E11_dbscan_table,"/Users/weiz/Desktop/R_code/e11/e11.alleles_e11_E11_E11_dbscan_table.txt",quote = F,sep = "\t")

load("/Users/weiz/Desktop/R_code/e11/radmixture-master/data/K12b.12.F.RData")
load("/Users/weiz/Desktop/R_code/e11/radmixture-master/data/K12b.alleles.RData")

setwd("/Users/weiz/Desktop/R_code/e11/23andme/")

human_name<-data.frame(file_name=list.files("/Users/weiz/Desktop/R_code/e11/23andme/"))


human_name_inter<-intersect(Pseudotime_data_table_EAS$sample_name,substring(human_name$file_name, 1, 7))
rownames(Pseudotime_data_table_EAS) <-Pseudotime_data_table_EAS$sample_name
Pseudotime_data_table<-Pseudotime_data_table_EAS[human_name_inter,]


K12B_dbscan_table <- data.frame(matrix(NA, ncol = 12, nrow = 0)) 
colnames(K12B_dbscan_table) <-c("Gedrosia","Siberian","Northwest_African","Southeast_Asian","Atlantic_Med","North_European","South_Asian","East_African","Southwest_Asian","East_Asian","Caucasus","Sub_Saharan")

K12B_stable_all<-  data.frame(matrix(NA, ncol = 12, nrow = 0)) 
colnames(K12B_stable_all) <-c("Gedrosia","Siberian","Northwest_African","Southeast_Asian","Atlantic_Med","North_European","South_Asian","East_African","Southwest_Asian","East_Asian","Caucasus","Sub_Saharan")

for (dbscan_num in unique(Pseudotime_data_table$dbscan))
  
{
  #dbscan_num <- "1"
  
  K12B_table <-  data.frame(matrix(NA, ncol = 12, nrow = 0)) 
  colnames(K12B_table) <-c("Gedrosia","Siberian","Northwest_African","Southeast_Asian","Atlantic_Med","North_European","South_Asian","East_African","Southwest_Asian","East_Asian","Caucasus","Sub_Saharan")
  
  grouped_person<-Pseudotime_data_table[which(Pseudotime_data_table$dbscan==dbscan_num),1]
  print(dbscan_num)
  for(person_name in grouped_person)
  {
    print(person_name)
    #person_name<-"HG02046"
    genotype <- read.delim(file = paste(person_name,"_my_23andme.1.txt",sep = ""),col.names = c("rsid", "chr", "site","genotype"),row.names = NULL,header = F)
    
    # Use K12B
    res <- tfrdpub(genotype, 12, K12b.alleles, K12b.12.F)
    # Use K12B
    ances <-  fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "K12b")
    K12B_table[person_name,]<-ances$q
  }
  write.table(K12B_table,paste("K12B_table_dbscan",dbscan_num,sep = "",".txt"),quote = F,sep = "\t")
  K12B_table_sum <- apply(K12B_table,2,sum)
  K12B_stable_all<-rbind(K12B_stable_all,K12B_table)
  K12B_dbscan_table[dbscan_num,]<-K12B_table_sum
  
}




plotPseudotime_data_table_EAS<-Pseudotime_data_table_EAS[match(rownames(K12B_stable_all),Pseudotime_data_table_EAS$sample_name),]

plotPseudotime_data_table_EAS_dim<- cbind(plotPseudotime_data_table_EAS,K12B_stable_all)

plotPseudotime_data_table_EAS_dim<-as.data.frame(plotPseudotime_data_table_EAS_dim)



for (area_K12B in colnames(K12B_stable_all))
  
{
  
  pdf(paste("Human trajectory plot ofplot_K12B_each_",area_K12B,".pdf",sep=""))
  plot<-ggplot(plotPseudotime_data_table_EAS_dim, aes(x=dim1, y=dim2,colour=(plotPseudotime_data_table_EAS_dim[,area_K12B])^2))+
    geom_point(size=3)+scale_colour_gradientn(colours = c("#253494","#ffffcc", "#dd1c77"),values = scales::rescale(c(0,0.5)))+labs(title=paste0("K12B ",area_K12B))+ theme(plot.title=element_text(size = 30,hjust=0.5),legend.title=element_blank())
  
  
  print(plot)
  dev.off()
}



plotbind_population_umap_DDR<-bind_population_umap_DDR[match(rownames(K12B_stable_all),bind_population_umap_DDR$Sample.name),]

plotbind_population_umap_DDR_dim<- cbind(plotbind_population_umap_DDR,K12B_stable_all)

plotbind_population_umap_DDR_dim<-as.data.frame(plotbind_population_umap_DDR_dim)



for (area_K12B in colnames(K12B_stable_all))
  
{
  
  pdf(paste("EAS ddrtree plot ofplot_K12B_each_",area_K12B,".pdf",sep=""))
  plot<-ggplot(plotbind_population_umap_DDR_dim, aes(x=Z1, y=Z2,colour=(plotbind_population_umap_DDR_dim[,area_K12B])))+
    geom_point(size=3)+scale_colour_gradientn(colours = c("#253494","#ffffcc", "#dd1c77"))+labs(title=paste0("K12B ",area_K12B))+ theme(plot.title=element_text(size = 30,hjust=0.5),legend.title=element_blank())
  
  
  print(plot)
  dev.off()
}


