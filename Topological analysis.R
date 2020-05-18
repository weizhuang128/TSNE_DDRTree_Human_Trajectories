
library("alphahull")
library("adehabitatHR")
library('maptools')
library(RCurl)



g_1000<-read.delim("/Users/weiz/Downloads/recoded_1000G.raw",header = T,sep = " ")
rownames(g_1000) <- g_1000[,2]
unprocessed_single_cell_data_total_transposed <-g_1000[,7:39188]
population<- read.delim("/Users/weiz/Downloads/igsr_samples.tsv",header = T)
population$Sample.name
population<-as.data.frame(population[which(population$Superpopulation.code=="EAS"),])

area_unprocessed_single_cell_data_total_transposed<-intersect(population$Sample.name,rownames(unprocessed_single_cell_data_total_transposed))

g_1000_area<-unprocessed_single_cell_data_total_transposed[area_unprocessed_single_cell_data_total_transposed,]


datadata<-as.data.frame(as.character(colnames(g_1000_area)) )
mutation_genotype<-datadata %>% separate("as.character(colnames(g_1000_area))", into=c("ID", "genotype"), "_")

colnames(g_1000_area) <-  mutation_genotype$ID


bind_population_umap<-read.delim("/Users/weiz/Desktop/R_code/e11/new/sup/bind_population_umap_DDR_EAS_.txt",header = T)

SNP_list <- unlist(read.delim("/Users/weiz4/Desktop/R_code/e11/SNP_List_500",header = F))
pdf("/Users/weiz/Desktop/R_code/e11/SNP_List_500_plot.pdf")

for (SNP in SNP_list)
  
{
  
  #SNP <- "rs6676961"
print(ggplot(bind_population_umap, aes(x=Z1, y=Z2,colour=as.factor( g_1000_area[,SNP])))+
    geom_point()+labs(title=paste0("E11 ",SNP))+theme(plot.title=element_text(size = 30,hjust=0.5),legend.title=element_blank()))
  
}
dev.off()


ashape.obj <- ashape(bind_population_umap[,c("Z1","Z2")], alpha = 2)

#Plotting the alpha
plot(ashape.obj, wlines= "del", col = c(4, 1, 2, 3), lwd = c(1,3))
plot(ashape.obj, wlines= "vor", col = c(4, 1, 2, 3), lwd = c(1,3))
plot(ashape.obj)

#home area analysis
xyt<-subset(bind_population_umap, select = c(Z1,Z2))
id<-data.frame(g_1000_area[,SNP])
locs1<-id
coordinates(locs1)<-xyt
class(locs1)
plot(locs1, col=as.data.frame(locs1)[,1]) 

cp<-mcp(locs1[,1], percent=50)
plot(cp)
plot(locs1, col=as.data.frame(locs1)[,1], add=T)
cp
kud<-kernelUD(locs1[,1])
homerange <- getverticeshr(kud,percent = 20)
class(homerange)
plot(homerange,col=1:3) 
kde.areas<- kernel.area(kud, percent=c(50,95))
kde.areas
vud<-getvolumeUD(kud)
image(vud[[1]],col = hcl.colors(1000, "terrain"))
xyzv <- as.image.SpatialGridDataFrame(vud[[2]])
contour(xyzv, add=TRUE)
image(vud[[2]], add=TRUE)
kud1 <- kernelUD(locs1[,1], h="href")
homerange1 <- getverticeshr(kud1, percent = 20)
plot(homerange1,border=1:3, lwd=6) # this better than before, doesn't it?
plot(locs1, col=as.data.frame(locs1)[,1], add=T)
core <- getverticeshr(kud1, percent = 50)
plot(core,border=1:3, lwd=4, lty = "dashed", add = T) 
contour(xyzv, add=TRUE)



