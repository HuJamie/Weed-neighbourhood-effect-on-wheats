####################08/July/2021, Jie Hu######################
######Venn diagram for shared clusters of weeds to wheat######
library(VennDiagram)
library(grid)
library(futile.logger)
library(lambda.r)
#library(parser) ## probably, we do not need this package
ngs.all<-as.data.frame(rbind(colSums(ngs.soil),colSums(ngs.root[,-1])))

ngs.venn<-colnames(ngs.all) ## res.up is your OTU table
ngs.venn1<-t(ngs.venn)

list1<-mat.or.vec(467,1) ##96 is the number of changed clusters, you can change it based on your data
list2<-mat.or.vec(467,1)

for (i in 1:dim(ngs.all)[2])
{
  if (ngs.all[1,i]!=0) list1[i]<-c(ngs.venn1[,i])
  if (ngs.all[2,i]!=0) list2[i]<-c(ngs.venn1[,i])
}

# Prepare a palette of colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel1") ## three is the group, can change it based on how many groups you have

venn.plot<-venn.diagram(
  x=list(list1,list2),
  category.names=c("Soil mycobiota","Wheat root mycobiota"),
  filename='#Venn_soil to root.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1500 , 
  width = 1600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = 'blue',
  
  # Numbers
  cex = .8,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135,140),
  cat.dist = c(0.055, 0.055, 0.085,0.095),
  cat.fontfamily = "sans",
  rotation = 1
)
