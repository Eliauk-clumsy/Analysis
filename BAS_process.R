LiverBAs <- read.csv('liverBAs.csv', header = TRUE,row.names=1,check.names = F)[-11,]
#PLS-DA分析----
library(mixOmics)
library(MASS)
library(lattice)
library(ggplot2)

YY<-rep(c("CON","LKO"),each=5)

png( 
  filename = "name.png", 
  width = 1400,           
  height = 1000,          
  units = "px",          
  bg = "transparent",          
  res = 300)             

plsda.datatm <-plsda(LiverBAs, YY, ncomp = 2)
plotIndiv(plsda.datatm, 
          ind.names = F,
          legend=TRUE,
          ellipse =TRUE,
          point.lwd = 1,
          style = "ggplot2",
          col=c("red","blue"),
          pch=16,
          title = "PLS-DA",
          legend.title = "Group")



background = background.predict(plsda.datatm, comp.predicted=2, dist = "max.dist") 
plotIndiv(plsda.datatm, comp = 1:2,
          group = YY, ind.names = F,
          legend = TRUE,  background = background,
          legend.title = "",title = "PLS-DA",size.xlabel = rel(1.3),
          size.ylabel = rel(1.3),  
          size.legend = rel(1.2),point.lwd = 1.5, col=c("red","steelblue3")
          )

dev.off()

library(ggplotify)
as.ggplot() #

#热图分析----
library(pheatmap)
library(ggplot2)
library(ggplotify)

Exp<-data.matrix(t(LiverBAs))

group_list<- data.frame(rep(c("CON","LKO"),each=5))

rownames(group_list)<- colnames(Exp)
colnames(group_list) <- "group"


heatmap_plot <- as.ggplot(pheatmap(Exp,scale="row",cluster_rows = T,
                                   cellwidth = 12, cellheight = 12,
                                   treeheight_row = 30, treeheight_col = 30,
                                   legend = T,cluster_cols = F,
                                   color = colorRampPalette(colors = c("blue","white","red"))(100),
                                   annotation_col =group_list,show_colnames = F)
)
heatmap_plot

ggsave(heatmap_plot,file="heatmap.tiff",
       width=12,height = 24,units ="cm",
       dpi = 300)

#堆叠柱状图----
library(plyr)
library(ggplot2)
library(RImagePalette)
library(imager)
library(scales)
library(ggplot2)
library(ggprism)
library(reshape)
library(ggalluvial)
library("ggsci")
library("ggplot2")
library("gridExtra")

Liver <- read.csv('liverBAs.csv', header = TRUE,check.names = F)[-11,]
colnames(Liver)[1]="samples"
Liver[,1]=paste0(c(rep(c("CON","LKO"),each=5)),1:5)

#换算百分比,各小组----
L1 <- melt(Liver,id.vars = 'samples')
pL1 <- ddply(L1,"samples",transform,percentage=value/sum(value)*100)
names(pL1)[1:2] <- c("group","X")  #修改列名
colnames(pL1)
write.csv(pL1,file = "长数据百分比.csv")

p1<-ggplot(pL1, aes( x = group,y= percentage,fill = X,
                      stratum = X, alluvium = X))+
  geom_stratum(width = 0.7, color='white')+
  geom_alluvium(alpha = 0.5,
                width = 0.7,
                color='white',
                size = 1,
                curve_type = "linear")+
  scale_y_continuous(expand = c(0,0))+
  labs(x="",y="Relative Abundance(%)",
       fill="group")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16, 
              base_line_size = 0.8, 
              axis_text_angle = 45)+ 
  theme(legend.position = 'top') 
p1

ggsave(p,file="各小鼠百分比图.tiff",
       width=24,height = 24,units ="cm",
       dpi = 300)

#换算百分比,大组----
avLiver <- read.csv('average_liverBAs.csv', header = TRUE,check.names = F)
colnames(avLiver)[1]="samples"
L2 <- melt(avLiver,id.vars = 'samples')
pL2 <- ddply(L2,"samples",transform,percentage=value/sum(value)*100)
names(pL2)[1:2] <- c("group","X")  #修改列名
colnames(pL2)
write.csv(pL2,file = "大组长数据百分比.csv")

p2<-ggplot(pL2, aes( x = group,y= percentage,fill = X,
                    stratum = X, alluvium = X))+
  geom_stratum(width = 0.7, color='white')+
  geom_alluvium(alpha = 0.5,
                width = 0.7,
                color='white',
                size = 1,
                curve_type = "linear")+
  scale_y_continuous(expand = c(0,0))+
  labs(x="",y="Relative Abundance(%)",
       fill="group")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16, 
              base_line_size = 0.8, 
              axis_text_angle = 0)+ 
  theme(legend.position = 'top') 
p2


ggsave(p1,file="组百分比图.tiff",
       width=24,height = 24,units ="cm",
       dpi = 300)

#换算百分比,各小组，结合与未结合----
conLiver <- read.csv('conliverBAs.csv', header = TRUE,check.names = F)
colnames(conLiver)[1]="samples"
L3 <- melt(conLiver,id.vars = 'samples')
pL3 <- ddply(L3,"samples",transform,percentage=value/sum(value)*100)
names(pL3)[1:2] <- c("group","X")  #修改列名
colnames(pL3)
#write.csv(pL3,file = "长数据百分比.csv")

p3<-ggplot(pL3, aes( x = group,y= percentage,fill = X,
                     stratum = X, alluvium = X))+
  geom_stratum(width = 0.7, color='white')+
  geom_alluvium(alpha = 0.5,
                width = 0.7,
                color='white',
                size = 1,
                curve_type = "linear")+
  scale_y_continuous(expand = c(0,0))+
  labs(x="",y="Relative Abundance(%)",
       fill="group")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16, 
              base_line_size = 0.8, 
              axis_text_angle = 45)+ 
  theme(legend.position = 'top') 
p3


ggsave(p3,file="组百分比图_结合未结合.tiff",
       width=12,height = 12,units ="cm",
       dpi = 300)

#换算百分比,大组，结合与未结合----
conbigLiver <- read.csv('conbigliverBAs.csv', header = TRUE,check.names = F)
colnames(conbigLiver)[1]="samples"
L4 <- melt(conbigLiver,id.vars = 'samples')
pL4 <- ddply(L4,"samples",transform,percentage=value/sum(value)*100)
names(pL4)[1:2] <- c("group","X")  #修改列名
colnames(pL4)
#write.csv(pL4,file = "长数据百分比.csv")

p4<-ggplot(pL4, aes( x = group,y= percentage,fill = X,
                     stratum = X, alluvium = X))+
  geom_stratum(width = 0.7, color='white')+
  geom_alluvium(alpha = 0.5,
                width = 0.7,
                color='white',
                size = 1,
                curve_type = "linear")+
  scale_y_continuous(expand = c(0,0))+
  labs(x="",y="Relative Abundance(%)",
       fill="group")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16, 
              base_line_size = 0.8, 
              axis_text_angle = 45)+ 
  theme(legend.position = 'top') 
p4


ggsave(p4,file="大组百分比图_结合未结合.tiff",
       width=12,height = 12,units ="cm",
       dpi = 300)


