library(nicheDE)
library(ggplot2)
#make niche-DE object
Merged_NDE = readRDS('NDE_liver1_subclone.rds')

#Figure 4E)get subclone plots
coord = Merged_NDE@coord
R = matrix(c(0,-1,1,0),2,2)
coord = t(R%*%as.matrix(t(coord)))
colnames(coord) = c('x','y')

subclone_0 = which(Merged_NDE@num_cells[,2]>0)
subclone_1 = which(Merged_NDE@num_cells[,5]>0)
subclone_2 = which(Merged_NDE@num_cells[,6]>0)
colors = rep('0',nrow(coord))
D = rep(1,nrow(coord))
D[c(subclone_0,subclone_1,subclone_2)] = 1
colors[subclone_0] = '0'
colors[subclone_1] = 'subclone_1'
colors[subclone_2] = 'subclone_2'

coord = data.frame(coord,D,colors)
ggplot(coord,aes(x,y,size=ifelse(D==0,NA, D),color = colors))+geom_point(size = 1.5)+
  scale_color_manual(colors, values = c('grey',"orange", "blue"))



# Subclone 1 marker genes in fibroblasts realtive to subclone 2
subclone1_marker = niche_DE_markers(Merged_NDE,'stromal','subclone_1','subclone_2',0.05)
View(subclone1_marker)

#Subclone 2 marker genes
subclone2_marker = niche_DE_markers(Merged_NDE,'stromal','subclone_2','subclone_1',0.05)
View(subclone2_marker)

#Figure 4H) Get processes corresponding to subclone 1 marker genes
subclone1_processes = enrichr(subclone1_marker[,1],databases = 'Reactome_2016')
View(subclone1_processes$Reactome_2016)

#Figure 4G) Get processes corresponding to subclone 2 marker genes
subclone2_processes = enrichr(subclone2_marker[,1],databases = 'Reactome_2016')
View(subclone2_processes$Reactome_2016)





