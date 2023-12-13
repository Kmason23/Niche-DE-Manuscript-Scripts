library(Matrix)
library(Seurat)
library(ggplot2)
library(abind)
library(enrichR)
library(nicheDE)

#read in niche-DE object
NDE = readRDS("integrated_2871_2873_nicheDE.rds")

#get niche-DE genes in fibroblasts when near PT cells
NDE_genes = get_niche_DE_genes(NDE,'I',index="Fibroblast",niche = "PT",pos = F,alpha = 0.05)
View(NDE_genes)

#get ligand-receptor processes in fibroblasts when near PT cells
NDE_lr = niche_LR_spot(NDE,ligand_cell = "PT",receptor_cell = "Fibroblast",
                       ligand_target_matrix = niche_net_ligand_target_matrix,
                       lr_mat = ramilowski_ligand_receptor_list, truncation_value = 5)
View(NDE_lr)

#perform pathway enrichment analysis in fibroblasts when near PT cells
NDE_processes = enrichr(NDE_genes[,1],databases = 'Reactome_2016')
View(NDE_processes$Reactome_2016)

