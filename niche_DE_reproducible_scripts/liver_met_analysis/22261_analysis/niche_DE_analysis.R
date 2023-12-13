library(nicheDE)
#read in niche-DE object
NDE = readRDS("22261_NDE.rds")
#get (tumor,fibroblast)- niche-DE genes
tumor_stromal_neg = get_niche_DE_genes(NDE,test.level = "I",index = "epithelial",niche = "stromal", positive = F)


