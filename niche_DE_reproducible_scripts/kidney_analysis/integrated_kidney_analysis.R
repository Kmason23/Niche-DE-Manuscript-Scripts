library(Matrix)
library(Seurat)
library(ggplot2)
library(abind)
library(enrichR)
library(nicheDE)

#read in sample names
sample = c("HK_2871","HK_2873")
#read niche-DE objects
NDE_2871 = readRDS("kidney_2871_nicheDE.rds")
NDE_2873 = readRDS("kidney_2873_nicheDE.rds")
#merge objects
NDE = MergeObjects(list(NDE_1 = NDE_2871,NDE_2 = NDE_2873))
#calculate effective niche
NDE = CalculateEffectiveNiche(NDE)

#register cluster
cl <- parallel::makeCluster(4,outfile = "outfile.txt")
doParallel::registerDoParallel(cl)

#run niche-DE
#start time of funtion
start_time <- Sys.time()
#run niche-De
NDE = niche_DE(NDE,cluster = cl,C = 300,M = 20,gamma = 0.8,print = T,Int = T,batch = T,self_EN = F)
#end time of function
end_time <- Sys.time()
#print how long it takes to run niche-DE
print(end_time-start_time)
#close cluster
doParallel::stopImplicitCluster()

saveRDS(NDE,"integrated_kidney_NDE.rds")

