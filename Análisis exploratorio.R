###
#Cachexia is a complex metabolic syndrome associated with an underlying illness (such as cancer) 
#and characterized by loss of muscle with or without loss of fat mass (Evans et al., 2008). 
#A total of 77 urine samples were collected being 47 of them patients with cachexia, 
#and 30 control patients (from the "specmine.datasets" R package)
###
#Eisner et al. 2010 (Learning to predict cancer-associated skeletal muscle wasting from 1H-NMR profiles of urinary metabolites)
#H-NMR spectra acquisition (= proton nuclear magnetic resonance spectroscopy) and targeted profiling of metabolites from urine samples
#Data are concentrations of these metabolites
#

#######
cachexiadata <- read.csv("https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2024-Cachexia/human_cachexia.csv")
head(cachexiadata, c(6,5))

dim(cachexiadata)

colnames(cachexiadata)

if(!(require("SummarizedExperiment"))){
  BiocManager::install("SummarizedExperiment")
}
library(SummarizedExperiment, quietly = T)


valores <- t(as.matrix(data.frame(cachexiadata[,-c(1,2)], row.names = cachexiadata[,1])))

columnasdatos <- data.frame(Group=cachexiadata[,2], row.names = cachexiadata[,1])

sumexpcachexia <- SummarizedExperiment(assays = valores, colData = columnasdatos)


if(!(require("POMA"))){
  BiocManager::install("POMA")
}
library(POMA, quietly = T)

valores <- data.frame(cachexiadata[,-c(1,2)], row.names = cachexiadata[,1])

columnasdatos <- data.frame(Patient.ID=cachexiadata[,1], Group=cachexiadata[,2])

sumexpcachexia2 <- PomaCreateObject(metadata = columnasdatos, features = valores)

sumexpcachexia2

head(assay(sumexpcachexia2), c(6,5))

head(assay(sumexpcachexia), c(6,5))

rm(list=c("cachexiadata","valores", "columnasdatos", "sumexpcachexia2"))

# NAs
if(mean(is.na(assay(sumexpcachexia)))==0){
  print("No hay missing values")
} else{
  cachexiamod <- PomaImpute(sumexpcachexia)
}


boxplot(assay(sumexpcachexia), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n", xlab="Patients", ylab="Metabolite values",
        main="Metabolite values for 77 patients in two groups")
axis(1, labels=F, at=c(1:77))

par(mfrow=c(2,2))
for(mod_val in c("auto_scaling","log_pareto","log_scaling","log")){
  boxplot(assay(PomaNorm(sumexpcachexia, method = mod_val)), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n",
          main=mod_val)}

par(mfrow=c(1,2))
boxplot(log2(assay(sumexpcachexia)), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n",
        main="Log-transformed values")
boxplot(scale(log2(assay(sumexpcachexia))), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n",
        main="Normalized log-transformed values")

par(mfrow=c(1,1))

sumexpcachexia_log <- PomaNorm(sumexpcachexia, method = "log")
#boxplot(assay(sumexpcachexia_log), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n",
        #main="Log-transformed values")

sumexpcachexia_log_norm <- sumexpcachexia_log
assay(sumexpcachexia_log_norm) <- scale(assay(sumexpcachexia_log_norm))
#boxplot(assay(sumexpcachexia_log_norm), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n",
        #main="Normalized log-transformed values")


#head(t(scale(assay(sumexpcachexia_log))),c(5,5))
#head(t(assay(sumexpcachexia_log_norm)),c(5,5))

#sumexpcachexia_log <- PomaNorm(sumexpcachexia, method = "log")

#head(assay(sumexpcachexia_log),c(6,5))
#head(log(assay(sumexpcachexia)),c(6,5))

#prcomplog2 <- prcomp(t(assay(sumexpcachexia_log_norm)))
#plot(prcomplog2$x[,1], prcomplog2$x[,2])

#prcomplog2 <- prcomp(t(assay(sumexpcachexia_log_norm)), scale. = F)
#plot(prcomplog2$x[,1], prcomplog2$x[,2])
#text(prcomplog2$x[,1], prcomplog2$x[,2], colData(sumexpcachexia_log_norm)[,1])
round(prcomplog2$sdev^2 / sum(prcomplog2$sdev^2) *100, 2)



library(mixOmics, quietly = T)

acp <- tune.pca(t(assay(sumexpcachexia_log_norm)), scale = F, ncomp = 20)
par(mfrow=c(1,2))
plot(acp, main="Proportion of explained variance for each PC")
plot(acp, main="Proportion of explained variance for each PC", type="b")
par(mfrow=c(1,1))

round(acp$cum.var*100,2)
round(acp$prop_expl_var$X*100,2)

acp <- tune.pca(t(assay(sumexpcachexia_log_norm)), scale = F, ncomp = 3)
round(acp$prop_expl_var$X*100,2)

par(mfrow=c(1,2))
plotIndiv(acp, ind.names = F, group = colData(sumexpcachexia_log_norm)[,1], ellipse = T, title = "ACP de 63 metabolitos para dos grupos de pacientes",
          legend = T, legend.title = "Grupo")
plotIndiv(acp, comp = c(2,3), ind.names = F, group = colData(sumexpcachexia_log_norm)[,1], ellipse = T, title = "ACP de 63 metabolitos para dos grupos de pacientes",
          legend = T, legend.title = "Grupo")


plotIndiv(acp, ind.names = F, group = colData(sumexpcachexia_log_norm)[,1], title = "ACP de 63 metabolitos para dos grupos de pacientes",
          legend = T, legend.title = "Grupo")
s.chull(acp$variates$X[,1:2], fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)

par(mfrow=c(1,1))
plot(acp$variates$X[,1:2])
s.chull(acp$variates$X[,1:2], fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)


plotVar(acp, cex = 3, title = "Contribución de los metabolitos en los componentes y correlaciones entre metabolitos")
plotVar(acp, comp = c(2,3), cex = 3, title = "Contribución de los metabolitos en los componentes y correlaciones entre metabolitos")

biplot(acp, group = colData(sumexpcachexia_log_norm)[,1], ind.names = F, var.names.size = 3, legend.title = "Grupo")
biplot(acp, comp = c(2,3), group = colData(sumexpcachexia_log_norm)[,1], ind.names = F, var.names.size = 3, legend.title = "Grupo")

selectVar(acp, comp = 1)$value
plotLoadings(acp, comp = 1)

selectVar(acp, comp = 2)$value
plotLoadings(acp, comp = 2)



acp_mod <- tune.spca(t(assay(sumexpcachexia_log_norm)), scale = F, ncomp = 3, 
                              folds = 5, 
                              test.keepX = 1:nrow(assay(sumexpcachexia_log_norm)), nrepeat = 25) 
acp_mod$choice.keepX

par(mfrow=c(1,1))
plot(acp_mod)

abs(selectVar(acp, comp = 1)$value)
select1 <- length((abs(selectVar(acp, comp = 1)$value)>0.3)[(abs(selectVar(acp, comp = 1)$value)>0.3)==T])

abs(selectVar(acp, comp = 2)$value)
select2 <- length((abs(selectVar(acp, comp = 2)$value)>0.3)[(abs(selectVar(acp, comp = 2)$value)>0.3)==T])

final.spca.multi <- spca(t(assay(sumexpcachexia_log_norm)), scale = F, keepX = c(select1, select2))
final.spca.multi$prop_expl_var$X * 100

plotIndiv(final.spca.multi, ind.names = F, group = colData(sumexpcachexia_log_norm)[,1], ellipse = T, title = "ACP de 63 metabolitos para dos grupos de pacientes",
          legend = T, legend.title = "Grupo")


dist_euc <- dist(t(assay(sumexpcachexia_log_norm)))
#heatmap(as.matrix(dist_euc), RowSideColors = rep(c("steelblue","tomato"),c(47,30)), 
         #ColSideColors = rep(c("steelblue","tomato"),c(47,30)))

dist_cor <- as.dist(1-cor(assay(sumexpcachexia_log_norm)))
#heatmap(as.matrix(dist_cor), RowSideColors = rep(c("steelblue","tomato"),c(47,30)), 
        #ColSideColors = rep(c("steelblue","tomato"),c(47,30)))


library(MASS, quietly = T)
samm <- sammon(dist_euc, trace=FALSE)
plot(samm$points, pch = rep(c(19,17),c(47,30)), col = rep(c("steelblue","tomato"),c(47,30)), 
     xlab="Component 1", ylab="Component 2", main="Sammon's non-linear mapping based on metabolite values for patients")
library(ade4, quietly = T)
s.chull(samm$points, fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("bottomleft", legend=unique(colData(sumexpcachexia_log_norm)[,1]), pch=c(19,17), col = c("steelblue","tomato"),title="Grupo")

samm <- sammon(dist_cor, trace=FALSE)
plot(samm$points, pch = rep(c(19,17),c(47,30)), col = rep(c("steelblue","tomato"),c(47,30)), 
     xlab="Component 1", ylab="Component 2", main="Sammon's non-linear mapping based on metabolite values for patients")
library(ade4, quietly = T)
s.chull(samm$points, fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("bottomleft", legend=unique(colData(sumexpcachexia_log_norm)[,1]), pch=c(19,17), col = c("steelblue","tomato"),title="Grupo")



mds <- cmdscale(dist_euc, k = 10, eig = TRUE)
valeigdist <- round(mds$eig^2 / sum(mds$eig^2) * 100, 2)
plot(mds$points, pch = rep(c(19,17),c(47,30)), col = rep(c("steelblue","tomato"),c(47,30)), 
     xlab=paste("Component 1", paste0("(",valeigdist[1], "%)")), ylab=paste("Component 2", paste0("(",valeigdist[2], "%)")), 
     main="Classical MDS based on metabolite values for patients")
s.chull(mds$points, fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("bottomright", legend=unique(colData(sumexpcachexia_log_norm)[,1]), pch=c(19,17), col = c("steelblue","tomato"),title="Grupo")

mds <- cmdscale(dist_cor, k = 10, eig = TRUE)
valeigdist <- round(mds$eig^2 / sum(mds$eig^2) * 100, 2)
plot(mds$points, pch = rep(c(19,17),c(47,30)), col = rep(c("steelblue","tomato"),c(47,30)), 
     xlab=paste("Component 1", paste0("(",valeigdist[1], "%)")), ylab=paste("Component 2", paste0("(",valeigdist[2], "%)")), 
     main="Classical MDS based on metabolite values for patients")
s.chull(mds$points, fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("bottomright", legend=unique(colData(sumexpcachexia_log_norm)[,1]), pch=c(19,17), col = c("steelblue","tomato"),title="Grupo")





plot(hclust(dist_euc,method="average"), labels = colData(sumexpcachexia_log_norm)[,1], hang=-1)

plot(hclust(dist_cor,method="average"), labels = colData(sumexpcachexia_log_norm)[,1], hang=-1)

###
metabolitos <- assay(sumexpcachexia_log_norm)

within_ss <- c()
for (i in 1:10) {
  within_ss[i] <- sum(kmeans(metabolitos, i, iter.max = 1000)$withinss)
}

plot(1:10,within_ss, type="b", xlab="Número de clusters",
     ylab="Suma de los cuadrados intra-grupos")

kmeans <- kmeans(metabolitos,4, iter.max = 1000)

library(factoextra, quietly = T)
fviz_cluster(kmeans, data = metabolitos, palette = c("steelblue", "tomato", "forestgreen", "brown"), 
             ellipse.type = "euclid", star.plot = TRUE)


library(made4, quietly = T)
datos <- assay(sumexpcachexia_log_norm)
colnames(datos) <- colData(sumexpcachexia_log_norm)[,1]

overview(datos)
k.coa <- ord(datos, type="coa")
plot.ord(k.coa)
plot.ord(k.coa, classvec = colData(sumexpcachexia_log_norm)[,1])