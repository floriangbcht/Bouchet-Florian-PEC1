###
#Cachexia is a complex metabolic syndrome associated with an underlying illness (such as cancer) 
#and characterized by loss of muscle with or without loss of fat mass (Evans et al., 2008). 
#A total of 77 urine samples were collected being 47 of them patients with cachexia, 
#and 30 control patients (from the "specmine.datasets" R package)
###

# El conjunto de datos se puede cargar a partir del paquete 'specmine.datasets' pero primero hay que tener instalado 'xcms' (que es de Bioconductor)

if(!(require("xcms"))){
  BiocManager::install("xcms")
}

library(xcms)

# El paquete 'specmine.datasets' ha sido quitado del repositorio CRAN, pero se puede instalar una versión más antigua a partir de los registros (archivos)
url <- "https://cran.r-project.org/src/contrib/Archive/specmine.datasets/specmine.datasets_0.0.2.tar.gz"
install.packages(url, repos = NULL, type = "source")

data("cachexia", package = "specmine.datasets")
head(cachexia)
cachexiaspecmine <- cachexia
head(cachexiaspecmine)
rm(cachexia)

#######
cachexiadata <- read.csv("https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2024-Cachexia/human_cachexia.csv")
head(cachexiadata, c(6,5))

dim(cachexiadata)
dim(cachexiaspecmine$data)

colnames(cachexiadata)
rownames(cachexiaspecmine$data)

rm(cachexiaspecmine)

if(!(require("SummarizedExperiment"))){
  BiocManager::install("SummarizedExperiment")
}
library(SummarizedExperiment)


if(!(require("POMA"))){
  BiocManager::install("POMA")
}
library(POMA)


valores <- t(as.matrix(data.frame(cachexiadata[,-c(1,2)], row.names = cachexiadata[,1])))

columnasdatos <- data.frame(Group=cachexiadata[,2], row.names = cachexiadata[,1])

sumexpcachexia <- SummarizedExperiment(assays = valores, colData = columnasdatos)

sumexpcachexia


if(!(require("POMA"))){
  BiocManager::install("POMA")
}
library(POMA)

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

head(assay(sumexpcachexia_log),c(6,5))
head(log(assay(sumexpcachexia)),c(6,5))

prcomplog2 <- prcomp(t(assay(sumexpcachexia_log)))
plot(prcomplog2$x[,1], prcomplog2$x[,2])

library(mixOmics)

acp <- tune.pca(t(assay(sumexpcachexia_log)), scale = T, ncomp = 10)
plot(acp, main="Proportion of explained variance for each PC")

round(acp$cum.var*100,2)
round(acp$prop_expl_var$X*100,2)

acp <- tune.pca(t(assay(sumexpcachexia_log)), scale = T, ncomp = 2)
round(acp$prop_expl_var$X*100,2)

plotIndiv(acp, ind.names = F, group = colData(sumexpcachexia_log)[,1], ellipse = T, title = "ACP de 63 metabolitos para dos grupos de pacientes",
          legend = T, legend.title = "Grupo")


plotVar(acp, cex = 3, title = "Contribución de los metabolitos en los componentes y correlaciones entre metabolitos")
plotVar(acp, var.names = F)

biplot(acp, group = colData(sumexpcachexia_log)[,1], ind.names = F, var.names.size = 3, legend.title = "Grupo")

selectVar(acp, comp = 1)$value
plotLoadings(acp, comp = 1)

selectVar(acp, comp = 2)$value
plotLoadings(acp, comp = 2)


set.seed(42)
acp_mod <- tune.spca(t(assay(sumexpcachexia_log)), scale = T, ncomp = 2, 
                              folds = 5, 
                              test.keepX = 10:nrow(assay(sumexpcachexia_log)), nrepeat = 50) 
acp_mod$choice.keepX

plot(acp_mod)


dist_euc <- dist(scale(t(assay(sumexpcachexia_log))))
heatmap(as.matrix(dist_euc), RowSideColors = rep(c("steelblue","tomato"),c(47,30)), 
         ColSideColors = rep(c("steelblue","tomato"),c(47,30)))

heatmap(as.matrix(cor(scale(assay(sumexpcachexia_log)))), RowSideColors = rep(c("steelblue","tomato"),c(47,30)), 
        ColSideColors = rep(c("steelblue","tomato"),c(47,30)))

require(MASS)
samm <- sammon(dist_euc, trace=FALSE)

plot(samm$points, pch = rep(c(19,17),c(47,30)), col = rep(c("steelblue","tomato"),c(47,30)), 
     xlab="Component 1", ylab="Component 2", main="Sammon's non-linear mapping of metabolite values for two patient groups")
library(ade4)
s.chull(samm$points, fac = factor(colData(sumexpcachexia_log)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("bottomright", legend=unique(colData(sumexpcachexia_log)[,1]), pch=c(19,17), col = c("steelblue","tomato"),title="Grupo")


mds <- cmdscale(dist_euc, k = 10, eig = TRUE)
round(mds$eig^2 / sum(mds$eig^2) * 100, 2)
plot(mds$points, pch = rep(c(19,17),c(47,30)), col = rep(c("steelblue","tomato"),c(47,30)), 
     xlab="Component 1", ylab="Component 2", main="Classical MDS of metabolite values for two patient groups")
s.chull(mds$points, fac = factor(colData(sumexpcachexia_log)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("topright", legend=unique(colData(sumexpcachexia_log)[,1]), pch=c(19,17), col = c("steelblue","tomato"),title="Grupo")


clust <- hclust(dist_euc, method="average")
plot(hclust(dist_euc,method="average"), labels = colData(sumexpcachexia_log)[,1], hang=-1)
