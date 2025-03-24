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

columnasdatos <- data.frame(Patient.type=cachexiadata[,2], row.names = cachexiadata[,1])

sumexpcachexia <- SummarizedExperiment(assays = valores, colData = columnasdatos)

sumexpcachexia


if(!(require("POMA"))){
  BiocManager::install("POMA")
}
library(POMA)

valores <- data.frame(cachexiadata[,-c(1,2)], row.names = cachexiadata[,1])

columnasdatos <- data.frame(Patient.ID=cachexiadata[,1], Patient.type=cachexiadata[,2])

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

# 
PomaBoxplots(sumexpcachexia, outcome = "Patient.type")
boxplot(assay(sumexpcachexia), col=rep(c("steelblue","tomato"),c(47,30)))

sumexpcachexia_norm <- PomaNorm(sumexpcachexia, method = "auto_scaling")
boxplot(scale(assay(sumexpcachexia)), col=rep(c("steelblue","tomato"),c(47,30)))
PomaBoxplots(sumexpcachexia_norm, outcome = "Patient.type")
PomaBoxplots(PomaNorm(sumexpcachexia, method = "log"), outcome = "Patient.type")

boxplot(log2(assay(sumexpcachexia)), col=rep(c("steelblue","tomato"),c(47,30)))

boxplot(scale(log2(assay(sumexpcachexia))), col=rep(c("steelblue","tomato"),c(47,30)))

head(scale(assay(sumexpcachexia)), c(6,5))

prcompscale <- prcomp(scale(t(assay(sumexpcachexia))), center = F, scale. = F)
head(prcompscale$x, c(6,5))

prcomplog2 <- prcomp(log2(t(assay(sumexpcachexia))), center = F, scale. = F)
head(prcomplog2$x, c(6,5))

prcomplog2 <- prcomp(scale(log2(t(assay(sumexpcachexia)))))
head(prcomplog2$x, c(6,5))

head(prcomp(scale(t(assay(sumexpcachexia))))$x, c(6,5))

plot(prcompscale$x[,1], prcompscale$x[,2])

plot(prcomplog2$x[,1], prcomplog2$x[,2])
