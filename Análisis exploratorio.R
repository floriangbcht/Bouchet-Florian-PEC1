#######

#Cachexia is a complex metabolic syndrome associated with an underlying illness (such as cancer) 
#and characterized by loss of muscle with or without loss of fat mass (Evans et al., 2008). 
#A total of 77 urine samples were collected being 47 of them patients with cachexia, 
#and 30 control patients (from the "specmine.datasets" R package)


#Eisner et al. 2010 (Learning to predict cancer-associated skeletal muscle wasting from 1H-NMR profiles of urinary metabolites)
#H-NMR spectra acquisition (= proton nuclear magnetic resonance spectroscopy) and targeted profiling of metabolites from urine samples
#Data are concentrations of these metabolites


## Un data sanity/integrity check se ha hecho previamente:
#Successfully passed sanity check!"
#[2] "Samples are not paired."
#[3] "2 groups were detected in samples."
#[6] "All data values are numeric."
#[7] "A total of 0 (0%) missing values were detected."

#######

### Descarga y exploración de los datos brutos 

## Cargamos los datos y miramos los primeros registros y columnas
cachexiadata <- read.csv("https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2024-Cachexia/human_cachexia.csv")
head(cachexiadata, c(6,5))
# Tenemos los pacientes (=individuos) en líneas y los metabolitos (=variables) en columnas
# La primera columna representa los identificadores de los pacientes y la segunda el grupo al cual pertenecen 
# El resto de columnas son las concentraciones de los metabolitos

## Miramos cuantos grupos hay
unique(cachexiadata[,2])
# Hay dos, los que tienen la patología (cachexic) y los que no (control)

## Verificamos las dimensiones de la tabla
dim(cachexiadata)
# Hay 65 - 2 = 63 metabolitos (hay que deducir la columna con el id y la con el grupo) y 77 pacientes

## Se nota que no solo hay 2 grupos de pacientes pero hay 3 tipos de identificadores "PIF", "NETL" y "NETCR"
## Puede ser interesante tener un segundo factor de agrupamiento que tenga en cuenta el tipo de paciente
identif <- strsplit(cachexiadata[,1], "_")
identif <- unlist(lapply(identif, FUN = function(x)(x[1])))
unique(identif)



### Creación del objeto de clase SummarizedExperiment para los análisis siguientes

## Instalamos el paquete SummarizedExperiment y lo cargamos
if(!(require("SummarizedExperiment", quietly = T))){
  BiocManager::install("SummarizedExperiment")
}
library(SummarizedExperiment, quietly = T)

## Preparamos una matríz que solo contiene los valores de metabolitos: tiene que ser en el formato líneas = variables y columnas = individuos
valores <- t(as.matrix(data.frame(cachexiadata[,-c(1,2)], row.names = cachexiadata[,1])))
# Los nombres de líneas son los nombres de los metabolitos

## Preparamos una dataframe con los metadatos de las columnas, es decir la información sobre los pacientes (aquí los grupos)
columnasdatos <- data.frame(Group=cachexiadata[,2], row.names = cachexiadata[,1])
# Una sola columna y los nombres de líneas son los identificadores de los pacientes

# No hay metadatos para los metabolitos puesto que no tenemos información especifica sobre ellos

## Ahora podemos crear el objeto de clase SummarizedExperiment
sumexpcachexia <- SummarizedExperiment(assays = valores, colData = columnasdatos)

sumexpcachexia
rm(list=c("cachexiadata","valores", "columnasdatos"))



### Análisis exploratorio de los datos: análisis univariante y pre-procesamiento de los datos

## Vamos a utilisar funciones del paquete POMA para el pre-procesamiento
if(!(require("POMA", quietly = T))){
  BiocManager::install("POMA")
}
library(POMA, quietly = T)

## Primero, miramos si hay valores faltantes (= missing values = NAs)
if(mean(is.na(assay(sumexpcachexia)))==0){
  print("No hay missing values")
}
# No tenemos NAs en esos datos

## Luego hay también que gestionar los ceros que pueden ser convertidos en NAs gracias a PomaImpute
# Por defecto PomaImpute quita las líneas que tienen NAs y quita columnas enteras si esas columnas tienen más que 20% de NAs
# Se puede especificar que convierta los ceros en NAs
sumexpcachexiamod <- PomaImpute(sumexpcachexia, zeros_as_na = T)
# No hay missing values y entonces no se han quitado datos. Significa que tampoco habían ceros
rm(sumexpcachexiamod)

## Miramos la distribución de los datos, es decir la distribución de los valores de metabolitos para cada paciente, para ver si hay diferencias de magnitudes en los valores etc.
boxplot(assay(sumexpcachexia), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n", xlab="Patients", ylab="Metabolite values",
        main="Metabolite values for 77 patients in two groups")
axis(1, labels=F, at=c(1:77))
# Hay efectivamente diferencias enormes de valores para un mismo paciente y también de uno a otro

## Comparamos distintas opciones de scaling/transformación de los datos: 
# Estandarización (normalization) clásica vs. escalamiento log Pareto (por defecto en PomaNorm) vs. escalamiento log clásico vs. transformación log
# Como que los datos son muy asimétricos seguramente hará falta hacer una transformación log
par(mfrow=c(2,2))
for(mod_val in c("auto_scaling","log_pareto","log_scaling","log")){
  boxplot(assay(PomaNorm(sumexpcachexia, method = mod_val)), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n",
          main=mod_val)}
# Confirmamos, una transformación log (aquí log natural, pero dará los mismo que log2) parece adecuada

## Aplicamos la transformación
sumexpcachexia_log <- PomaNorm(sumexpcachexia, method = "log")

## Ahora miramos si después de la transformación log una estandárización sería útil
par(mfrow=c(1,2))
boxplot(assay(sumexpcachexia_log), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n",
        main="Log-transformed values")
boxplot(scale(assay(sumexpcachexia_log)), col=rep(c("steelblue","tomato"),c(47,30)), xaxt="n",
        main="Normalized log-transformed values")
par(mfrow=c(1,1))
# Claramente, hace falta

## Estandárizamos también los datos
sumexpcachexia_log_norm <- sumexpcachexia_log
assay(sumexpcachexia_log_norm) <- scale(assay(sumexpcachexia_log_norm))

rm(list=c("sumexpcachexia","sumexpcachexia_log"))

## Recordamos los tipos de pacientes y los grupos: podemos mirar la proporción de cada tipo de pacientes por grupo
table(identif, colData(sumexpcachexia_log_norm)[,1])
# Los tres tipos están representados en los dos grupos, será interesante para otros análisis



### Análisis exploratorio de los datos: Análisis multivariante

## Vamos a usar funciones del paquete mixOmics que permite, entre otros, hacer ACPs
if(!(require("mixOmics", quietly = T))){
  BiocManager::install("mixOmics")
}
library(mixOmics, quietly = T)

## Hacemos un ACP (seleccionamos de antemano un máximo de 20 PCs): recordar traspuestar la matríz para tener los individuos en líneas y variables en columnas y recordar NO hacer un scale porque ya lo hicimos (pero centrar los datos sí puesto que ahora la matriz es traspuesta) 
acp <- tune.pca(t(assay(sumexpcachexia_log_norm)), scale = F, ncomp = 20)

## Miramos la proporción de varianza explicada por cada PC
par(mfrow=c(1,2))
plot(acp, main="Proportion of explained variance for each PC")
plot(acp, main="Proportion of explained variance for each PC", type="b")
par(mfrow=c(1,1))
round(acp$prop_expl_var$X*100,2)
# El drop es continuo, no es muy marcado. La varianza se reparte en muchos PCs y por tanto los primeros PCs no explican mucha varianza...

## No vamos a quedarnos con docenas de PCs, nos quedarémos con 3 solos
acp <- tune.pca(t(assay(sumexpcachexia_log_norm)), scale = F, ncomp = 3)
round(acp$prop_expl_var$X*100,2)[1:3]
# Los tres primeros PCs explican aproximademente 11, 8, 7% de la varianza, respectivamente

## Miramos los gráficos con los scores, primero diferenciando los pacientes por grupo (cachexic vs. control)
plotIndiv(acp, ind.names = F, group = colData(sumexpcachexia_log_norm)[,1], ellipse = T, title = "PCA of 63 metabolites for two patient groups",
          legend = T, legend.title = "Group")
plotIndiv(acp, comp = c(2,3), ind.names = F, group = colData(sumexpcachexia_log_norm)[,1], ellipse = T, title = "PCA of 63 metabolites for two patient groups",
          legend = T, legend.title = "Group")
# Hay mucho solape claro, pero menos en el plot PC1 vs. PC2 que en el plot PC2 vs. PC3. El PC3 no sirve mucho para separar los dos grupos

## Las elipses pueden ser un poco confundientes, entonces vamos a mirar el plot PC1 vs. PC2 de otra manera, usando convex hulls
library(ade4)
plot(acp$variates$X[,1:2], pch=19, col = rep(c("steelblue","tomato"),c(47,30)))
s.chull(acp$variates$X[,1:2], fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("bottomright", legend=levels(factor(colData(sumexpcachexia_log_norm)[,1])), pch=19, col = c("steelblue","tomato"),title="Grupo")
# Más comprensible así. Pero confirmamos que hay mucho solape

## Ahora miramos si pasa algo al mirar el ACP por tipo de pacientes (PIF vs. NETL vs. NETCR)
plotIndiv(acp, ind.names = F, group = identif, ellipse = T, title = "PCA of 63 metabolites for three patient types",
          legend = T, legend.title = "Type")
plot(acp$variates$X[,1:2], pch=19)
s.chull(acp$variates$X[,1:2], fac = factor(identif), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato","darkgrey"), add.plot = T)
legend("topleft", legend=levels(factor(identif)), pch=19, col = c("steelblue","tomato","darkgrey"),title="Type")
# Hay una distinción entre los tres tipos, de hecho se nota que tenemos básicamente los PIF vs. los NETL+NETCR
# Pero recordamos que cada tipo de paciente está representado en cada grupo (cachexic o control), o sea que la distinción por tipo de pacientes no refleja la distinción por grupo
# Tenemos un efecto batch debido a si el paciente es PIF o si es NETL/NETCR!

## Miramos como las variables se correlan entre ellas, y con los PCs, mirando la separación entre pacientes por tipo puesto que por grupo hay mucho solape (es decir, el efecto batch parece más interesante de estudiar)
biplot(acp, group = identif, ind.names = F, var.names.size = 3, legend.title = "Type")
# Hay muchos metabolitos por lo que es dificil comentar las correlaciones entre ellos
# Ningúna variable tiene correlación estricta con uno o otro PC: Es decir, ningún vector propio es bien parallelo a uno o otro PC
# Aunque no hayan agrupamientos claros de los metabolitos podemos decir que los más a la izquierda pueden contribuir a tener scores más negativos del PC1 y los más a la dercha a tener scores más positivos del PC1
# Básicamente eso metabolitos más a la izquierda y más a la derecha pueden contribuir a la separación más marcada entre los pacientes PIF y los pacientes NETL+NETCR

## Miramos directamente los valores de los vectores propios, mediante un gráfico
par(mfrow=c(1,2))
plotLoadings(acp, comp = 1, size.title = 1)
plotLoadings(acp, comp = 2, size.title = 1)
# No hay variables (metabolitos) que se destacan especificamente más que las demás
# Quizá valdría la pena seleccionar un subconjunto de variables y ver si los vectores propios cambian (y por tanto, últimamente, la variabilidad explicada por cada PC)
# Sería coma una filtración de las variables, que al final puede formar parte del procesamiento de los datos...

## Testamos si hay un número óptimo de variables a incluir que podría maximizar la variabilidad explicada por los PCs (aquí 3 PCs)
# Testamos con todas las variables para 3 PCs, y repetimos 25 veces la validación cruzada (máximo 50 pero ya con 30 deberíamos obtener un resultado coherente)
# Básciamente esta función tune.spca permite obtener el número de variables a incluir que maximiza la correlación entre ellas y los PCs, para cada PC, empezado con el número más bajo posible de variables
acp_mod <- tune.spca(t(assay(sumexpcachexia_log_norm)), scale = F, ncomp = 3, 
                              folds = 5, 
                              test.keepX = 1:nrow(assay(sumexpcachexia_log_norm)), nrepeat = 30) 
acp_mod$choice.keepX
# Puesto que el resultado es un número alto de variables, cerca del número original de variables que tenemos, no vale la pena cambiar nada


## Además de un ACP podemos hacer análisis basados en distancias (que permiten también reducir las dimensiones y visualizar las relaciones entre individuos)
# Creamos dos sets de distancias (distancias entre pacientes), uno con distancias euclidianas y otro con distancias basadas en correlaciones de Pearson
# Distancias basadas en la correlación se usan bastante en genomica (datos de expresión), y quizá puedan aportar algo aquí
dist_euc <- dist(t(assay(sumexpcachexia_log_norm)))
dist_cor <- as.dist(1-cor(assay(sumexpcachexia_log_norm)))

## Primero, haremos un non-metric multidimensional scaling, usando el mapping no lineal de Sammon, para contrastar el ACP (que reduce las dimensiones de manera lineal)
library(MASS, quietly = T)
# Con distancias euclidianas
samm <- sammon(dist_euc, trace=FALSE)
par(mfrow=c(1,2))
plot(samm$points, pch = rep(c(19,17),c(47,30)), col = rep(c("steelblue","tomato"),c(47,30)), xlab='', ylab='',
     main="Sammon's NLM of 77 patients from metabolite values", sub="Euclidean distances", cex.sub=1.5)
s.chull(samm$points, fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("bottomleft", legend=levels(factor(colData(sumexpcachexia_log_norm)[,1])), pch=c(19,17), col = c("steelblue","tomato"),title="Group")
# Notamos el solape marcado entre los dos grupos. Hay más variabiliad en el grupo cachexic
plot(samm$points, pch=19, , xlab='', ylab='', main="Sammon's NLM of 77 patients from metabolite values", sub="Pearson distances", cex.sub=1.5)
s.chull(samm$points, fac = factor(identif), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato","darkgrey"), add.plot = T)
legend("bottomleft", legend=levels(factor(identif)), pch=19, col = c("steelblue","tomato","darkgrey"),title="Type")
# Tenemos otra vez los pacientes PIF por un lado y los pacientes NETCR y NETL por otro lado (efecto batch). Hay más variabiliad en los pacientes PIF
# El paciente con coordenadas aprox (-2, -10) podría representar un outlier...

# Con distancias basadas en correlaciones de Pearson
samm2 <- sammon(dist_cor, trace=FALSE)
plot(samm2$points, pch = rep(c(19,17),c(47,30)), col = rep(c("steelblue","tomato"),c(47,30)), xlab='', ylab='',
     main="Sammon's NLM 77 patients from metabolite values", sub="Euclidean distances", cex.sub=1.5)
s.chull(samm2$points, fac = factor(colData(sumexpcachexia_log_norm)[,1]), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato"), add.plot = T)
legend("bottomleft", legend=levels(factor(colData(sumexpcachexia_log_norm)[,1])), pch=c(19,17), col = c("steelblue","tomato"),title="Group")
# El solape aún es marcado pero la variabilidad es más parecida ente los dos grupos
plot(samm2$points, pch=19, xlab='', ylab='', main="Sammon's NLM 77 patients from metabolite values", sub="Pearson distances", cex.sub=1.5)
s.chull(samm2$points, fac = factor(identif), xax=1, optchull = 1, clabel = 0, col = c("steelblue","tomato","darkgrey"), add.plot = T)
legend("bottomleft", legend=levels(factor(identif)), pch=19, col = c("steelblue","tomato","darkgrey"),title="Type")
# Pasa lo mismo con los tres tipos de pacientes. Aún tenemos los pacientes PIF separados de los NETCR+NETL (efecto batch).
# El paciente con coordenadas aprox (-0.2, 0.20) podría representar el mismo outlier que en el análisis basado en distancias euclidianas

## Podemos también hacer un análisis de conglomerados (método UPGMA) usando esas distancias y especificando para cada paciente su grupo y su tipo
doublegroup <- paste(identif, rep(c("cachexic","control"),c(47,30)))
plot(hclust(dist_euc,method="average"), labels = doublegroup, hang=-1, xlab='', cex=0.75,
     main="UPGMA cluster of 77 patients from metabolite values", sub="Euclidean distances", cex.sub=1.5)
plot(hclust(dist_cor,method="average"), labels = doublegroup, hang=-1, xlab='', cex=0.75,
     main="UPGMA cluster of 77 patients from metabolite values", sub="Pearson distances", cex.sub=1.5)
# Las topologías de los dos clusteres son bastante similares. En general, los individuos de un mismo tipo (ej., PIF) y de un mismo grupo (ej., cachexic) se agrupan entre ellos y también con otros individuos del mismo grupo
# Pero se nota también que individuos de un mismo tipo (ej., NETL) pueden agruparse aunque no sean del mismo grupo, en clusteres más grandes
# Esto confirma que hay una cierta diferenciación entre los dos grupos (cachexic vs. control) pero que también hay cierto grado de diferenciación ente los tipos de pacientes y que por tanto oculta las diferencias entre los dos grupos (efecto batch)

## Finalmente, podemos hacer un análisis de correspondencia para ir más allá y relacionar los pacientes con los metabolitos
# Usamos la función ord para el análisis. A la función le roporcionamos directamente la matriz del objeto SummarizedExperiment sin transponerla
library(made4, quietly = T)
ca <- ord(assay(sumexpcachexia_log_norm), type="coa")
round(ca$ord$eig^2 / sum(ca$ord$eig^2) *100,2)
# Los primeros cPCs explican más varianza que los primeros PCs del ACP que hicimos más arriba (PC1 15% más y PC2 7% más) pero tampoco explican mucho

# Miramos los resultados pintando los pacientes por grupo (cachexic vs. control)
plot.ord(ca, classvec = colData(sumexpcachexia_log_norm)[,1])
# Miramos también los resultados pintando los pacientes por tipo (PIF vs. NETL vs. NETCR)
plot.ord(ca, classvec = identif)
# En ambos casos no hay agrupamientos claros respecto a las variables (los metabolitos)
# Cuando pintamos por grupo vemos que hay aún más solapamiento que en el ACP y el mapping no lineal previos
# Cuando pintamos por tipo de nuevo vemos el efecto batch visto previamente. 
# En este segundo caso podemos ver un cierto corte entre los metabolitos aunque no este muy claro: Un grupo a la izquierda (Sucrose, Acetate, Succinate etc.), y uno a la derecha (Acetone, Tartrate, Pyruvate etc.)
# Podemos decir que máss o menos el grupo de la izquierda corresponde a los pacientes PIF y el de la derecha corresponde a los pacientes NETL+NETCR

## Para acabar podemos intentar ver si existen realmente grupos (clusteres) de metabolitos en nuestros datos usando K-means
# Primero miramos cual es el número óptimo de clusteres usando el método del codo (para ver con cual número de clusteres la variabilidad intragrupos es miníma)
within_ss <- c()
for (i in 1:10) {
  within_ss[i] <- sum(kmeans(assay(sumexpcachexia_log_norm), i, iter.max = 1000)$withinss)
}
par(mfrow=c(1,1))
plot(1:10,within_ss, type="b", xlab="Number of clusters",
     ylab="Within-group sum of squares", main="Within-group variability depending on cluster number")
# Parece que es 4. Si repetimos el proceso varias veces siempre saldrá que con K = 4 tenemos el valor de variabilidad intragrupos miníma (o sea, a partir de este K este valor ya no baja significativamente)

# Aplicamos la función
k_means <- kmeans(assay(sumexpcachexia_log_norm), 4, iter.max = 1000)

# Hacemos un gráfico con los clusters usando el paquete factoextra
library(factoextra, quietly = T)
fviz_cluster(k_means, data = assay(sumexpcachexia_log_norm), palette = c("steelblue", "tomato", "forestgreen", "brown"), 
             ellipse.type = "euclid", star.plot = TRUE, main="Clusters of metabolites from k-means analysis")
# Efectivamente tenemos 4 clusteres con poco solape aunque hay cercanía, el último siendo muy pequeño y más separado de los otros
# Podría ser interesante a partir de este resultado intentar ver como el grupo y el tipo de pacientes explica esta variabilidad (hacerlo al reves que en el análisis de correspondencia previo)



