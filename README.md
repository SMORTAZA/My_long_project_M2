#README
#Author : SMORTAZA
#Date : date du rapport

##Intégration de données transcriptomiques et protéomiques dans un modèle d'hépatocytes humains

```{r}
#source("http://bioconductor.org/biocLite.R")
#biocLite("oligo")
#source("https://bioconductor.org/biocLite.R")
#biocLite("pd.hugene.2.0.st")
#source("https://bioconductor.org/biocLite.R")
#biocLite("hugene20sttranscriptcluster.db")
```

Lecture des fichiers.CEL pour obtention d'un table data regroupant toutes les données

```{r}
#verification du repertoire courant
getwd()
#changement de repertoire
setwd("./data")
#verification du changement de repertoire
getwd()
#chargement de la librairie pour analyse array
library(oligo)
#liste des fichiers .CEL
celFiles <- list.celfiles()
#Lecture des fichiers .CEL
affyRaw <- read.celfiles(celFiles)
#chargement de la librairie du array utilisé
library(pd.hugene.2.0.st)
#RMA sur cet array pour normalisation des données
eset <- rma(affyRaw)
#ecriture des resultats dans un fichier data
write.exprs(eset,file="./../results/data.txt")
#Sortir des dossiers contenant les données (fichiers .CEL)
setwd("./..")
#Verifier ce changement de repertoire
getwd()
#enregistrement des data dans un objet data
data <- read.table('./results/data.txt')
#verification de l'objet data
dim(data)
#Visualisation de la table
View(data)
```

Modification des noms de colonnes

```{r}
#nom des colonnes
colnames(data)
ech_names <- c('Baseline', 'Ctrl', 'HNO3', 'M1', 'M2')
n = 1 #n° de replicat
ind = 0 #n° d'echantillon ds le replicat
for (i in seq(length(colnames(data)))){
  col_name <- colnames(data)[i]
  ind = ind + 1
  colnames(data)[i] <- paste(ech_names[ind], '_', n, sep="")
  if (ind == 5){
   n = n + 1
   ind = 0
  }
}
#Vérification du nom des colonnes
colnames(data)
#Enregistrement de la matrice de résultats avec les bons noms de colonnes dans un RData
save(data, file = "./results/data.RData")
```

Modification des noms de lignes (remplacement par des noms de genes)

```{r}
load('./results/data.RData')
library(hugene20sttranscriptcluster.db)
#hugene20sttranscriptcluster()
#permet d'avoir un dataframe Annot contenant les annotations
Annot <- data.frame(ACCNUM=sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse=", "))
#Enregistrement des résultats d'annotation dans un Rdata
save(data, file = "./results/Annot.RData")
#Creation d'une colonne contenant également le nom des probes
data$probe_id <- rownames(data)
#Enregistrement du tableau modifié
save(data, file = "./results/data.RData")
```

Loading and normalizing

Now, I have a matrix containing the corrected and summarized data at probe-set level:
5 baseline - columns 1, 6, 11, 16, 21
5 controle - columns 2, 7, 12, 17, 22
5 HNO3 - columns 3, 8, 13, 18, 23
5 M1 - columns 4, 9, 14, 19, 24
5 M2 - columns 5, 10, 15, 20, 25
This is also indicated by the column names of the matrix 'data'


```{r}
#Chargement de la matrice
load('./results/data.RData')
#Verification de la normalisation de la matrice avec density plot
plot(density(data[,1], na.rm = TRUE), main = "Normalized data", xlab = "Intensity")
for (i in seq(2,ncol(data)-2)){
  points(density(data[,i], na.rm = TRUE), type = "l", col = rainbow(20)[i])
}
```

Singular Value Decomposition and outlier assessment

In order to evaluate the quality of microarray data there are many things that can be done. One such thing is the density plot, that I've just made. 
Another good way of evaluating my data is to do a Singular Value Decomposition (SVD), which is kind of the same as a Principal Component Analysis (PCA), just designed to take matrice that have many more rows than columns. An SVD analysis kind of breaks my data into 25 components (= number of samples). The first component contains the most variation, and the last the least (none). 

```{r}
real_data = data[,c(-26,-27)]
m <- ncol(real_data)
X <- real_data - rowMeans(real_data)
svd <- svd(X)
V <- svd$v
S <- matrix(0,m,m)
#pourcentage de variance expliquée par chaque composante principale
diag(S) <- svd$d
#scree plot de la variance expliquée par chaque composante principale
plot(svd$d, type = 'b', col = 'blue', lwd = 3, main = 'Explained Variance of Principal Component', xlab = 'Principal Component', ylab = 'Explained Variance (in %)')
#nom des echantillons (conditions)
pnames <- colnames(real_data)
#couleur selon echantillons (conditions de experience)
cat_col <- rep(c("blue", "darkgreen", "red", "orange", "black"), 5)
#par(mfrow = c(1,2))
#heatmap des 25 composantes principales
image( S%*%t(V), main="Singular Value Decomposition",
ylab = "Samples ( 1-25 )", xlab = "The 25 Components" )
#SVD
plot(V[,2], V[,3], col = "white", main = "Singular Value
Decomposition",
xlab = paste("Second Component (",round(svd$d[2]),"%)"), ylab = paste("Third Component (",round(svd$d[3]),"%)"))
text(V[,2], V[,3], pnames, col = cat_col)
```

The Two-Way ANOVA [entre conditions X et Y]

The two-way ANOVA uses all the data in a single test. This enables us to not only find the genes that are expressed differently in X compared to Y, but also the genes that are expressed differently between the two cell types. On top of this we will also get a list of genes that exhibit an interaction effect between the two factors. 