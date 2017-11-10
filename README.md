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
save(data, file = "./results/data.RData")
```

```{r}
load('./results/data.RData')

```

```{r}
#lignes de code à ne pas prendre en compte (brouillon)
line_num <- rownames(data)
data$Ensembl_gene_ID <- rep(NA)
a <- names(xx)
for (i in seq(length(line_num))){
  if (which(a == line_num[i]) == integer(0)){
    print('ok')
  }
  #data$Ensembl_gene_ID[i] <- xx[which(a == line_num[i])]
  #print(line_num[i])
}

## Bimap interface:
library(hugene20sttranscriptcluster.db)
x <- hugene20sttranscriptclusterENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
# Get the Ensembl gene IDs for the first five genes
xx[1:5]
# Get the first one
xx[[1]]
}
```


