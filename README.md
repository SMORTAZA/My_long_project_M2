#README
#Author : SMORTAZA
#Date : date du rapport

##Intégration de données transcriptomiques et protéomiques dans un modèle d'hépatocytes humains

```{r}
#source("http://bioconductor.org/biocLite.R")
#biocLite("oligo")
#source("https://bioconductor.org/biocLite.R")
#biocLite("pd.hugene.2.0.st")
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
write.exprs(eset,file="data.txt")
#enregistrement des data dans un objet data
data <- read.table('data.txt')
#verification de l'objet data
dim(data)
#Visualisation de la table
View(data)
#Sortir des dossiers contenant les données (fichiers .CEL)
setwd("./..")
#Verifier ce changement de repertoire
getwd()
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
colnames(data)
```

Modification des noms de lignes (remplacement par des noms de genes)

```{r}
```
