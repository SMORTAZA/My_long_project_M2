#####README

#####Date : 08/01/2018

####Intégration de données transcriptomiques et protéomiques dans un modèle d'hépatocytes humains

#### Indications 

- Lancer les lignes de commandes dans R (ou RStudio) en se plaçant au niv du répertoire My_project
(Plus tard, toutes les lignes de commande seront mis dans un seul script à part pour pouvoir lancer ensuite une ligne de commande - répertoire src/ pas à utiliser pour l'instant)
- Tous les résultats obtenus (figures ...) sont dans le répertoire results/

### Importation et Normalisation des données

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

### Visualisation des données

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

### Ttests (calcul des pvalues et fold change)

Prendre uniquement les colonnes qui nous intéressent pour les tests statistiques

```{r}
#rearrangement des col utiles de data
w_data=data[,c(1,6,11,16,21,2,7,12,17,22,3,8,13,18,23,4,9,
               14,19,24,5,10,15,20,25)]
#Enregistrement du nouveau tableau
save(w_data, file = "./results/w_data.RData")
#ouverture de ces données
load('./results/data.RData')
load('./results/w_data.RData')
```
Création de design pour chaque comparaison de conditions:
- 1°) Contrôle vs HNO3
- 2°) HNO3 vs M1
- 3°) HNO3 vs M2
- 4°) M1 vs M2
- 5°) Baseline vs Contrôle

```{r}
#creation du design
design <- matrix(c(rep(c(1,0,0,0,0),5),rep(c(0,1,0,0,0),5),
                   rep(c(0,0,1,0,0),5),rep(c(0,0,0,1,0),5),
                   rep(c(0,0,0,0,1),5)), ncol=5,byrow=TRUE)
colnames(design) <- c("Baseline","Ctrl","HNO3","M1","M2")
rownames(design) <- colnames(w_data)
design <- data.frame(design)
```

(Volcanoplot pr determiner un seuil de significativité pour pvalue) 
A volcano plot shows the connection between the pvalues and the log2 of the fold change, compared to the same to the same analysis of the same data randomized. the permutation should be a balanced randomization of the colums of the experiment. After the permutation we run the ttest again and calculate a new log2 fold change (M) value. Finally, we plot the two outputs on top of each other in order to compare.
So, in order to make a volcanoplot, we need to obtain the random or permutated pvalues for each dataset. To get those, we need to perform a ttest on the same data, but with permutated labels. Similarly for the "permutated" fold-change. Finally, we'll plot the "permutated" values vs permutated fold-change (red dots) on top of the "real" data.

Création du "design permuté", correspondant aux données simulées

```{r}
#avoir le nb de lignes et de colonnes de la matrice design
dim(design)
#creation de la matrice permutated de la meme dimension que la matrice design
#permutated_design = matrix(nrow = 25, ncol = 5)
#for (i in seq(dim(design)[2])){
#  print(i)
#  print(design[,i])
#  random_lables <- sample(design[,i])
#  print(random_lables)
#  permutated_design[,i] <- random_lables
#}
#colnames(permutated_design) <- c("Baseline","Ctrl","HNO3","M1","M2")
#rownames(permutated_design) <- colnames(w_data)
#permutated_design <- data.frame(permutated_design)
#enregistrement des données car comme sample, alors change tt le tps
#save(permutated_design, file = "./results/permutated_design.RData")
#chargement des données
load('./results/permutated_design.RData')
#obtention des col permutées fom w_data
Baseline_permutated <- w_data[,c(permutated_design$Baseline==1)]
Ctrl_permutated <- w_data[,c(permutated_design$Ctrl==1)]
HNO3_permutated <- w_data[,c(permutated_design$HNO3==1)]
M1_permutated <- w_data[,c(permutated_design$M1==1)]
M2_permutated <- w_data[,c(permutated_design$M2==1)]
```

Obtention des pvalues avec un ttest.

```{r}
require(genefilter)
#en chargeant w_data, son typeof est list o
#typeof(data)
#modifier le typeof en double
w_data <- data.matrix(w_data)
#verification de la modif
#typeof(data)
#calcul des pvalues
#1°) Données réelles
t.pval.CtrlvsHNO3 <- rowttests(w_data[,seq(6,15)],
                               factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
t.pval.HNO3vsM1 <- rowttests(w_data[,seq(11,20)],
                             factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
t.pval.HNO3vsM2 <- rowttests(w_data[,c(11,12,13,14,15,21,22,23,24,25)],
                             factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
t.pval.M1vsM2 <- rowttests(w_data[,c(16,17,18,19,20,21,22,23,24,25)],
                           factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
t.pval.BaselinevsCtrl <- rowttests(w_data[,seq(1,10)],
                           factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
#2°) Données simulées
#permutated pval
perm_t.pval.CtrlvsHNO3 <- rowttests(cbind(Ctrl_permutated,HNO3_permutated), 
                                    factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
perm_t.pval.HNO3vsM1 <- rowttests(cbind(HNO3_permutated,M1_permutated), 
                                  factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
perm_t.pval.HNO3vsM2 <- rowttests(cbind(HNO3_permutated,M2_permutated), 
                                  factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
perm_t.pval.M1vsM2 <- rowttests(cbind(M1_permutated,M2_permutated), 
                                factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
perm_t.pval.BaselinevsCtrl <- rowttests(cbind(Baseline_permutated,Ctrl_permutated), 
                                factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
```

Obtention des fold change (fc)

```{r}
#to know how genes are affected (up or down), 
#calculate the log2 of the fold-change
#1°) Données réelles
fc.CtrlvsHNO3 <- rowMeans(w_data[,design$Ctrl==1])-
  rowMeans(w_data[,design$HNO3==1])
fc.HNO3vsM1 <- rowMeans(w_data[,design$HNO3==1])-
  rowMeans(w_data[,design$M1==1])
fc.HNO3vsM2 <- rowMeans(w_data[,design$HNO3==1])-
  rowMeans(w_data[,design$M2==1])
fc.M1vsM2 <- rowMeans(w_data[,design$M1==1])-
  rowMeans(w_data[,design$M2==1])
#2°) Données simulées
perm_fc.CtrlvsHNO3 <- rowMeans(w_data[,permutated_design$Ctrl==1])-
  rowMeans(w_data[,permutated_design$HNO3==1])
perm_fc.HNO3vsM1 <- rowMeans(w_data[,permutated_design$HNO3==1])-
  rowMeans(w_data[,permutated_design$M1==1])
perm_fc.HNO3vsM2 <- rowMeans(w_data[,permutated_design$HNO3==1])-
  rowMeans(w_data[,permutated_design$M2==1])
perm_fc.M1vsM2 <- rowMeans(w_data[,permutated_design$M1==1])-
  rowMeans(w_data[,permutated_design$M2==1])
perm_fc.BaselinevsCtrl <- rowMeans(w_data[,permutated_design$Baseline==1])-
  rowMeans(w_data[,permutated_design$Ctrl==1])
```

### Volcanoplots

Pour choisir le seuil, il faut que sous ce cut-off, il y ait plus de points bleus (données réelles) que de points rouges (données simulées) car c'est là où les gènes vont être le plus significatif.

```{r}

plot(fc.CtrlvsHNO3, t.pval.CtrlvsHNO3, main = "Volcano Plot\nCtrl vs HNO3", 
     log = "y",xlab = "M(log2 fold change)", ylab = "p-value", pch = 20, 
     col = "blue")
points(perm_fc.CtrlvsHNO3, perm_t.pval.CtrlvsHNO3, type = "p", pch = 20, 
       col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"),pch=20,
       bg="white")
#seuil t1
t1 <- 0.001
#verification du seuil
length(which(t.pval.CtrlvsHNO3 < t1)) #28 pvalues
length(which(perm_t.pval.CtrlvsHNO3 < t1)) #26 pvalues
abline(h= t1)
#plot du volcano de la seconde comparaison
plot(fc.HNO3vsM1, t.pval.HNO3vsM1, main = "Volcano Plot\nHNO3 vs M1", 
     log = "y",xlab = "M(log2 fold change)", ylab = "p-value", pch = 20,
     col = "blue")
points(perm_fc.HNO3vsM1, perm_t.pval.HNO3vsM1, type = "p", pch = 20,
       col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"),pch=20,
       bg="white")
#seuil t2
t2 <- 0.004
#verification du seuil --> deja ici bcp trop de genes seront selectionnés
length(which(t.pval.HNO3vsM1 < t2)) # 290 pvalues
length(which(perm_t.pval.HNO3vsM1 < t2)) # 21 pvalues
abline(h=t2)
#plot du volcano de la troisieme comparaison
plot(fc.HNO3vsM2, t.pval.HNO3vsM2, main = "Volcano Plot\nHNO3 vs M2", 
     log = "y",xlab = "M(log2 fold change)", ylab = "p-value", pch = 20, 
     col = "blue")
points(perm_fc.HNO3vsM2, perm_t.pval.HNO3vsM2, type = "p", pch = 20, 
       col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"),pch=20,
       bg="white")
#seuil t3
t3 <- 0.004
#verification du seuil --> bcp de genes seront selectionnés
length(which(t.pval.HNO3vsM2 < t3)) #300 pvalues
length(which(perm_t.pval.HNO3vsM2 < t3)) #157 pvalues
abline(h=t3)
#plot du volcano de la quatrieme comparaison
plot(fc.M1vsM2, t.pval.M1vsM2, main = "Volcano Plot\nM1 vs M2", 
     log = "y", xlab = "M(log2 fold change)", ylab = "p-value", 
     pch = 20, col = "blue")
points(perm_fc.M1vsM2, perm_t.pval.M1vsM2, type = "p", pch = 20, 
       col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"),pch=20,
       bg="white")
#seuil t4
t4 <- 0.02
#verification du seuil --> bcp de genes seront selectionnés
length(which(t.pval.M1vsM2 < t4)) # 668 pvalues
length(which(perm_t.pval.M1vsM2 < t4)) # 29 pvalues
abline(h=t4)
#plot du volcano de la cinquieme comparaison
plot(fc.BaselinevsCtrl, t.pval.BaselinevsCtrl, main = "Volcano Plot\n
     Baseline vs Ctrl", log = "y",xlab = "M(log2 fold change)", 
     ylab = "p-value", pch = 20, col = "blue")
points(perm_fc.BaselinevsCtrl, perm_t.pval.BaselinevsCtrl, type = "p", 
       pch = 20, col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"),pch=20,
       bg="white")
#seuil t5
t5 <- 0.02
#verification du seuil --> bcp de genes seront selectionnés
length(which(t.pval.BaselinevsCtrl < t5)) # 2315 pvalues
length(which(perm_t.pval.BaselinevsCtrl < t5)) # 18 pvalues
abline(h=t5)
```

### ANOVA

The Two-Way ANOVA [entre conditions X et Y]

The two-way ANOVA uses all the data in a single test. This enables us to not only find the genes that are expressed differently in X compared to Y, but also the genes that are expressed differently between the two cell types. On top of this we will also get a list of genes that exhibit an interaction effect between the two factors. 
Cet anova est effectué avec le package limma, utilisé justement pour ce type de données. 

```{r}
#importation de la librairie limma
library(limma)
#conditions à comparer
contrast.matrix <- makeContrasts(Ctrl-HNO3, HNO3-M1, HNO3-M2, M1-M2, 
                                 Baseline-Ctrl, levels=design)
#paramètres du modèle ajustés en fct de nos données
fit <- lmFit(log2(w_data), design)
#recherche de differences significatives
eBayesResultat <- eBayes(contrasts.fit(fit, contrast.matrix))
#pvalues
CtrlvsHNO3 <- eBayesResultat$p.value[,1]
HNO3vsM1 <- eBayesResultat$p.value[,2]
HNO3vsM2 <- eBayesResultat$p.value[,3]
M1vsM2 <- eBayesResultat$p.value[,4]
BaselinevsCtrl <- eBayesResultat$p.value[,5]
#rechercher les genes significatifs
top1 <- which(sort(CtrlvsHNO3) <= t1)
top2 <- which(HNO3vsM1 <= t2) 
top3 <- which(HNO3vsM2 <= t3)
top4 <- which(M1vsM2 <= t4)
top5 <- which(BaselinevsCtrl <= t5)
#connaitre le nombre de genes significatifs
length(top1)
length(top2)
length(top3)
length(top4)
length(top5)
```

2 types de Heatmap pour chacune des cinq conditions avec les genes significatifs sélectionnés:
- 1°) avec les probes
- 2°) avec les noms des gènes (par contre, beaucoup de NA, donc pas informatif)

```{r}
#heatmap1 avec les genes significatifs
#avec les probes
matrix <- as.matrix(top1)
new <- as.matrix(w_data[matrix[,1],seq(6,15)])
heatmap(new, main = "Ctrl vs HNO3", xlab = "samples             ", 
        ylab="significant genes")
#avec les genes
new_genes <- as.matrix(w_data[matrix[,1],seq(6,15)])
rownames(new_genes) <- data$geneSymbol[matrix[,1]]
heatmap(new_genes, main = "Ctrl vs HNO3")
#heatmap2 avec les genes significatifs
#avec les probes
matrix <- as.matrix(top2)
new <- as.matrix(w_data[matrix[,1],seq(11,20)])
heatmap(new, main = "HNO3 vs M1", xlab = "samples           ", 
        ylab="significant genes")
#avec les genes
new_genes <- as.matrix(w_data[matrix[,1],seq(11,20)])
rownames(new_genes) <- data$geneSymbol[matrix[,1]]
heatmap(new_genes, main = "HNO3 vs M1")
#heatmap3 avec les genes significatifs
#avec les probes
matrix <- as.matrix(top3)
new <- as.matrix(w_data[matrix[,1],c(11,12,13,14,15,21,22,23,24,25)])
heatmap(new, main = "HNO3 vs M2", xlab = "samples            ", 
        ylab="significant genes")
#avec les genes
new_genes <- as.matrix(w_data[matrix[,1],c(11,12,13,14,15,21,22,23,24,25)])
rownames(new_genes) <- data$geneSymbol[matrix[,1]]
heatmap(new_genes, main = "HNO3 vs M2")
#heatmap4 avec les genes significatifs
#avec les probes
matrix <- as.matrix(top4)
new <- as.matrix(w_data[matrix[,1],c(16,17,18,19,20,21,22,23,24,25)])
heatmap(new, main = "M1 vs M2", xlab = "samples", ylab="significant genes")
#avec les genes
new_genes <- as.matrix(w_data[matrix[,1],c(16,17,18,19,20,21,22,23,24,25)])
rownames(new_genes) <- data$geneSymbol[matrix[,1]]
heatmap(new_genes, main = "M1 vs M2")
#heatmap5 avec les genes significatifs
#avec les probes
matrix <- as.matrix(top5)
new <- as.matrix(w_data[matrix[,1],seq(1,10)])
heatmap(new, main = "Baseline vs Ctrl", xlab = "samples             ", 
        ylab="significant genes")
#avec les genes
new_genes <- as.matrix(w_data[matrix[,1],seq(1,10)])
rownames(new_genes) <- data$geneSymbol[matrix[,1]]
heatmap(new_genes, main = "Baseline vs Ctrl")
```

### Enrichissement Biologique avec Pathways et GOterms

Pour chaque condition, les probes et noms de gènes sont enregistrés dans des fichiers textes pour pouvoir les utiliser par la suite

```{r}
#noms des genes pour la premiere comparaison
n = 0
gene_names_1 <- vector(length = length(top1_sondes))
for (i in top1_sondes){
  n = n + 1 
  #print(i)
  gene_symbol <- data$geneSymbol
  gene_symbol <- as.matrix(gene_symbol)
  gene_name <- gene_symbol[which(data$probe_id == i)]
  gene_names_1[n] <- gene_name
  #enregistrer ds un fichier le nom des genes 
  capture.output(gene_name, file="./results/Tests_statistiques/Significant_genes/CtrlvsHNO3.txt", append = TRUE)
  #enregistrer les probes
  capture.output(i, file="./results/Tests_statistiques/Significant_probes/probes_CtrlvsHNO3.txt", append = TRUE)
}
#noms des genes pour la deuxième comparaison
n = 0
gene_names_2 <- vector(length = length(top2_sondes))
for (i in top2_sondes){
  n = n + 1 
  #print(i)
  gene_symbol <- data$geneSymbol
  gene_symbol <- as.matrix(gene_symbol)
  gene_name <- gene_symbol[which(data$probe_id == i)]
  gene_names_2[n] <- gene_name
  #enregistrer ds un fichier le nom des genes 
  capture.output(gene_name, file="./results/Tests_statistiques/Significant_genes/HNO3vsM1.txt", append = TRUE)
  #enregistrer les probes
  capture.output(i, file="./results/Tests_statistiques/Significant_probes/probes_HNO3vsM1.txt", append = TRUE)
}
#noms des genes pour la troisième comparaison
n = 0
gene_names_3 <- vector(length = length(top3_sondes))
for (i in top3_sondes){
  n = n + 1 
  #print(i)
  gene_symbol <- data$geneSymbol
  gene_symbol <- as.matrix(gene_symbol)
  gene_name <- gene_symbol[which(data$probe_id == i)]
  gene_names_3[n] <- gene_name
  #enregistrer ds un fichier le nom des genes 
  capture.output(gene_name, file="./results/Tests_statistiques/Significant_genes/HNO3vsM2.txt", append = TRUE)
  #enregistrer les probes
  capture.output(i, file="./results/Tests_statistiques/Significant_probes/probes_HNO3vsM2.txt", append = TRUE)
}
#noms des genes pour la quatrieme comparaison
n = 0
gene_names_4 <- vector(length = length(top4_sondes))
for (i in top4_sondes){
  n = n + 1 
  #print(i)
  gene_symbol <- data$geneSymbol
  gene_symbol <- as.matrix(gene_symbol)
  gene_name <- gene_symbol[which(data$probe_id == i)]
  gene_names_4[n] <- gene_name
  #enregistrer ds un fichier le nom des genes 
  capture.output(gene_name, file="./results/Tests_statistiques/Significant_genes/M1vsM2.txt", append = TRUE)
  #enregistrer les probes
  capture.output(i, file="./results/Tests_statistiques/Significant_probes/probes_M1vsM2.txt", append = TRUE)
}
#noms des genes pour la cinquieme comparaison
n = 0
gene_names_5 <- vector(length = length(top5_sondes))
for (i in top5_sondes){
  n = n + 1 
  #print(i)
  gene_symbol <- data$geneSymbol
  gene_symbol <- as.matrix(gene_symbol)
  gene_name <- gene_symbol[which(data$probe_id == i)]
  gene_names_5[n] <- gene_name
  #enregistrer ds un fichier le nom des genes 
  capture.output(gene_name, file="./results/Tests_statistiques/Significant_genes/BaselinevsCtrl.txt", append = TRUE)
  #enregistrer les probes
  capture.output(i, file="./results/Tests_statistiques/Significant_probes/probes_BaselinevsCtrl.txt", append = TRUE)
}
#rmq = les fichiers ont ete modifiés à la main par la suite.

```

- DAVID et Panther et Gorilla ont été utilisés (voir dans le répertoire results/)