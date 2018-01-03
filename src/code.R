##########################################
# Chargement et visualisation des données#
##########################################
#fusionner avec le code du readme
###########################
#Normalisation des données#
###########################
#fusionner avec le code du readme
######################################################
#Identification des gènes différentiellement exprimés#
######################################################
#vérifier le tableau contenant toutes les données
dim(data)
#rearrangement des col utiles de data
w_data=data[,c(1,6,11,16,21,2,7,12,17,22,3,8,13,18,23,4,9,14,19,24,5,10,15,20,25)]
#Enregistrement du nouveau tableau
save(w_data, file = "./results/w_data.RData")
#ouverture de ces données
load('./results/data.RData')
load('./results/w_data.RData')
#creation du design
design <- matrix(c(rep(c(1,0,0,0,0),5),rep(c(0,1,0,0,0),5),rep(c(0,0,1,0,0),5),rep(c(0,0,0,1,0),5),rep(c(0,0,0,0,1),5)), ncol=5,byrow=TRUE)
colnames(design) <- c("Baseline","Ctrl","HNO3","M1","M2")
rownames(design) <- colnames(w_data)
design <- data.frame(design)

########################################################
#Pour avoir un seuil pour la significativité des pvalues
########################################################
#####the t-test#####
require(genefilter)
#en chargeant w_data, son typeof est list o
#typeof(data)
#modifier le typeof en double
w_data <- data.matrix(w_data)
#verification de la modif
#typeof(data)
#calcul des pvalues
t.pval.CtrlvsHNO3 <- rowttests(w_data[,seq(6,15)],factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
t.pval.HNO3vsM1 <- rowttests(w_data[,seq(11,20)],factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
t.pval.HNO3vsM2 <- rowttests(w_data[,c(11,12,13,14,15,21,22,23,24,25)],factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
#regarder les top10 genes selon les pvalues
sort(t.CtrlvsHNO3)[1:10]
sort(t.HNO3vsM1)[1:10]
sort(t.HNO3vsM2)[1:10]
#annotation
#require(annotate)
#library(pd.hugene.2.0.st)
#AffyID <- rownames(w_data)
#ne donne que des NA ?
#GeneNames <- getSYMBOL( AffyID , "hugene20sttranscriptcluster" )
#to know how genes are affected (up or down), calculate the log2 of the fold-change
fc.CtrlvsHNO3 <- rowMeans(w_data[,design$Ctrl==1])-rowMeans(w_data[,design$HNO3==1])
fc.HNO3vsM1 <- rowMeans(w_data[,design$HNO3==1])-rowMeans(w_data[,design$M1==1])
fc.HNO3vsM2 <- rowMeans(w_data[,design$HNO3==1])-rowMeans(w_data[,design$M2==1])
#Volcanoplot pr determiner un seuil de significativité pour pvalue
# a volcano plot shows the connection between the pvalues and the log2 of the fold change, 
#compared to the same to the same analysis of the same data randomized. the permutation 
#should be a balanced randomization of the colums of the experiment. 
#After the permutation we run the ttest again and calculate a new log2
#fold change (M) value. Finally, we plot the two outputs on top of each
#other in order to compare.
#So, in order to make a volcanoplot, we need to obtain the random or 
#permutated pvalues for each dataset. To get those, we need to perform a ttest on the
#same data, but with permutated labels. Similarly for the "permutated"
#fold-change. Finally, we'll plot the "permutated" values vs permutated
#fold-change (red dots) on top of the "real" data.

#avoir le nb de lignes et de colonnes de la matrice design
dim(design)
#creation de la matrice permutated de la meme dimension que la matrice design
permutated_design = matrix(nrow = 25, ncol = 5)
for (i in seq(dim(design)[2])){
  print(i)
  print(design[,i])
  random_lables <- sample(design[,i])
  print(random_lables)
  permutated_design[,i] <- random_lables
}
colnames(permutated_design) <- c("Baseline","Ctrl","HNO3","M1","M2")
rownames(permutated_design) <- colnames(w_data)
permutated_design <- data.frame(permutated_design)
#obtention des col permutées fom w_data
Ctrl_permutated <- w_data[,c(permutated_design$Ctrl==1)]
HNO3_permutated <- w_data[,c(permutated_design$HNO3==1)]
M1_permutated <- w_data[,c(permutated_design$M1==1)]
M2_permutated <- w_data[,c(permutated_design$M2==1)]
#permutated pval
perm_t.pval.CtrlvsHNO3 <- rowttests(cbind(Ctrl_permutated,HNO3_permutated), factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
perm_t.pval.HNO3vsM1 <- rowttests(cbind(HNO3_permutated,M1_permutated), factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
perm_t.pval.HNO3vsM2 <- rowttests(cbind(HNO3_permutated,M2_permutated), factor(c(1,1,1,1,1,2,2,2,2,2)))$p.value
#permutated fc
perm_fc.CtrlvsHNO3 <- rowMeans(w_data[,permutated_design$Ctrl==1])-rowMeans(w_data[,permutated_design$HNO3==1])
perm_fc.HNO3vsM1 <- rowMeans(w_data[,permutated_design$HNO3==1])-rowMeans(w_data[,permutated_design$M1==1])
perm_fc.HNO3vsM2 <- rowMeans(w_data[,permutated_design$HNO3==1])-rowMeans(w_data[,permutated_design$M2==1])
#plot du volcano de la premiere comparaison -> pr choisir le seuil, faut que sous ce cut-off, il y ait 
#plus de points bleus (données reelles) que de points rouges (données simulées) car c'est là où
#les gènes vont être le plus significatif.
plot(fc.CtrlvsHNO3, t.pval.CtrlvsHNO3, main = "Volcano Plot\nCtrl vs HNO3", log = "y",xlab = "M(log2 fold change)", ylab = "p-value", pch = 20, col = "blue")
points(perm_fc.CtrlvsHNO3, perm_t.pval.CtrlvsHNO3, type = "p", pch = 20, col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"),pch=20,bg="white")
abline(h=0.00025)
#plot du volcano de la seconde comparaison
plot(fc.HNO3vsM1, t.pval.HNO3vsM1, main = "Volcano Plot\nHNO3 vs M1", log = "y",xlab = "M(log2 fold change)", ylab = "p-value", pch = 20, col = "blue")
points(perm_fc.HNO3vsM1, perm_t.pval.HNO3vsM1, type = "p", pch = 20, col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"),pch=20,bg="white")
abline(h=0.01)
#plot du volcano de la troisieme comparaison
plot(fc.HNO3vsM2, t.pval.HNO3vsM2, main = "Volcano Plot\nHNO3 vs M2", log = "y",xlab = "M(log2 fold change)", ylab = "p-value", pch = 20, col = "blue")
points(perm_fc.HNO3vsM2, perm_t.pval.HNO3vsM2, type = "p", pch = 20, col = "red")
legend("topleft", col=c("red","blue"), legend=c("perm", "real"),pch=20,bg="white")
abline(h=0.04)

#############################mise a part pour l'instant
#importation de la librairie limma
library(limma)
#conditions à comparer
contrast.matrix <- makeContrasts(Ctrl-HNO3, HNO3-M1, HNO3-M2, levels=design)
#paramètres du modèle ajustés en fct de nos données
fit <- lmFit(log2(w_data), design)
#recherche de differences significatives
eBayesResultat <- eBayes(contrasts.fit(fit, contrast.matrix))
#head des pvalues
CtrlvsHNO3 = eBayesResultat$p.value[,1]
HNO3vsM1 = eBayesResultat$p.value[,2]
HNO3vsM2 = eBayesResultat$p.value[,3]


