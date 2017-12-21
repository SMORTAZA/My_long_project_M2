##########################################
# Chargement et visualisation des données#
##########################################

###########################
#Normalisation des données#
###########################

######################################################
#Identification des gènes différentiellement exprimés#
######################################################
#vérifier le tableau contenant toutes les données
dim(data)
#rearrangement des col utiles de data
w_data=data[,c(1,6,11,16,21,2,7,12,17,22,3,8,13,18,23,4,9,14,19,24,5,10,15,20,25)]
#Enregistrement du nouveau tableau
save(w_data, file = "./results/w_data.RData")
#creation du design
design <- matrix(c(rep(c(1,0,0,0,0),5),rep(c(0,1,0,0,0),5),rep(c(0,0,1,0,0),5),rep(c(0,0,0,1,0),5),rep(c(0,0,0,0,1),5)), ncol=5,byrow=TRUE)
colnames(design) <- c("Baseline","Ctrl","HNO3","M1","M2")
rownames(design) <- colnames(w_data)

#importation de la librairie limma
library(limma)
#conditions à comparer
contrast.matrix <- makeContrasts(Ctrl-HNO3, HNO3-M1, HNO3-M2, levels=design)
#paramètres du modèle ajustés en fct de nos données
fit <- lmFit(log2(w_data), design)
#recherche de differences significatives
eBayesResultat <- eBayes(contrasts.fit(fit, contrast.matrix))
#head des pvalues
head(eBayesResultat$p.value)