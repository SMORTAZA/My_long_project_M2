# Notes pour l'avancement du projet long

##10/10/17
- rdv Karine Audouze
- téléchargement des données de Martine Aggerbeck

##23/10/2017
- Les données à analyser se trouvent dans Raw Data_juillet 2017_HepaRG_phtalates et metals_UPD_E16_093-2bis. 
- pb = enlever le fichier word du dossier contenant les data. Mais erreur à résoudre

##25/10/2017
- https://bioconductor.org/packages/release/bioc/html/xps.html
- installation de xps pr voir => ne s'installe pas correctement
- essayer avec oligo
https://bioconductor.org/packages/release/bioc/html/oligo.html
```{r}
- le pb vient de la librairie Affy et non de list.celfiles()
Erreur : 

The affy package is not designed for this array type.
Please use either the oligo or xps package.
```
- http://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor#Open_CEL_files_from_newer_Affymetrix_Arrays_.28HTA.2C_Gene_ST....29_using_oligo

##26/10/2017
- http://homer.ucsd.edu/homer/basicTutorial/affymetrix.html
- installation d'un package pour l'annotation (n'a pas servi)
- Affymetrix Gene ST Arrays 
- (https://www.bioconductor.org/help/workflows/arrays/#affymetrix-gene-st-arrays)
- pour choisir l'array specifique des données, aller sur le site : 
 https://www.bioconductor.org/packages/release/data/annotation/
- Dans l'article, (Affymetrix GeneChip HuGene 2_0-st) en page 9.
https://bioconductor.org/packages/release/data/annotation/html/pd.hugene.2.0.st.html

##30/10/2017
- creation de nvx dossiers 
- http://rug.mnhn.fr/semin-r/PDF/semin-R_data_JPedraza_100608.pdf
pr enregistrement la dataframe intermédiaire

##31/10/2017
- page 6 de https://bioconductor.org/packages/release/data/annotation/manuals/hugene20sttranscriptcluster.db/man/hugene20sttranscriptcluster.db.pdf (site donné par Karine AUDOUZE pour l'annotation des gènes)
- http://www.ensembl.org/info/data/biomart/biomart_r_package.html pr aller de id ensembl à noms de genes

##02/11/2017
- pb au niv annotation genes pr avoir identifiant EnsemblGeneID (pas pr certains ?)
Error in if (which(a == line_num[i]) == integer(0)) { : 
  l'argument est de longueur nulle

##07/11/2017
- pas de nom de genes pour la plupart des 'probes' -> ?

##10/11/2017
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

##13/11/2017
> length(which(data$geneSymbol == "NA"))
[1] 27018
> length(which(data$geneSymbol != "NA"))
[1] 26599

##14/11/2017
- 53617 probes au total ds la bdd
> length(contents(hugene20sttranscriptclusterSYMBOL))
[1] 53617

##28/11/2017
- envoi de mail pour analyse des SVD/PCA
- note : choix de anova plutot que de ttest car plus puissant
- essayer de comprendre comment marche Annova 2 ways.

##04/12/2017
- comprendre anova

##05/12/2017
- creer des design pour chaque condition ... ?
- demande aide à un bioinformaticien de C3BI pour faire anova, m'a conseillé d'utiliser package limma, plus adapté à mes données

##21/12/2017
- http://silico.biotoul.fr/site/index.php/InfoBio_TD_transcriptome
- http://bioinfo.unil.ch/tp/bachelor_3rd_year/tp1/
- ce dernier site detaille bien les etapes a suivre. dc reussi à obtenir des pvalues pour les trois comparaisons. reste à faire sort() dessus, puis de faire volcanoplot.
- commennt avoir le fold change pr faire le volcanoplot ?
