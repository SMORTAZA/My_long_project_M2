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
