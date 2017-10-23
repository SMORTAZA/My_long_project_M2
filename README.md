#README

##Intégration de données transcriptomiques et protéomiques dans un modèle d'hépatocytes humains

```{R}
library(affy)
setwd("/home/sdv/m2bi/smortaza/Bureau/Projet_Long/Raw Data_juillet 2017_HepaRG_phtalates plus metals_UPD_E16_093-2bis")
getwd()
data <- ReadAffy(filenames=sort(list.celfiles("celfiles", full=TRUE)))
#setwd("/home/sdv/m2bi/smortaza/Bureau/Projet_Long")
#getwd()
```
Erreur : 

The affy package is not designed for this array type.
Please use either the oligo or xps package.