#README

##Intégration de données transcriptomiques et protéomiques dans un modèle d'hépatocytes humains

```{r}
library(affy)
getwd()
setwd("/home/sdv/m2bi/smortaza/Bureau/Projet_Long/Raw Data_juillet 2017_HepaRG_phtalates plus metals_UPD_E16_093-2bis")
getwd()
data <- ReadAffy(filenames=sort(list.celfiles()))
#setwd("/home/sdv/m2bi/smortaza/Bureau/Projet_Long")
#getwd()
```

```{r}
celpath = "/home/sdv/m2bi/smortaza/Bureau/Projet_Long/Raw Data_juillet 2017_HepaRG_phtalates plus metals_UPD_E16_093-2bis"
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)
```