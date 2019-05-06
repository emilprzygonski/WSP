setwd("C:/Users/superstudent/Desktop/lab_wsp")
#Po zainstalowaniu programu R należy zainstalować Bioconductor poleceniem:
source("http://bioconductor.org/biocLite.R")#sciąga skrypt instalacyjny
biocLite() #instaluje domyślny zestaw pakietów
# lista wszystkich dostępnych pakietów http://www.bioconductor.org/packages/release/BiocViews.html #___Software
biocLite("BiocUpgrade") # uaktualnia, istnieje również możliwość aktualizacji
# PRZYKŁADOWE ustawienie ścieżki katalogu roboczego uwaga używamy "/" jak w linuxie
library(Biobase)
biocLite('affy') # standardowo nie jest zainstalowana ta biblioteka biocLite('affy')
library('affy')
biocLite('gahgu95av2.db')
library(gahgu95av2.db) # standardowo nie jest zainstalowana ta biblioteka biocLite(gahgu95av2.db)
biocLite("gplots")
install.packages(gplots)
library(gplots) # heatmap.2
### tworzenie pełnego obiektu ExpressionSet (zgodnie z instrukcją ExpressionSetIntroduction.pdf)
### Wczytanie danych wykorzystując plik opisu danych datasetA.txt, który może zawierać opis kolumn na początku pliku
exampleFile = system.file("extdata", "pData.txt", package="Biobase")
data = read.table("datasetA_scans.txt", header = TRUE, sep = "\t")
data=data[c(1:5,244:248),] # do testów lepiej pracować na małej liczbie mikromacierzy
opis = read.AnnotatedDataFrame("datasetA_scans.txt", sep="\t", header=TRUE, row.names=4,
                               stringsAsFactors = F) 
opis=opis[c(1:5,244:248)] 
sampleNames(opis) = paste(sampleNames(opis), ".CEL", sep="") 
data_Affy=ReadAffy(filenames=sampleNames(opis), verbose=TRUE) 
data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="") 
data_Affy@annotation=paste("ga",data_Affy@annotation,sep="") 
RMA=rma(data_Affy) 
dataRMA=exprs(RMA) 

#tworzenie obiektu ExpressionSet
experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane tesowe",
                 abstract = "Przyklad",url = "http://www.bioconductor.org", other = list(notes = "inne"))
ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opis,
             experimentData=experiment,annotation="gahgu95av2.db") 
expr_sort=sort(rowMeans(exprs(ExprSet)),index.return=T) 
feat_num=dim(ExprSet)[1] 
cutoff=round(dim(ExprSet)[1]*0.025)
ind_clear=expr_sort$ix[c(1:cutoff,(feat_num-cutoff):feat_num)]
ExprSet=ExprSet[-ind_clear,]
PCA_model=prcomp(t(exprs(ExprSet)))
summary(PCA_model) 
PCA_model$x 
adeno=which(pData(ExprSet)$CLASS=='ADENO') #wyszukujemy próby ADENO
squamous=which(pData(ExprSet)$CLASS=='SQUAMOUS') # analogicznie do SQUAMOUS
colors=ifelse(pData(ExprSet)$CLASS=='ADENO', 'red', 'blue') #dobieramy kolory
plot(PCA_model$x[,1:2], col=colors, main='PCA') #wykres przetworzonych danych w nowym układziewspółrzędnych (bazie własnej)
barplot(PCA_model$sdev[1:5]/sum(PCA_model$sde),main='PCA') 