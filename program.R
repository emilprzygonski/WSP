source("https://bioconductor.org/biocLite.R")
ELOTUARTUR=1

library(Biobase)
biocLite('affy')
library('affy')

biocLite('gahgu95av2.db')
library(gahgu95av2.db) # standardowo nie jest zainstalowana ta biblioteka biocLite(gahgu95av2.db)
biocLite("gplots")
library(gplots) 

library('gahgu95av2.db')
library('org.Hs.eg.db')
library('gahgu95av2cdf')

setwd("C:\Users\superstudent\Desktop\WSP2019") 

exampleFile = system.file("extdata", "pData.txt", package="Biobase")
data = read.table("datasetA_scans.txt", header = TRUE, sep = "\t")
data=data[c(1:5,244:248),]
opis = read.AnnotatedDataFrame("datasetA_scans.txt", sep="\t", header=TRUE, row.names=4,stringsAsFactors = F) 
opis=opis[c(1:5,244:248)]
sampleNames(opis) = paste(sampleNames(opis), ".CEL", sep="")
data_Affy=ReadAffy(filenames=sampleNames(opis), verbose=TRUE)
data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="")
data_Affy@annotation=paste("ga",data_Affy@annotation,sep="")

<<<<<<< HEAD
hello=1
proba <- 3
=======
>>>>>>> c953499ad8c592e4250a3ecedbd1698542bf2ce4

RMA=rma(data_Affy) 

dataRMA=exprs(RMA) 


experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane tesowe",abstract = "Przyklad",url = "http://www.bioconductor.org", other = list(notes = "inne"))

ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opis,experimentData=experiment,annotation="gahgu95av2.db") 

expr_sort=sort(rowMeans(exprs(ExprSet)),index.return=T) #sortujemy średnie ekspresje 

feat_num=dim(ExprSet)[1] 

cutoff=round(dim(ExprSet)[1]*0.025) #wyznaczamy ilość sond do usunięcia po obu stronach
ind_clear=expr_sort$ix[c(1:cutoff,(feat_num-cutoff):feat_num)]
ExprSet=ExprSet[-ind_clear,]


PCA_model=prcomp(t(exprs(ExprSet))) 
summary(PCA_model) 
PCA_model$x


adeno=which(pData(ExprSet)$CLASS=='ADENO') #wyszukujemy próby ADENO
squamous=which(pData(ExprSet)$CLASS=='SQUAMOUS') # analogicznie do SQUAMOUS
colors=ifelse(pData(ExprSet)$CLASS=='ADENO', 'red', 'blue') #dobieramy kolory
plot(PCA_model$x[,1:2], col=colors, main='PCA')

barplot(PCA_model$sdev[1:5]/sum(PCA_model$sde),main='PCA') 


symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) 

genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME))

entrezy=mget(rownames(dataRMA),as.environment(as.list(gahgu95av2ENTREZID)))

entrezy_nazwy=mget(rownames(dataRMA),as.environment(as.list(gahgu95av2GENENAME)))

macierz=data.frame(unlist(entrezy),unlist(entrezy_nazwy))
names(macierz)[1]="entrez_id"
names(macierz)[2]="entrez_nazwy"

ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opis,experimentData=experiment,annotation="gahgu95av2.db",featureData=AnnotatedDataFrame(macierz))
