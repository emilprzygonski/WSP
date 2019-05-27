source("https://bioconductor.org/biocLite.R")

library(Biobase)
biocLite('affy')
library('affy')

biocLite("gplots")
library(gplots) 

setwd("C:/Users/superstudent/Downloads/projekt R") 

exampleFile = system.file("extdata", "pData.txt", package="Biobase")
data = read.table("datasetA_scans.txt", header = TRUE, sep = "\t")
data=data[c(1:5,244:248),]
opis = read.AnnotatedDataFrame("datasetA_scans.txt", sep="\t", header=TRUE, row.names=4,stringsAsFactors = F) 
opis=opis[c(1:5,244:248)]
sampleNames(opis) = paste(sampleNames(opis), ".CEL", sep="")
data_Affy=ReadAffy(filenames=sampleNames(opis), verbose=TRUE)
data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="")
data_Affy@annotation=paste("ga",data_Affy@annotation,sep="")


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

#funkcja dodaj?ca entrez id do ExprSet
updated_ExprSet=function(ExprSet, dataRMA){
  
  biocLite('gahgu95av2.db')
  library(gahgu95av2.db) # standardowo nie jest zainstalowana ta biblioteka biocLite(gahgu95av2.db)
  library('gahgu95av2.db')
  library('org.Hs.eg.db')
  library('gahgu95av2cdf')
  
  symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) 
  
  genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME))
  
  entrezy=mget(rownames(dataRMA),as.environment(as.list(gahgu95av2ENTREZID)))
  
  entrezy_nazwy=mget(rownames(dataRMA),as.environment(as.list(gahgu95av2GENENAME)))
  
  macierz=data.frame(unlist(entrezy),unlist(entrezy_nazwy))
  names(macierz)[1]="entrez_id"
  names(macierz)[2]="entrez_nazwy"
  
  ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opis,experimentData=experiment,annotation="gahgu95av2.db",featureData=AnnotatedDataFrame(macierz))
  
  return(ExprSet)
}

ExprSet=updated_ExprSet(ExprSet, dataRMA)
######
summary_table=function(ExprSet, method, sort_criterion, col_nr, sep){
  adeno=which(pData(ExprSet)$CLASS=='ADENO') 
  squamous=which(pData(ExprSet)$CLASS=='SQUAMOUS') 
  expr=exprs(ExprSet)
  sr_adeno=rowMeans(expr[,adeno])
  sr_squamous=rowMeans(expr[,squamous])
  # fold change
  FC=sr_adeno/sr_squamous
  #test t
  statistic=apply(expr,1,function(x) t.test(x[adeno],x[squamous])$statistic)
  #przydatne polecenia apply vapply sapply..
  pval=apply(expr,1,function(x) t.test(x[adeno],x[squamous])$p.val) # p warto?ci
  #Korekta na wielokrotne testowanie
  p_val_skorygowane=p.adjust(pval, method = method) # metoda korekcji pobierana jako parametr funkcji
  #SYMBOLE GEN?W
  symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) #mo?na jako lista albo wektor
  #NAZWY GEN?W
  genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME)) # podobnie lista przekszta?cona na wektor
  TAB=array(dim=c(dim(expr)[1],9))
  colnames(TAB)=c("FerrariID","Symbol","description","fold change","?rednia w gr.ADENO","?rednia w
                  gr.NORMAL","t-statistics",
                  "pvalue","corrected p-value")
  #tabela dla wszystkich sond, przed selekcj? mo?na ju? wcze?niej zainicjowa? i bezpo?rednio wpisywa?
  TAB[,1]=featureNames(ExprSet)
  TAB[,2]=symbol
  TAB[,3]=genNames
  TAB[,4]=FC
  TAB[,5]=sr_adeno
  TAB[,6]=sr_squamous
  TAB[,7]=statistic
  TAB[,8]=pval
  TAB[,9]=p_val_skorygowane
  head(TAB) #sprawdzenie jak wygl?da nasza tabela
  #wyb?r sond do tabeli wynikowej zabezpieczenie mo?na zrobi? na podanie z?ej warto?ci (ujemnej)
  # je?eli mamy u?amek to warto?? >1 to ilo?? sond zwracamy =1 to wszystkie
  #je?eli warto??
  if (sort_criterion>1)
  {
    ind_sort=sort(p_val_skorygowane,index=TRUE)$ix #interesuj? mnie indeksy
    TAB=TAB[ind_sort[1:sort_criterion],] # tablica TAB zawiera dane dla n sond o najni?szym pvalue
  } # tutaj dane s? posortowane co w niczym nie przeszkadza
  if (sort_criterion<1)
  {
    ind_sort=which(p_val_skorygowane<sort_criterion)
    TAB=TAB[ind_sort,]
  }
  # czy w?a?ciwe sortowanie, 0 brak sortowania
  if (col_nr!=0){
    ind_sort=sort(TAB[col_nr,],index=TRUE)$ix
    TAB=TAB[ind_sort,]
  }
  # zapisanie tablicy wynikowej
  write.table(TAB, file="data_table.txt",sep=sep,row.names = FALSE) #separator pobierany jako paramer
  #tworzenie heatmapy
  #wartosci ekspresji dla wybranych do tabeli gen?w
  ekspr_wybrane=expr[ind_sort,] # tablica ekspresji tylko dla wybranych zestaw?w sond
  png("heatmap.png") # tworzymy rysunek mapy ciep?a
  heatmap.2(ekspr_wybrane) # rozbudowana wersja mapy ciep?a: mo?na poustawia? r??ne parametry
  dev.off() #zamykamy okno
  return(TAB)
} #koniec funkcji

#przyk?adowe wywo?anie
summary_table(ExprSet, method='holm', sort_criterion=15, col_nr=5, sep=',') #method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")