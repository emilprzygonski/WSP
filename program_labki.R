source("https://bioconductor.org/biocLite.R")

biocLite('Biobase')
library(Biobase)
biocLite('affy')
library('affy')

biocLite("gplots")
library(gplots) 

biocLite('gahgu95av2.db')
library(gahgu95av2.db) # standardowo nie jest zainstalowana ta biblioteka biocLite(gahgu95av2.db)
library('gahgu95av2.db')
biocLite('org.Hs.eg.db')
library('org.Hs.eg.db')
biocLite('gahgu95av2cdf')
library('gahgu95av2cdf')

setwd("C:/Users/Wiki/Documents/projekt R") 

data = read.table("datasetA_scans.txt", header = TRUE, sep = "\t")
opisy = read.AnnotatedDataFrame("datasetA_scans.txt", sep="\t", header=TRUE, row.names=4,stringsAsFactors = F) 
experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))

# data1=data[1:10,]
# data1=na.omit(data1)
# opis=opisy[1:10]
# opis=na.omit(opis)
sampleNames(opis) = paste(sampleNames(opisy), ".CEL", sep="")
data_Affy=ReadAffy(filenames=sampleNames(opisy), verbose=TRUE)
data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="")
data_Affy@annotation=paste("ga",data_Affy@annotation,sep="")

RMA=rma(data_Affy)
dataRMA=exprs(RMA) 


# x = seq(11,length(data[,1]),10)

# for (i in x){
# data1=data[i:(i+9),]
# data1=na.omit(data1)
# opis1=opisy[i:(i+9)]
# opis1=na.omit(opis1)
# 
# sampleNames(opis1) = paste(sampleNames(opis1), ".CEL", sep="")
# data_Affy=ReadAffy(filenames=sampleNames(opis1), verbose=TRUE)
# data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="")
# data_Affy@annotation=paste("ga",data_Affy@annotation,sep="")
# 
# RMA1=rma(data_Affy)
# dataRMA1=exprs(RMA1) 
# dataRMA=cbind(dataRMA, dataRMA1)
# }

ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opisy,experimentData=experiment,annotation="gahgu95av2.db") 

#funkcja dodajaca entrez id do ExprSet
updated_ExprSet=function(ExprSet, dataRMA){
  
  symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) 
  
  genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME))
  
  entrezy=unlist(mget(rownames(dataRMA),as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
  
  entrezy_nazwy=unlist(mget(rownames(dataRMA),as.environment(as.list(gahgu95av2GENENAME)),ifnotfound=NA))
 
  
  macierz=data.frame(unlist(featureNames(ExprSet)),unlist(symbol), unlist(genNames),unlist(entrezy), unlist(entrezy_nazwy))
  names(macierz)[1]="próbka"
  names(macierz)[2]="symbol"
  names(macierz)[3]="nazwa"
  names(macierz)[4]="entrez_id"
  names(macierz)[5]="entrez_nazwa"
  
  ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opisy,experimentData=experiment,annotation="gahgu95av2.db",featureData=AnnotatedDataFrame(macierz))
  
  return(ExprSet)
}

ExprSet=updated_ExprSet(ExprSet, dataRMA)

install.packages("openxlsx")
library("openxlsx")

cechy=ExprSet@featureData
cechy=cechy@data
wb <- createWorkbook()
addWorksheet(wb, "Geny")
writeData(wb, "Geny", cechy, startCol = 1, startRow = 1, rowNames = FALSE, colNames=TRUE)



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
