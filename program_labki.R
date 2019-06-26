source("https://bioconductor.org/biocLite.R")

BiocInstaller::biocLite(pkgs=c("limma", "Biobase", "affy", "gplots","gahgu95av2.db","org.Hs.eg.db","gahgu95av2cdf"), ask=FALSE)
library(limma)
library(Biobase)
library('affy')
library(gplots) 
library(gahgu95av2.db) # standardowo nie jest zainstalowana ta biblioteka biocLite(gahgu95av2.db)
library('gahgu95av2.db')
library('org.Hs.eg.db')
library('gahgu95av2cdf')

setwd("C:/Users/Wiki/Documents/projekt R") 

data = read.table("datasetA_scans.txt", header = TRUE, sep = "\t")
opisy = read.AnnotatedDataFrame("datasetA_scans.txt", sep="\t", header=TRUE, row.names=4,stringsAsFactors = F) 
experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))

# data1=data[1:10,]
# data1=na.omit(data1)
# opis=opisy[1:10]
# opis=na.omit(opis)
sampleNames(opisy) = paste(sampleNames(opisy), ".CEL", sep="")
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
updated_ExprSet=function(ExprSet){
  
  ekspresje=exprs(ExprSet)
  
  symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) 
  
  genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME))
  
  entrezy=unlist(mget(rownames(ekspresje),as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
  
  entrezy_nazwy=unlist(mget(rownames(ekspresje),as.environment(as.list(gahgu95av2GENENAME)),ifnotfound=NA))
 
  
  macierz=data.frame(unlist(featureNames(ExprSet)),unlist(symbol), unlist(genNames),unlist(entrezy), unlist(entrezy_nazwy))
  names(macierz)[1]="próbka"
  names(macierz)[2]="symbol"
  names(macierz)[3]="nazwa"
  names(macierz)[4]="entrez_id"
  names(macierz)[5]="entrez_nazwa"
  
  opisy=ExprSet@phenoData
  experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))
  ExprSet= new("ExpressionSet", expr=ekspresje, phenoData = opisy,experimentData=experiment,annotation="gahgu95av2.db",featureData=AnnotatedDataFrame(macierz))
  
  return(ExprSet)
}

ExprSet=updated_ExprSet(ExprSet)

#eksport do pliku Excela
install.packages("openxlsx")
library("openxlsx")
installr::install.rtools() #R.3.3.x or later
Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")
cechy=ExprSet@featureData
cechy=cechy@data
write.xlsx(cechy, "dane_geny.xlsx", asTable = TRUE)

#analiza cech różnicujących - test t lub Fold Change
summary_table_FC=function(ExprSet, correction_method, sort_method, sort_criterion, sep){
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
  p_val_skorygowane=p.adjust(pval, method = correction_method) # metoda korekcji pobierana jako parametr funkcji
  #SYMBOLE GEN?W
  symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) #mo?na jako lista albo wektor
  #NAZWY GEN?W
  genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME)) # podobnie lista przekszta?cona na wektor
  
  entrezy=unlist(mget(rownames(expr),as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
  
  TAB=array(dim=c(dim(expr)[1],11))
  colnames(TAB)=c("Indeks","FerrariID","Symbol","description","entrez id","fold change","srednia w gr.ADENO","srednia w gr.NORMAL","t-statistics","pvalue","corrected p-value")
  #tabela dla wszystkich sond, przed selekcj? mo?na ju? wcze?niej zainicjowa? i bezpo?rednio wpisywa?
  TAB[,1]=c(1:dim(expr)[1])
  TAB[,2]=featureNames(ExprSet)
  TAB[,3]=symbol
  TAB[,4]=genNames
  TAB[,5]=entrezy
  TAB[,6]=FC
  TAB[,7]=sr_adeno
  TAB[,8]=sr_squamous
  TAB[,9]=statistic
  TAB[,10]=pval
  TAB[,11]=p_val_skorygowane
  #wyb?r sond do tabeli wynikowej zabezpieczenie mo?na zrobi? na podanie z?ej warto?ci (ujemnej)
  # je?eli mamy u?amek to warto?? >1 to ilo?? sond zwracamy =1 to wszystkie
  #je?eli warto??
  if (sort_method==1){
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
    ind_sort=sort(as.double(TAB[,11]),index=TRUE)$ix
    TAB=TAB[ind_sort,]
  }
  
  if (sort_method==2){
    if (sort_criterion>5)
    {
      ind_sort=sort(FC,decreasing = TRUE, index=TRUE)$ix #interesuj? mnie indeksy
      TAB=TAB[ind_sort[1:sort_criterion],] # tablica TAB zawiera dane dla n sond o najni?szym pvalue
    } # tutaj dane s? posortowane co w niczym nie przeszkadza
    if (sort_criterion<5)
    {
      ind_sort=which(FC>sort_criterion)
      TAB=TAB[ind_sort,]
    }
    # czy w?a?ciwe sortowanie, 0 brak sortowania
    ind_sort=sort(as.double(TAB[,6]),decreasing = TRUE,index=TRUE)$ix
    TAB=TAB[ind_sort,]
  }
  
  # design<- model.matrix(~c(rep(1,1000),rep(0,100)))
  # fit = lmFit(expr , design)
  # fit = eBayes(fit)
  # tt = topTable(fit, coef = 2, adjust.method = "BH", sort.by = "p", number = 10)
  
  # zapisanie tablicy wynikowej
  write.table(TAB, file="summary_table.txt",sep=sep,row.names = FALSE) #separator pobierany jako paramer
  return(TAB)
} 

#przyk?adowe wywo?anie
#correction_method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#sort_method: 1 - corrected p_val, 2 - Fold Change
TAB_geny_roznicujace=summary_table_FC(ExprSet, correction_method='holm', sort_method=1 ,sort_criterion=0.01, sep=', ') #correction_method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")

#tworzenie heatmapy
#wartosci ekspresji dla wybranych do tabeli gen?w
expr=exprs(ExprSet)
ind=TAB_geny_roznicujace[,1]
ekspr_wybrane=expr[as.double(ind),] # tablica ekspresji tylko dla wybranych zestaw?w sond
png("heatmap.png") # tworzymy rysunek mapy ciep?a
heatmap.2(ekspr_wybrane) # rozbudowana wersja mapy ciep?a: mo?na poustawia? r??ne parametry
dev.off() #zamykamy okno

#wszystkie oraz różnicujące geny do analizy ścieżek sygnalnych np n podstawie entrez id
wszystkie_geny=ExprSet@featureData
wszystkie_geny=wszystkie_geny@data
wszystkie_geny_entrez_id=wszystkie_geny$entrez_id

roznicujce_geny_entrez_id=TAB_geny_roznicujace_ttest[,4]
