source("https://bioconductor.org/biocLite.R")

BiocInstaller::biocLite(pkgs=c("limma", "Biobase", "affy", "gplots","gahgu95av2.db","org.Hs.eg.db","gahgu95av2cdf"), ask=FALSE)
library(limma)
library(Biobase)
library('affy')
library(gplots) 
library('gahgu95av2.db')
library('org.Hs.eg.db')
library('gahgu95av2cdf')

create_ExprSet=function(){
  
  #funkcja tworząca obiekt typu ExpressionSet na podstawie plików .CEL, 
  #których nazwy i dane znajdują się w pliku z adnotacjami. Plik ten powinien zostać wskazany przez użytkownika.
  
  require(Biobase)
  require(affy)
  require(tcltk)
  require('gahgu95av2.db')
  
  #wybór folderu z plikami .CEL
  choose.dir(getwd(), "Wybierz folder z plikami .CEL")
  #wybór pliku z adnotacjami
  file=tk_choose.files(caption = "Wybierz plik .txt z adnotacjami")
  #wczytnie danych
  data = read.table(file, header = TRUE, sep = "\t")
  opisy = read.AnnotatedDataFrame(file, sep="\t", header=TRUE, row.names=4,stringsAsFactors = F) 
  experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))
  
  #nazwy plików.CEL
  sampleNames(opisy) = paste(sampleNames(opisy), ".CEL", sep="")
  #wczytanie plików i tworzenie macierzy ekspresji
  data_Affy=ReadAffy(filenames=sampleNames(opisy), verbose=TRUE)
  data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="")
  data_Affy@annotation=paste("ga",data_Affy@annotation,sep="")
  
  #normalizacja 
  RMA=rma(data_Affy)
  dataRMA=exprs(RMA) 
  
  #tworzenie nowego obiektu ExpressionSet
  ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opisy,experimentData=experiment,annotation="gahgu95av2.db") 
  
  return(ExprSet)
}

#wywołanie
ExprSet=create_ExprSet() 

#wczytanie pliku tekstowego z nazwami sond wybranego przez użytkownika
wczytaj_sondy=function(){
  
  require(tcltk)
  
  #wczytanie pliku tekstowego z nazwami sond wybranego przez użytkownika
  sondy=lapply(read.delim(tk_choose.files(caption = "Wybierz plik .txt z nazwami sond"), header = FALSE, stringsAsFactor = FALSE), as.character)$V1
  
  return(sondy)
}

#wywołanie
sondy=wczytaj_sondy()

#wczytanie Expressionset (tworzenie z gotowej macierzy ekspresji)
create_ExprSet2=function(){
  
  require(Biobase)
  require(affy)
  require(tcltk)
  require('gahgu95av2.db')
  
  exprsFile =read.table(tk_choose.files(caption = "Wybierz plik .txt z macierzą ekspresji"))
  #wybór pliku z adnotacjami
  phenoFile = tk_choose.files(caption = "Wybierz plik .txt z adnotacjami")
  phenoFile = read.AnnotatedDataFrame(phenoFile, sep="\t", header=TRUE, row.names=4,stringsAsFactors = F)
  experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))

  #tworzenie nowego obiektu ExpressionSet
  ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opisy,experimentData=experiment,annotation="gahgu95av2.db") 

  return(ExprSet)
}

ExprSet=create_ExprSet2()

updated_ExprSet=function(ExprSet){
  #funkcja dodajaca symbole genów, ich nazwy oraz entrez id do ExprSet
  
  require(Biobase)
  require(affy)
  require('gahgu95av2cdf')
  require('gahgu95av2.db')
  require('org.Hs.eg.db')
  
  if (is.character(ExprSet)=='TRUE'){ #gdy użytkownik wczytuje jedynie litę z nazwami sond - zwracany ExprSet z zerowymi ekspresjami
    
    symbol=unlist(mget(ExprSet,env=gahgu95av2SYMBOL)) 
    
    genNames=unlist(mget(ExprSet,env=gahgu95av2GENENAME))
    
    entrezy=unlist(mget(ExprSet,as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
    
    entrezy_nazwy=unlist(mget(ExprSet,as.environment(as.list(gahgu95av2GENENAME)),ifnotfound=NA))
    
    
    macierz=data.frame(unlist(ExprSet),unlist(symbol), unlist(genNames),unlist(entrezy), unlist(entrezy_nazwy))
    names(macierz)[1]="próbka"
    names(macierz)[2]="symbol"
    names(macierz)[3]="nazwa"
    names(macierz)[4]="entrez_id"
    names(macierz)[5]="entrez_nazwa"
    
    ekspresje=matrix(0L, nrow=length(sondy),ncol=1) 
    
    ExprSet= new("ExpressionSet",expr=ekspresje, annotation="gahgu95av2.db",featureData=AnnotatedDataFrame(macierz))
    
    return(ExprSet)
    
  } else {
    
    ekspresje=exprs(ExprSet)
    
    symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) 
    
    genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME))
    
    entrezy=unlist(mget(featureNames(ExprSet),as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
    
    entrezy_nazwy=unlist(mget(featureNames(ExprSet),as.environment(as.list(gahgu95av2GENENAME)),ifnotfound=NA))
    
    
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
}

#wywołanie
up_ExprSet=updated_ExprSet(ExprSet)

#eksport informacji o genach do pliku Excela
zapis_xlsx=function(up_ExprSet){
  
  require(openxlsx)
  require(tcltk)
  
  #eksport informacji o genach do pliku Excela - wymaga intalacji Rtools
  #installr::install.rtools() #R.3.3.x or later
  #wskazanie pliku Rtools zip.exe, np. "C:/Rtools/bin/zip.exe"
  Sys.setenv("R_ZIPCMD" = tk_choose.files(caption = "Wskaż plik zip.exe w folderze Rtools>bin"))
  cechy=ExprSet@featureData
  cechy=cechy@data
  write.xlsx(cechy, "dane_geny.xlsx", asTable = TRUE)
}

#wywołanie
sondy=zapis_xlsx(up_ExprSet)

summary_table=function(ExprSet, class1, class2, correction_method, FoldChangeCutoff, sort_criterion, number, sep){
  #analiza cech różnicujących - test t i Fold Change
  
  require(Biobase)
  require(affy)
  require('gahgu95av2cdf')
  require('gahgu95av2.db')
  require(limma)
  
  adeno=which(pData(ExprSet)$CLASS==class1) 
  squamous=which(pData(ExprSet)$CLASS==class2) 
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
  
  #histogrm wszystkich p-wartości po korekcji zapisywany do pliku png
  png("histogram.png") # tworzymy rysunek histogramu
  p=as.double(TAB[,11])
  hist(p, main="Histogram skorygowanych p-wartości", xlab = "p-wartość",col="lightblue")
  dev.off() #zamykamy okno
  
  #wybor sond do tabeli wynikowej 
  #sortowanie po skorygowanych p-wartościach
  #jezeli sort_criterion to ulamek to poziom ufnosci p-wartosci, jeśli nie to komunikat z błędem np 0.05
  #jezeli liczba zwracanych wyników (number) >1 to ilosc sond do zwrócenia, jezeli =0 to wszystkie
  # jezeli FoldChangeCutoff <1 - zmniejszona ekspresja w stosunku do kontroli, >1 nadekspresja, =1 brak zmian 
  #użytkownik podaje liczbę >=1, np.2 (2 razy większa lub mniejsza ekspresja w stosunku do kontroli) 
  
  if (FoldChangeCutoff>=1){
    if (sort_criterion>1 | sort_criterion<=0)
    {
      require(tcltk)
      msgBox <- tkmessageBox(title = "Uwaga!",
                             message = "Błędna wartość zadanego poziomu ufności! Proszę podać liczbę z zakresu 0 do 1.",icon = "info", type = "ok")
    }
    
    if (sort_criterion<=1 & sort_criterion>0)
    {
      ind_sort=which(p_val_skorygowane<sort_criterion)
      TAB=TAB[ind_sort,]
    }
    
    if (sort_criterion==1)
    {
      ind_sort=sort(p_val_skorygowane,index=TRUE)$ix
      TAB=TAB[ind_sort,]
    }
    #sprawdzenie, czy wlasciwe sortowanie oraz wybór n sond z zadanym FoldChange
    if (number>0){
      ind_sort=sort(as.double(TAB[,11]),index=TRUE)$ix
      TAB=TAB[ind_sort[1:number],]
      ind_FC=which(as.double(TAB[,6])>FoldChangeCutoff | as.double(TAB[,6])<(FoldChangeCutoff)^(-1))
      TAB=TAB[ind_FC,]
      
    }
    if (number==0){
      ind_sort=sort(as.double(TAB[,11]),index=TRUE)$ix
      TAB=TAB[ind_sort,]
      ind_FC=which(as.double(TAB[,6])>FoldChangeCutoff | as.double(TAB[,6])<(FoldChangeCutoff)^(-1))
      TAB=TAB[ind_FC,]
    }
    if (number<0){
      require(tcltk)
      msgBox <- tkmessageBox(title = "Uwaga!",
                             message = "Błędna liczba rząanych wyników! Proszę podać liczbę większą od 0 lub wybrać 0 w celu wyświetlenia wszystkich wyników.",icon = "info", type = "ok")
    }
  } else {
    require(tcltk)
    msgBox <- tkmessageBox(title = "Uwaga!",
                           message = "Błędny poziom odcięcia Fold Change! Proszę podać liczbę większą lub równą 1.",icon = "info", type = "ok")
  }
  
  # zapisanie tablicy wynikowej
  write.table(TAB, file="summary_table.txt",sep=sep,row.names = FALSE) #separator pobierany jako paramer
  return(TAB)
} 

#przykładowe wywołanie - użytkownik wybiera:
#class1, class2 - użytkownik musi znać oznaczenia w phenodata$CLASS dla grup które chce porównywać - 2 pola do wpiania przez użytkownika
#correction_method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none") - wybór przez użytkownika
#FoldChangeCutoff: użytkownik podaje liczbę >=1, np.2 (2 razy większa lub mniejsza ekspresja w stosunku do kontroli)
#sort_criterion: użytkownik podaje poziom ufności dla p-wartości
#number - liczba wyników do wyświetlenia zadna przez użytkownika (gdy =0 to wyświetlamy wszystkie)
#sep - separtor dla pliku txt z wynikami - pole do wpisania przez użytkownika
TAB_geny_roznicujace=summary_table(ExprSet, class1='ADENO', class2='SQUAMOUS', correction_method='holm', FoldChangeCutoff=1.2 ,sort_criterion=0.01, number=0, sep=', ')

#tworzenie heatmapy
heatmapa=function(ExprSet, TAB_geny_roznicujace){
  
  require(Biobase)
  require(gplots)
  
  #wartosci ekspresji dla wybranych do tabeli gen?w
  expr=exprs(ExprSet)
  ind=TAB_geny_roznicujace[,1]
  ekspr_wybrane=expr[as.double(ind),] # tablica ekspresji tylko dla wybranych zestaw?w sond
  nr=dim(ekspr_wybrane)[1]
  nc=dim(ekspr_wybrane)[2]
  png("heatmap.png", width=1000, height=800) # tworzymy rysunek mapy ciep?a 
  # rozbudowana wersja mapy ciep?a
  #oś x to obserwacje, a oś y to sondy
  par(oma=c(1,1,1,1))
  par(mar=c(2,2,2,2))
  # definiowanie własnego układu heatmapy
  # tabela 3x3 ze zdefiniowaną lokalizacją elementów heatmapy
  mylmat = rbind(c(0,3,0), #tytuł
                 c(0,1,2), #dendrogramy
                 c(0,4,0)) #legenda kolorów
  mylwid = c(0.2,5,0.1) #szerokość każdej z trzech sekcji kolumnowych plota
  mylhei = c(0.4,5,0.8) #wysokość każdej z trzech sekcji kolumnowych plota
  heatmap.2(ekspr_wybrane, 
            main = "Ekspresja genów różnicujących", # heat map title
            xlab = "obserwacje",
            ylab = "sondy",
            col = bluered, #paleta kolorów
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            lmat=mylmat, lwid=mylwid, lhei=mylhei,
            margins =c(8,8),     # widens margins around plot 
            dendrogram="none",     #dendrogram ukazuje związki między elementami na podstawie przyjętego kryterium
            Colv="NA",
            cexRow = 0.07 + 1/log10(nr), #rozmiary czcionki nazw wierszy i kolumn
            cexCol = 0.07 + 1/log10(nc),
            offsetRow = 0.2,    #odległość nazw obserwacji od heatmapy
            offsetCol = 0.2)    #odległość nazw sond od heatmapy
  dev.off() #zamykamy okno
}

#wywołanie
heat_mapa=heatmapa(ExprSet, TAB_geny_roznicujace)

------------------------------------------------------------------------------------------------------


#wczytanie wlasnej listy z genami do analizy sciezek sygnalnych - własna lista lub zapisany plik 
wczytaj_liste=function(){
  
  require(tcltk)
  
  #wczytanie pliku tekstowego z nazwami sond wybranego przez użytkownika
  lista_genow=lapply(read.delim(tk_choose.files(caption = "Wybierz plik txt. z listą genów do analizy scieżek sygnalnych"), header = FALSE, stringsAsFactor = FALSE), as.character)$V1
  
  return(lista_genow)
}

#wywołanie
lista_genow=wczytaj_liste()


#wszystkie oraz różnicujące geny do analizy ścieżek sygnalnych np n podstawie entrez id
wszystkie_geny=ExprSet@featureData
wszystkie_geny=wszystkie_geny@data
wszystkie_geny_entrez_id=wszystkie_geny$entrez_id

roznicujce_geny_entrez_id=TAB_geny_roznicujace[,4]
