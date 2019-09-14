#wczytanie Expressionset (tworzenie z gotowej macierzy ekspresji)
create_ExprSet2=function(csv_file,txt_file){
  
  require(Biobase)
  require(affy)
  require(tcltk)
  require('gahgu95av2.db')
  

  s<-as.matrix(csv_file)
  #wczytnie danych
  phenoFile = read.AnnotatedDataFrame(txt_file, sep="\t", header=TRUE, row.names=4,stringsAsFactors = F)
  experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))
  
  #nazwy plików.CEL
  sampleNames(phenoFile) = paste(sampleNames(phenoFile), ".CEL", sep="")
  #tworzenie nowego obiektu ExpressionSet
  ExprSet= new("ExpressionSet", expr=s, phenoData = phenoFile,experimentData=experiment,annotation="gahgu95av2.db") 
  
  return(ExprSet)
}

updated_ExprSet=function(ExprSet){
  #funkcja dodajaca symbole genow, ich nazwy oraz entrez id do ExprSet
  
  require(Biobase)
  require(affy)
  require('gahgu95av2cdf')
  require('gahgu95av2.db')
  require('org.Hs.eg.db')
  
  require(openxlsx)
  require(tcltk)
  

    ekspresje=exprs(ExprSet)
    
    symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) 
    
    genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME))
    
    entrezy=unlist(mget(featureNames(ExprSet),as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
    
    entrezy_nazwy=unlist(mget(featureNames(ExprSet),as.environment(as.list(gahgu95av2GENENAME)),ifnotfound=NA))
    
    
    macierz=data.frame(unlist(featureNames(ExprSet)),unlist(symbol), unlist(genNames),unlist(entrezy), unlist(entrezy_nazwy))
    names(macierz)[1]="probka"
    names(macierz)[2]="symbol"
    names(macierz)[3]="nazwa"
    names(macierz)[4]="entrez_id"
    names(macierz)[5]="entrez_nazwa"
    
    opisy=ExprSet@phenoData
    experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))
    ExprSet= new("ExpressionSet", expr=ekspresje, phenoData = opisy,experimentData=experiment,annotation="gahgu95av2.db",featureData=AnnotatedDataFrame(macierz))

    cechy=ExprSet@featureData
    cechy=cechy@data
    return(cechy)
  
}

dla_gsa<-function(ExprSet){
  #funkcja dodajaca symbole genow, ich nazwy oraz entrez id do ExprSet
  
  require(Biobase)
  require(affy)
  require('gahgu95av2cdf')
  require('gahgu95av2.db')
  require('org.Hs.eg.db')
  
  require(openxlsx)
  require(tcltk)
  
  
  ekspresje=exprs(ExprSet)
  
  symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) 
  
  genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME))
  
  entrezy=unlist(mget(featureNames(ExprSet),as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
  
  entrezy_nazwy=unlist(mget(featureNames(ExprSet),as.environment(as.list(gahgu95av2GENENAME)),ifnotfound=NA))
  
  
  macierz=data.frame(unlist(featureNames(ExprSet)),unlist(symbol), unlist(genNames),unlist(entrezy), unlist(entrezy_nazwy))
  names(macierz)[1]="probka"
  names(macierz)[2]="symbol"
  names(macierz)[3]="nazwa"
  names(macierz)[4]="entrez_id"
  names(macierz)[5]="entrez_nazwa"
  
  opisy=ExprSet@phenoData
  experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))
  ExprSet= new("ExpressionSet", expr=ekspresje, phenoData = opisy,experimentData=experiment,annotation="gahgu95av2.db",featureData=AnnotatedDataFrame(macierz))
  return(ExprSet)
}

dla_sond<-function(sondy){
  require(Biobase)
  require(affy)
  require('gahgu95av2cdf')
  require('gahgu95av2.db')
  require('org.Hs.eg.db')
  
  require(openxlsx)
  require(tcltk)
  
  symbol=unlist(mget(as.character(sondy$V1),env=gahgu95av2SYMBOL)) 
  
  genNames=unlist(mget(as.character(sondy$V1),env=gahgu95av2GENENAME))
  
  entrezy=unlist(mget(as.character(sondy$V1),as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
  
  entrezy_nazwy=unlist(mget(as.character(sondy$V1),as.environment(as.list(gahgu95av2GENENAME)),ifnotfound=NA))
  
  
  cechy=data.frame(unlist(sondy),unlist(symbol), unlist(genNames),unlist(entrezy), unlist(entrezy_nazwy))
  names(cechy)[1]="probka"
  names(cechy)[2]="symbol"
  names(cechy)[3]="nazwa"
  names(cechy)[4]="entrez_id"
  names(cechy)[5]="entrez_nazwa"
  
  # ekspresje=matrix(0L, nrow=length(sondy),ncol=1) 
  # 
  # ExprSet= new("ExpressionSet",expr=ekspresje, annotation="gahgu95av2.db",featureData=AnnotatedDataFrame(macierz))
  # 
  # cechy=ExprSet@featureData
  # cechy=cechy@data
  
  return(cechy)
  
}

summary_table=function(ExprSet, class1, class2, correction_method, FoldChangeCutoff, sort_criterion, number, sep){
  #analiza cech roznicujacych - test t i Fold Change
  
  require(Biobase)
  require(affy)
  require('gahgu95av2cdf')
  require('gahgu95av2.db')
  require(limma)
  
  adeno=which(pData(ExprSet)$CLASS==class1) 
  squamous=which(pData(ExprSet)$CLASS==class2) 
  expr=exprs(ExprSet)
  klasy_wektor<-array(dim=c(dim(expr)[1],1))
  sr_adeno=rowMeans(expr[,adeno])
  sr_squamous=rowMeans(expr[,squamous])
  # fold change
  FC=sr_adeno/sr_squamous
  #test t
  statistic=apply(expr,1,function(x) t.test(x[adeno],x[squamous])$statistic)
  #przydatne polecenia apply vapply sapply..
  pval=apply(expr,1,function(x) t.test(x[adeno],x[squamous])$p.val) # p wartosci
  #Korekta na wielokrotne testowanie
  p_val_skorygowane=p.adjust(pval, method = correction_method) # metoda korekcji pobierana jako parametr funkcji
  #SYMBOLE GENOW
  symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) #mozna jako lista albo wektor
  #NAZWY GENOW
  genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME)) # podobnie lista przeksztalcona na wektor
  
  entrezy=unlist(mget(rownames(expr),as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
  
  TAB=array(dim=c(dim(expr)[1],12))
  colnames(TAB)=c("Indeks","FerrariID","Symbol","description","entrez id","fold change","srednia w gr.ADENO","srednia w gr.NORMAL","t-statistics","pvalue","corrected p-value","klasy")
  #tabela dla wszystkich sond, przed selekcja mozna juz wczesniej zainicjowac i bezposrednio wpisywac
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
  TAB[,12]=klasy_wektor
  
  #histogrm wszystkich p-wartosci po korekcji zapisywany do pliku png
  png("histogram.png") # tworzymy rysunek histogramu
  p=as.double(TAB[,11])
  hist(p, main="Histogram skorygowanych p-wartości", xlab = "p-wartość",col="lightblue")
  dev.off() #zamykamy okno
  
  #wybor sond do tabeli wynikowej 
  #sortowanie po skorygowanych p-wartosciach
  #jezeli sort_criterion to ulamek to poziom ufnosci p-wartosci, jesli nie to komunikat z blędem np 0.05
  #jezeli liczba zwracanych wynikow (number) >1 to ilosc sond do zwrocenia, jezeli =0 to wszystkie
  # jezeli FoldChangeCutoff <1 - zmniejszona ekspresja w stosunku do kontroli, >1 nadekspresja, =1 brak zmian 
  #uzytkownik podaje liczbę >=1, np.2 (2 razy wieksza lub mniejsza ekspresja w stosunku do kontroli) 
  
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
    #sprawdzenie, czy wlasciwe sortowanie oraz wybor n sond z zadanym FoldChange
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
                           message = "Bledny poziom odciecia Fold Change! Prosze podac liczbe wieksza lub rowna 1.",icon = "info", type = "ok")
  }
  
  klasy_wektor<-array(dim=c(dim(TAB)[1],1))
  klasy_wektor[1]<-class1
  klasy_wektor[2]<-class2
  TAB[,12]=klasy_wektor
  
  # # zapisanie tablicy wynikowej
  # write.table(TAB, file="summary_table.txt",sep=sep,row.names = FALSE) #separator pobierany jako paramer
  return(TAB)
} 


heatmapa=function(ExprSet, TAB_geny_roznicujace){
  
  require(Biobase)
  require(gplots)
  
  #tworzenie heatmapy
  #wartosci ekspresji dla wybranych do tabeli genow
  expr=exprs(ExprSet)
  ind=TAB_geny_roznicujace[,1]
  ekspr_wybrane=expr[as.double(ind),] # tablica ekspresji tylko dla wybranych zestawow sond
  nr=dim(ekspr_wybrane)[1]
  nc=dim(ekspr_wybrane)[2]
  png("heatmap.png", width=1000, height=800) # tworzymy rysunek mapy ciepla 
  # rozbudowana wersja mapy ciepla
  #os x to obserwacje, a os y to sondy
  par(oma=c(1,1,1,1))
  par(mar=c(2,2,2,2))
  # definiowanie wlasnego ukladu heatmapy
  # tabela 3x3 ze zdefiniowana lokalizacja elementow heatmapy
  mylmat = rbind(c(0,3,0), #tytul
                 c(0,1,2), #dendrogramy
                 c(0,4,0)) #legenda kolorow
  mylwid = c(0.2,5,0.1) #szerokosc kazdej z trzech sekcji kolumnowych plota
  mylhei = c(0.4,5,0.8) #wysokosc kazdej z trzech sekcji kolumnowych plota
  heatmap.2(ekspr_wybrane, 
            main = "Ekspresja genów różnicujących", # heat map title
            xlab = "obserwacje",
            ylab = "sondy",
            col = bluered, #paleta kolorow
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            lmat=mylmat, lwid=mylwid, lhei=mylhei,
            margins =c(8,8),     # widens margins around plot 
            dendrogram="none",     #dendrogram ukazuje zwiazki między elementami na podstawie przyjetego kryterium
            Colv="NA",
            cexRow = 0.07 + 1/log10(nr), #rozmiary czcionki nazw wierszy i kolumn
            cexCol = 0.07 + 1/log10(nc),
            offsetRow = 0.2,    #odległosc nazw obserwacji od heatmapy
            offsetCol = 0.2)    #odległosc nazw sond od heatmapy
  dev.off() #zamykamy okno
}

MSigDBr_import<-function(gs_category){
  #funkcja pobierajaca z bazy MSigDB geneset z wybranej przez uzytkownika kategorii oraz zapisujaca
  #przetworzony geneset do pliku .gmt
  
  require(msigdbr)
  require(GSEABase)
  
  if (gs_category==1){
    Gene_Sets=msigdbr(species='Homo sapiens', category='H')
  }
  else if (gs_category==2){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C1')
  }
  else if (gs_category==3){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C2')
  }
  else if (gs_category==4){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C2', subcategory='CGP')
  }
  else if (gs_category==5){
    Gene_Sets=msigdbr(species='Homo sapiens',category='C2', subcategory='CP')
  }
  else if (gs_category==6){
    Gene_Sets=msigdbr(species='Homo sapiens',category='C2', subcategory='CP:BIOCARTA')
  }
  else if (gs_category==7){
    Gene_Sets=msigdbr(species='Homo sapiens',category='C2', subcategory='CP:KEGG')
  }
  else if (gs_category==8){
    Gene_Sets=msigdbr(species='Homo sapiens',category='C2', subcategory='CP:REACTOME')
  }
  else if (gs_category==9){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C3')
  }
  else if (gs_category==10){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C3', subcategory='MIR')
  }
  else if (gs_category==11){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C3', subcategory='TFT')
  }
  else if (gs_category==12){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C4')
  }
  else if (gs_category==13){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C4', subcategory='CGN')
  }
  else if (gs_category==14){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C4', subcategory='CM')
  }
  else if (gs_category==15){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C5')
  }
  else if (gs_category==16){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C5', subcategory='BP')
  }
  else if (gs_category==17){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C5', subcategory='CC')
  }
  else if (gs_category==18){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C5', subcategory='MF')
  }
  else if (gs_category==19){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C6')
  }
  else if (gs_category==20){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C7')
  }
  
  geneset.names=unique(Gene_Sets$gs_name)
  # genesets<-vector("list",length(geneset.names))
  # for(i in 1:length(geneset.names)){
  #   nazwa<-geneset.names[i]
  #   geny<-Gene_Sets$gene_symbol[Gene_Sets$gs_name==nazwa]
  #   genesets[[i]] <- GeneSet(setName=nazwa, geneIds=as.character(geny), geneIdType=SymbolIdentifier())
  # }
  # 
  # Names_Col=GeneSetCollection(genesets)
  # toGmt(Names_Col, "GeneSetNames.gmt")
  
  genesets_e<-vector("list",length(geneset.names))
  for(i in 1:length(geneset.names)){
    nazwa<-geneset.names[i]
    geny<-Gene_Sets$entrez_gene[Gene_Sets$gs_name==nazwa]
    genesets_e[[i]] <- GeneSet(setName=nazwa, geneIds=as.character(geny), geneIdType=EntrezIdentifier())
  }
  
  Entrez_Col=GeneSetCollection(genesets_e)
  toGmt(Entrez_Col, "EntrezIDs.gmt")
}


przeprowadzanie_GSA=function(wybor_genow, up_ExprSet, TAB_geny_roznicujace,entrezid){
  #jako wej�cie funkcja przyjmuje:
  #wybor_genow - wszystkie z ExprSet (1) lub tylko roznicujace (2) (wybor opcji)
  #analizowany i zaktualizowany objekt ExprSet
  #tabele genow roznicujacych
  
  require(GSA)
  require(tcltk)
  
  #wczytanie pliku tekstowego z list� �ciezek- to w serwer.R
  
  if (wybor_genow == 1){
    wszystkie_geny=up_ExprSet@featureData
    wszystkie_geny=wszystkie_geny@data
    entrez_names=wszystkie_geny$entrez_id
    entrez_names <- entrez_names[!is.na(entrez_names)]
    entrez_names <- as.character(entrez_names)
    
    ekspresje=exprs(up_ExprSet)  
    
    klasy=pData(up_ExprSet)$CLASS
    uniq_klasy=unique(klasy)
    for (i in 1:length(uniq_klasy)){
      klasy[klasy==uniq_klasy[i]]=i
    }
    
    genenames<-wszystkie_geny$symbol
  }
  
  if (wybor_genow == 2) {
    entrez_names=TAB_geny_roznicujace[,5]
    
    ekspresje=exprs(up_ExprSet)
    indeksy=as.numeric(TAB_geny_roznicujace[,1])
    ekspresje=ekspresje[indeksy,]
    
    klasy=TAB_geny_roznicujace[1:2,12]
    cl1=which(pData(up_ExprSet)$CLASS==klasy[1]) 
    cl2=which(pData(up_ExprSet)$CLASS==klasy[2]) 
    klasy=as.vector(array(dim=c(max(cl1,cl2),1)))
    klasy[cl1]=1
    klasy[cl2]=2
    klasy <- klasy[!is.na(klasy)]
    genenames<-TAB_geny_roznicujace[,3]
  }
  
  
  x<-unname(as.array(ekspresje))
  y<-klasy
  genesety<-entrezid$genesets
  genesetsnames<-entrezid[["geneset.names"]]
  
  GSA.obj<-GSA(x,y, genenames=entrez_names, genesets=genesety,  resp.type="Multiclass", nperms=100)
  
  wynik<-GSA.listsets(GSA.obj, geneset.names = genesetsnames, FDRcut = 1)
  wynik<-wynik[["positive"]]
  
  GSA.plot(GSA.obj)
  return(wynik)
}



