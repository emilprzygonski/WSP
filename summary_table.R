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
  pval=apply(expr,1,function(x) t.test(x[adeno],x[squamous])$p.val) # p wartosci
  #Korekta na wielokrotne testowanie
  p_val_skorygowane=p.adjust(pval, method = correction_method) # metoda korekcji pobierana jako parametr funkcji
  #SYMBOLE GENOW
  symbol=unlist(mget(featureNames(ExprSet),env=gahgu95av2SYMBOL)) 
  #NAZWY GENOW
  genNames=unlist(mget(featureNames(ExprSet),env=gahgu95av2GENENAME)) 
  
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
  #histogram wszystkich p-wartości po korekcji zapisywany do pliku png
  png("histogram.png") # tworzymy rysunek histogramu
  p=as.double(TAB[,11])
  hist(p, main="Histogram skorygowanych p-wartosci", xlab = "p-wartosc",col="lightblue")
  dev.off() #zamykamy okno
  
  #wybor sond do tabeli wynikowej 
  #sortowanie po skorygowanych p-wartościach
  #jezeli sort_criterion to ulamek to poziom ufnosci p-wartosci, jesli nie to komunikat z bledem np 0.05
  #jezeli liczba zwracanych wynikow (number) >1 to ilosc sond do zwrocenia, jezeli =0 to wszystkie
  # jezeli FoldChangeCutoff <1 - zmniejszona ekspresja w stosunku do kontroli, >1 nadekspresja, =1 brak zmian 
  #uzytkownik podaje liczbe >=1, np.2 (2 razy wieksza lub mniejsza ekspresja w stosunku do kontroli) 
  
  if (FoldChangeCutoff>=1){
    if (sort_criterion>1 | sort_criterion<=0)
    {
      require(tcltk)
      msgBox <- tkmessageBox(title = "Uwaga!",
                             message = "Bledna wartosc zadanego poziomu ufnosci! Prosze podac liczbe z zakresu 0 do 1.",icon = "info", type = "ok")
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
                             message = "Bledna liczba rzaanych wynikow! Prosze podac liczbe wieksza od 0 lub wybrac 0 w celu wyswietlenia wszystkich wynikow.",icon = "info", type = "ok")
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
  
  # zapisanie tablicy wynikowej
  write.table(TAB, file="summary_table.txt",sep=sep,row.names = FALSE) #separator pobierany jako paramer
  return(TAB)
} 

#przykladowe wywolanie - uzytkownik wybiera:
#class1, class2 - uzytkownik musi znac oznaczenia w phenodata$CLASS dla grup ktore chce porownywac - 2 pola do wpiania przez uzytkownika
#correction_method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none") - wybor przez uzytkownika
#FoldChangeCutoff: uzytkownik podaje liczbe >=1, np.2 (2 razy wieksza lub mniejsza ekspresja w stosunku do kontroli)
#sort_criterion: uzytkownik podaje poziom ufnosci dla p-wartosci
#number - liczba wynikow do wyswietlenia zadna przez uzytkownika (gdy =0 to wyswietlamy wszystkie)
#sep - separtor dla pliku txt z wynikami - pole do wpisania przez uzytkownika
TAB_geny_roznicujace=summary_table(ExprSet, class1='ADENO', class2='SQUAMOUS', correction_method='holm', FoldChangeCutoff=1.2 ,sort_criterion=0.01, number=0, sep=', ')
