przeprowadzanie_GSA=function(wybor_genow, up_ExprSet, TAB_geny_roznicujace){
  #jako wejœcie funkcja przyjmuje:
  #wybor_genow - wszystkie z ExprSet (1) lub tylko roznicujace (2) (wybor opcji)
  #analizowany i zaktualizowany objekt ExprSet
  #tabele genow roznicujacych
  
  require(GSA)
  require(tcltk)
  
  #wczytanie pliku tekstowego z list¹ œciezek
  entrez_ids=tk_choose.files(caption = "Proszê wybraæ plik .gmt zawieraj¹cy listê œcie¿ek z entrez id dla genów")
  entrezid=GSA.read.gmt(entrez_ids)
  
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

GSAwynik=przeprowadzanie_GSA(1, up_ExprSet, TAB_geny_roznicujace)

