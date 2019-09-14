updated_ExprSet=function(ExprSet){
  #funkcja dodajaca symbole genow, ich nazwy oraz entrez id do ExprSet
  
  require(Biobase)
  require(affy)
  require('gahgu95av2cdf')
  require('gahgu95av2.db')
  require('org.Hs.eg.db')
  
  if (is.character(ExprSet)=='TRUE'){ #gdy uzytkownik wczytuje jedynie lite z nazwami sond - zwracany ExprSet z zerowymi ekspresjami
    
    symbol=unlist(mget(ExprSet,env=gahgu95av2SYMBOL)) 
    
    genNames=unlist(mget(ExprSet,env=gahgu95av2GENENAME))
    
    entrezy=unlist(mget(ExprSet,as.environment(as.list(gahgu95av2ENTREZID)),ifnotfound=NA))
    
    entrezy_nazwy=unlist(mget(ExprSet,as.environment(as.list(gahgu95av2GENENAME)),ifnotfound=NA))
    
    
    macierz=data.frame(unlist(ExprSet),unlist(symbol), unlist(genNames),unlist(entrezy), unlist(entrezy_nazwy))
    names(macierz)[1]="probka"
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
}

#wywolanie
up_ExprSet=updated_ExprSet(ExprSet)