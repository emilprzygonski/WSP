#wczytanie Expressionset (tworzenie z gotowej macierzy ekspresji)
create_ExprSet2=function(){
  
  require(Biobase)
  require(affy)
  require(tcltk)
  require('gahgu95av2.db')
  
  
  exprsFile = read.csv(tk_choose.files(caption = "Wybierz plik .txt z macierzą ekspresji"), row.names = 1, header= TRUE)
  s<-as.matrix(exprsFile)
  #wybor pliku z adnotacjami
  adnotFile = tk_choose.files(caption = "Wybierz plik .txt z adnotacjami")
  #wczytnie danych
  phenoFile = read.AnnotatedDataFrame(adnotFile, sep="\t", header=TRUE, row.names=4,stringsAsFactors = F)
  experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))
  
  #nazwy plików.CEL
  sampleNames(phenoFile) = paste(sampleNames(phenoFile), ".CEL", sep="")
  #tworzenie nowego obiektu ExpressionSet
  ExprSet= new("ExpressionSet", expr=s, phenoData = phenoFile,experimentData=experiment,annotation="gahgu95av2.db") 
  
  return(ExprSet)
}

ExprSet=create_ExprSet2()