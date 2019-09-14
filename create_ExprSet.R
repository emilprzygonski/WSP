create_ExprSet=function(){
  
  #funkcja tworzaca obiekt typu ExpressionSet na podstawie plikow .CEL, 
  #ktorych nazwy i dane znajduja sie w pliku z adnotacjami. Plik ten powinien zostac wskazany przez uzytkownika.
  
  require(Biobase)
  require(affy)
  require(tcltk)
  require('gahgu95av2.db')
  
  #wybor folderu z plikami .CEL (dziala tylko gdy pliki znajduja sie w folderze z ta funkcja)
  choose.dir(getwd(), "Wybierz folder z plikami .CEL")
  #wybor pliku z adnotacjami
  file=tk_choose.files(caption = "Wybierz plik .txt z adnotacjami")
  #wczytnie danych
  data = read.table(file, header = TRUE, sep = "\t")
  opisy = read.AnnotatedDataFrame(file, sep="\t", header=TRUE, row.names=4,stringsAsFactors = F) 
  experiment = new("MIAME", name = "Dane mikromacierzowe",lab = "IO",title = "dane mikromacierzowe",url = "http://www.bioconductor.org", other = list(notes = "inne"))
  
  #nazwy plikow.CEL
  sampleNames(opisy) = paste(sampleNames(opisy), ".CEL", sep="")
  #wczytanie plikow i tworzenie macierzy ekspresji
  data_Affy=ReadAffy(filenames=sampleNames(opisy), verbose=TRUE)
  data_Affy@cdfName=paste("ga",data_Affy@cdfName,sep="")
  data_Affy@annotation=paste("ga",data_Affy@annotation,sep="")
  
  #normalizacja 
  RMA=rma(data_Affy)
  dataRMA=exprs(RMA) 
  
  #tworzenie nowego obiektu ExpressionSet
  ExprSet= new("ExpressionSet", expr=dataRMA, phenoData = opisy,experimentData=experiment,annotation="gahgu95av2.db") 
  x= exprs(ExprSet)
  write.csv(x, "ekspresje.csv", col.names=TRUE, row.names=TRUE)
  
  return(ExprSet)
}

#wywolanie
ExprSet=create_ExprSet() 