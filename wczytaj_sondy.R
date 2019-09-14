wczytaj_sondy=function(){
  
  require(tcltk)
  
  #wczytanie pliku tekstowego z nazwami sond wybranego przez uzytkownika
  sondy=lapply(read.delim(tk_choose.files(caption = "Wybierz plik .txt z nazwami sond"), header = FALSE, stringsAsFactor = FALSE), as.character)$V1
  
  return(sondy)
}

#wywolanie
sondy=wczytaj_sondy()