zapis_xlsx=function(up_ExprSet){
  
  require(openxlsx)
  require(tcltk)
  
  #eksport informacji o genach do pliku Excela - wymaga intalacji Rtools
  #installr::install.rtools() #R.3.3.x or later
  #wskazanie pliku Rtools zip.exe, np. "C:/Rtools/bin/zip.exe"
  #Sys.setenv("R_ZIPCMD" = tk_choose.files(caption = "Wskaz plik zip.exe w folderze Rtools>bin"))
  cechy=up_ExprSet@featureData
  cechy=cechy@data
  write.xlsx(cechy, "dane_geny.xlsx", asTable = TRUE)
}

#wywolanie
sondy=zapis_xlsx(up_ExprSet)