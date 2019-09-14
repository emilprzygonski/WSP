heatmapa=function(ExprSet, TAB_geny_roznicujace){
  
  require(Biobase)
  require(gplots)
  
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
            main = "Ekspresja genow roznicujacych", # heat map title
            xlab = "obserwacje",
            ylab = "sondy",
            col = bluered, #paleta kolorow
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            lmat=mylmat, lwid=mylwid, lhei=mylhei,
            margins =c(8,8),     # widens margins around plot 
            dendrogram="none",     #dendrogram ukazuje zwiazki miedzy elementami na podstawie przyjetego kryterium
            Colv="NA",
            cexRow = 0.07 + 1/log10(nr), #rozmiary czcionki nazw wierszy i kolumn
            cexCol = 0.07 + 1/log10(nc),
            offsetRow = 0.2,    #odleglosc nazw obserwacji od heatmapy
            offsetCol = 0.2)    #odleglosc nazw sond od heatmapy
  dev.off() #zamykamy okno
}

#wywolanie
heat_mapa=heatmapa(ExprSet, TAB_geny_roznicujace)