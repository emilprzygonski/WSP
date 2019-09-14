options(shiny.maxRequestSize=40*1024^2)
#install.packages("shiny", type="binary")
library('shiny')
library(ggplot2)
library(shinyTime)
library(log4r)
library('htmlwidgets')
library(shinythemes)



require(tcltk)
require(openxlsx)


loggerDebug <- create.logger()
logfile(loggerDebug) <- 'data/debugData.log'
level(loggerDebug) <- 'DEBUG'

loggerServer <- create.logger()
logfile(loggerServer) <- 'data/serverData.log'
level(loggerServer) <- 'DEBUG'


source("funkcje/ExprSet.R",local = TRUE, echo = FALSE)
shinyServer(function(input, output, session) {
  
  #UPLOAD FILES = inputs
  
  csv <-reactive({
    inFile <- input$input1
    if (is.null(inFile))
      return(NULL)
    df <- read.csv(inFile$datapath, row.names=1, header = TRUE, sep=",")
    return(df)
  })
    
  output$plik1<- renderTable({
    df <- csv()
    if(is.null(df))
      return(NULL)
    return(df)
  })
  
  txt <-reactive({
    inFile2 <- input$input2$datapath
    if(is.null(inFile2))
      return(NULL)
    return(inFile2)
  })
  
  
  txt2 <-reactive({
    inFile3 <- input$input2
    if (is.null(inFile3))
      return(NULL)
    df <- read.table(inFile3$datapath, header=TRUE, fill=TRUE)
    return(df)
  })
  
  sondy<-reactive({
    inFile4<-input$input3
    if(is.null(inFile4))
      return(NULL)
    df<-read.delim(inFile4$datapath, header = FALSE, stringsAsFactor = FALSE)
    return(df)
  })
  
  #inaczej wczytanie - read.table() zamiast lapply(read.delim())
  sondy2<-reactive({
    inFile4<-input$input3
    if(is.null(inFile4))
      return(NULL)
    df<- read.table(inFile4$datapath)
    return(df)
  })
  
  gsa_sciezki<-reactive({
    File<-input$input_gsa$datapath
    if(is.null(File))
      return(NULL)
    return(File)
  })
  
  output$tabela_sondy<-renderTable({
    df4<-sondy2()
    if(is.null(df4))
      return(NULL)
    return(df4)
  })
  
  output$plik2<- renderTable({
    df5 <- txt2()
    if(is.null(df5))
      return(NULL)
    return(df5)
  })
  
  output$gsa<-renderTable({
    wynik<-gsa_sciezki()
    if(is.null(gsa_sciezki()))
      return(NULL)
    return(wynik)
  })
  
  #WYNIKI - CZESC BIOINFORMATYCZNA 
 
  wynik_macierz<-reactive({
    if(is.null(csv()))
      return(NULL)
    if(is.null(txt()))
      return(NULL)
    create_ExprSet2(csv(),txt())
  })
  
  wynik_macierz2<-reactive({
    if(is.null(wynik_macierz()))
      return(NULL)
    updated_ExprSet(wynik_macierz())
  })
  
  wynik_dla_gsa<-reactive({
    if(is.null(wynik_macierz()))
      return(NULL)
    dla_gsa(wynik_macierz())
  })
  
  output$tabela_cechy_macierz<- renderTable({
    if(is.null(wynik_macierz2()))
      return(NULL)
    tabela <- wynik_macierz2()
    return(tabela)
  })
  
  wynik_sondy<-reactive({
    if(is.null(sondy2()))
      return(NULL)
    dla_sond(sondy2())
  })
  
  output$tabela_cechy_sondy<- renderTable({
    if(is.null(wynik_sondy()))
      return(NULL)
    tabela <- wynik_sondy()
    return(tabela)
  })
  
  tabela<-reactive({
    
    if(is.null(wynik_macierz()))
      return(NULL)
    
    class1<-req(input$class1)
    class2<-req(input$class2)
    correction_method<-req(input$correction_method)
    FoldChangeCutoff<-req(input$FoldChangeCutoff)
    sort_criterion<-req(input$sort_criterion)
    number<-req(input$number)
    sep<-req(input$sep)
    
    summary_table(wynik_macierz(),class1,class2,correction_method,FoldChangeCutoff,sort_criterion,number,sep)
    
  })
  
  output$summary_table<-renderTable({
    if(is.null(wynik_macierz()))
      return(NULL)
    
    # class1<-req(input$class1)
    # class2<-req(input$class2)
    # correction_method<-req(input$correction_method)
    # FoldChangeCutoff<-req(input$FoldChangeCutoff)
    # sort_criterion<-req(input$sort_criterion)
    # number<-req(input$number)
    # sep<-req(input$sep)
    
    #summary_table(wynik_macierz(),class1,class2,correction_method,FoldChangeCutoff,sort_criterion,number,sep)
    
    tabela <- tabela()
    return(tabela)
  })

  
  #utworzenie heatmapy po kliknieciu przycisku 'Utworz heatmape' i zapisanie w folderze projektu
  observeEvent(input$create_heatmap, {
    if(is.null(csv()))
      return(NULL)
    if(is.null(txt()))
      return(NULL)
    if(is.null(wynik_macierz()))
      return(NULL)
    if(is.null(tabela()))
      return(NULL)
    heatmapa<-heatmapa(wynik_macierz(),tabela())
  })
  
  
  

  #import sciezek sygnalowych ze strony 
  
  observeEvent(input$import_z, {
    if(is.null(input$kategoria))
      return(NULL)
    kategoria<-req(input$kategoria)
    import<-MSigDBr_import(kategoria)
  })
  
 
  #wykonanie analizy GSA
  observeEvent(input$wykonaj_gsa,{
    
    wybor_genow<-req(input$kategoria2)
    up_ExprSet<-wynik_dla_gsa()
    TAB_geny_roznicujace<-tabela()
    
    entrezid=GSA.read.gmt(gsa_sciezki())
    
    if(is.null(wybor_genow))
      return(p("Prosze wybrac plik z lista sciezek"))
    
    if(is.null(up_ExprSet))
      return(NULL)
    
    if(is.null(TAB_geny_roznicujace))
      return(NULL)
    
    require(GSA)
    
    
    gsa<-przeprowadzanie_GSA(wybor_genow,up_ExprSet,TAB_geny_roznicujace,entrezid)
    gsa<-reactive({
      gsa_wynik<-gsa
      return(gsa_wynik)
    })
  
  })
  
  output$wynik_analizy<-renderTable({
    wynik<-gsa()
    return(wynik)
  })
  
  
  #DOWNLOAD RESULTS
  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste(input$input1, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(csv(), file, row.names = FALSE)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste(input$input2, ".txt", sep = "")
    },
    content = function(file) {
      write.table(txt(), file, row.names = FALSE)
    }
  )
  
  #pobieranie wyniku dla macierzy ekspresji
  output$downloadData3 <- downloadHandler(
    
    # filename = function() {"dane_geny.xlsx"},
    # content = function(file) {
    #   fname <- paste(file,"xlsx",sep=".")
    #   wb <- loadWorkbook(fname,create = TRUE)
    #   createSheet(wb,"geny")
    #   writeWorksheet(wb,data=wynik_macierz2(),sheet = "geny")
    #   saveWorkbook(wb)
    #   file.rename(fname,file)
    # },
    # contentType="application/xlsx"
    filename = function() { "dane_geny.xlsx" },
    content = function(file) {
      
      tempFile <- tempfile(fileext = ".xlsx")
      write.xlsx(wynik_macierz2(), tempFile)
      file.rename(tempFile, file)
      file.copy(filename,file)
      
    }
  )
  
  
  #pobieranie wyniku dla sond
  # output$downloadData4 <- downloadHandler(
  #   
  #   filename = function() {"dane_geny_sondy.xlsx"},
  #   content = function(file) {
  #     fname <- paste(file,"xlsx",sep=".")
  #     wb <- loadWorkbook(fname,create = TRUE)
  #     createSheet(wb,"geny")
  #     writeWorksheet(wb,data=wynik_sondy(),sheet = "geny")
  #     saveWorkbook(wb)
  #     file.rename(fname,file)
  #   },
  #   contentType="application/xlsx"
  # )
  
  
  #pobieranie wyniku dla sond
  #proba
  #  content=function(file){
  #  wb <- createWorkbook()
  #  addWorksheet(wb, sheetName = "Geny", gridLines = TRUE)
  #  intestazione <- paste0("jakis tekst", "TOSCANA", ", ", "TOSCANA", ".")
  #  writeData(wb, 1, x = intestazione)
  #  writeDataTable(wb, sheet = 1, startRow = 3, x = wynik_sondy(), colNames = TRUE)
  #  saveWorkbook(wb, file)
  # }

  
  output$downloadData4 <- downloadHandler(
    filename = function() { "dane_geny_sondy.xlsx" },
    content = function(file) {
      
      tempFile <- tempfile(fileext = ".xlsx")
      write.xlsx(wynik_sondy(), tempFile)
      file.rename(tempFile, file)
      file.copy(filename,file)
      
    }
  )
  
  
  output$pobierz_summary_table<-downloadHandler(
    filename = function() {
      #paste(input$tabela(), ".txt", sep = "")
      "summary_table.txt"
    },
    content = function(file) {
      write.table(tabela(), file, row.names = FALSE)
    }
  )
  
  
})