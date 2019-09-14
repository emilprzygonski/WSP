library(plotly)
library(shinyFiles)
library(data.table)
library(shinythemes)
inline = function (x) {
  tags$div(style="display:inline-block;", x)
}

shinyUI(fluidPage(theme=shinytheme("sandstone"),
  titlePanel("ANALIZA SCIEZEK SYGNALOWYCH"),
  tabsetPanel(
    tabPanel("Informacje o aplikacji",
             br(),
             wellPanel(
             h4("Aplikacja umozliwia selekcje 
                zestawow sond roznicujacych dla dwoch grup danych mikromacierzowych
                oraz analize sciezek sygnalowych."),
             br(),
             h3("ANALIZA ROZNICUJACA", style = "color:purple"),
             h5("W celu wykonania", tags$b("analizy roznicujacej"), "nalezy wybrac plik 
                z macierza ekspresji oraz plik z adnotacjami w panelu
                'Identyfikacja sygnatur genowych dla obiektu ExpressionSet'."),
             h5("W panelu w poszczegolnych zakladkach uzytkownik ma mozliwosc 
                podgladu wczytanych danych oraz tabeli cech, ktora zawiera 
                wynik konwersji oznaczen sond Affymetrix popularnych mikromacierzy 
                na symbole genow, nazwy genow, EntrezId, opisowe.
                W zakladce 'Pobierz' jest mozliwosc pobrania wynikowej tabeli cech."),
             h5("Kolejnym krokiem w analizie roznicujacej jest przejscie do panelu 
                'Selekcja zestawow sond roznicujacych'. 
                W pierwszej zakladce, po wybraniu ww. panelu, zostaje wyswietlona 
                tabela zawierajaca wynik analizy roznicujacej. 
                Zakladka 'Pobierz wynik' umozliwia pobranie tabeli. 
                Z kolei w zakladce 'Heatmapa' jest mozliwosc pobrania wygenerowanej 
                na podstawie analizy roznicujacej heatmapy do folderu, 
                w ktorym znajduja sie pliki aplikacji, tj. serwer.R oraz ui.R."),
             br(),
             h3("ANALIZA SCIEZEK SYGNALOWYCH", style = "color:purple"),
             h5("Aby wykonac", tags$b("analize sciezek sygnalowych"), "nalezy przejsc do panelu 
                'Analiza sciezek sygnalowych', a nastepnie, po lewej stronie, 
                wybrac kategorie. Kolejny krok polega na imporcie sciezek sygnalowych 
                z kolekcji msigdb. Aby automatycznie pobrac sciezki sygnalowe, 
                nalezy kliknac przycisk 'Importuj sciezki sygnalowe'. 
                W zakladce 'Analiza sciezek sygnalowych za pomoca metody GSA' 
                nalezy wczytac liste sciezek klikajac przycisk 'Wybierz plik z sciezkami', 
                a nastepnie wybrac przycisk 'Wykonaj analize GSA'."))
    ),
    tabPanel("Identyfikacja sygnatur genowych dla obiektu ExpressionSet",
             br(),
        
             sidebarPanel(
              
               br(), 
              
               fileInput("input2", "Wybierz plik z adnotacjami:", multiple=FALSE,accept="txt"),
               
               br(),
               fileInput('input1', 'Wybierz plik z macierza ekspresji:',
                         accept = c( 'text/csv','text/comma-separated-values','.csv')
               )
               ),

             mainPanel(
               tabsetPanel(
                 tabPanel("Wczytana macierz ekspresji",
                          tableOutput("plik1")
                         ),
                 tabPanel("Wczytane adnotacje",
                          tableOutput("plik2")
                          ),
                 tabPanel("Wynikowa tabela cech",
                          br(),br(),
                          tableOutput("tabela_cechy_macierz")
                          ),
                 tabPanel("Pobierz",
                          br(),br(),br(),
                          downloadButton("downloadData1", "Pobierz plik z macierza ekspresji .csv"),
                          br(),br(),br(),
                          downloadButton("downloadData2", "Pobierz plik z adnotacjami .txt"),
                          br(),br(),br(),
                          downloadButton("downloadData3", "Pobierz tabele cech .xlsx")
                         )
               ))),
    tabPanel("Selekcja zastawow sond roznicujacych",
             br(),
             sidebarPanel(
               h4("Wprowadz wymagane dane do analizy roznicujacej"),
               h6("Uwaga: W polach 'Pierwsza klasa obserwacji do analizy roznicujacej' oraz 'Druga klasa obserwacji do analizy roznicujacej' nalezy wpisac klasy obserwacji uzywajac cudzyslow"),
               br(),
               textInput(inputId="class1", label="Pierwsza klasa obserwacji do analizy roznicujacej:", value = 'ADENO', width = NULL,
                         placeholder = NULL),
               textInput(inputId="class2", label="Druga klasa obserwacji do analizy roznicujacej:", value = 'SQUAMOUS', width = NULL,
                         placeholder = NULL),
               selectInput(inputId="correction_method", label="metoda korekcji:", 
                           choices=list("holm","hochberg","hommel",
                                        "bonferroni","BH","BY",
                                        "fdr","none"), 
                           selected = 'holm', multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
               numericInput(inputId="FoldChangeCutoff", label="ekspresja w stosunku do kontroli:", value=1.2, min = 1, max = NA, step = 0.1,
                            width = NULL),
               numericInput(inputId="sort_criterion", label="poziom ufnosci dla p-wartosci:", value=0.01, min = 0, max = NA, step = 0.01,
                            width = NULL),
               numericInput(inputId="number", label="liczba wynikow do wyswietlenia:", value=0, min = 0, max = NA, step = NA,
                            width = NULL),
               radioButtons(inputId="sep", label="separator:", 
                            choices = list(", ", ";"),
                            selected = ', ',inline = FALSE, width = NULL, choiceNames = NULL,
                            choiceValues = NULL)
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Selekcja zestawu sond roznicujacych",
                          tableOutput("summary_table")
                          
                          ),
                 tabPanel("Pobierz wynik",
                          br(),
                          downloadButton("pobierz_summary_table", "Pobierz zestaw sond roznicujacych")
                          ),
                 tabPanel("heatmapa",
                          br(),
                          h5("Po nacisnieciu przycisku 'Utworz i pobierz heatmape', heatmapa zostanie zapisana w folderze, w ktorym znajduja sie pliki aplikacji Shiny."),
                          actionButton("create_heatmap", "Utworz heatmape")
                          #imageOutput("heatmapa_image")
                          )
               )
             )
           ),
    tabPanel("Analiza sciezek sygnalowych",
             br(),br(),
             sidebarPanel(
             h4("IMPORT SCIEZEK SYGNALOWYCH", style = "color:purple"),
             h5("W celu zaimportowania sciezek sygnalowych z kolekcji msigdb, 
                wybierz kategorie dostepne w bazie danych"),
             br(),
             selectInput(inputId="kategoria", label="Wybierz kategorie", 
                         choices=list("kategoria H"= 1,"kategoria C1" = 2,"kategoria C2" = 3,
                                      "kategoria C2, subkategoria CGP" = 4,
                                      "kategoria C2, subkategoria CP" = 5,
                                      "kategoria C2, subkategoria CP:BIOCARTA" = 6,
                                      "kategoria C2, subkategoria CP:KEGG" = 7,
                                      "kategoria C2, subkategoria CP:REACTOME" = 8,
                                      "kategoria C3" = 9,
                                      "kategoria C3, subkategoria MIR" =10,
                                      "kategoria C3, subkategoria TFT" = 11,
                                      "kategoria C4" = 12,
                                      "kategoria C4, subkategoria CGN" = 13,
                                      "kategoria C4, subkategoria CM" = 14,
                                      "kategoria C5" = 15,
                                      "kategoria C5, subkategoria BP" = 16,
                                      "kategoria C5, subkategoria CC" = 17,
                                      "kategoria C5, subkategoria MF" = 18,
                                      "kategoria C6" = 19,
                                      "kategoria C7" = 20
                                      ), 
                         selected = 'kategoria H', multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
             br(),
             h4("ANALIZA GSA",style = "color:purple"),
             h5("Aby przeprowadzic", tags$b("analize GSA"), "wczytaj ponizej
                plik zawierajacy sciezki sygnalowe (zapisany z importu EntrezIDs.gmt lub wlasny plik 
                - przykldowy pobrany recznie z bazy - msigdb_entrez_id.gmt)
                oraz wybierz, dla jakich genow ma zostac przeprowadzona analiza."),
             h5("Nastepnie wybierz Panel 'Analiza sciezek sygnalowych za pomoca metody GSA' 
                i kliknij przycisk 'Przeprowadz analize GSA'"),
             br(),
             br(),
             fileInput('input_gsa', 'Wybierz plik ze sciezkami sygnalowymi:', multiple=FALSE,accept=".gmt"),
             selectInput(inputId="kategoria2", label="Wybierz geny z obiektu ExpressionSet:", 
                         choices=list("wszystkie geny"=1,"tylko roznicujace geny"=2),
                         selected = "wszystkie geny", multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
             br(),
             h6("Uwaga: Aby wykonac analize GSA, niezbedne jest wczesniejsze uzyskanie 
             wynikowej tabeli cech w Panelu 'IDENTYFIKACJA SYGNATUR GENOWYCH DLA OBIEKTU 
             EXPRESSIONSET',w zakladce 'WYNIKOWA TABELA CECH', a takze uzyskanie wyniku 
             selekcji sond roznicujacych poprzez wybranie Panelu 
             'SELEKCJA ZESTAWOW SOND ROZNICUJACYCH'")
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Import sciezek sygnalowych",
                          br(),br(),
                          h5("Aplikacja umozliwia zaimportowanie sciezek sygnalowych 
                             z kolekcji msigdb oraz zapisanie ich do folderu, w ktorym 
                             znajduja sie pliki aplikacji (serwer.R, ui.R)."),
                             
                          h5("W celu zaimportowania listy sciezek sygnalowych, 
                             kliknij przycisk 'Importuj sciezki sygnalowe z kolekcji msigdb'"),
                          
                          h5("UWAGA! Proce moze chwile potrwac, prosimy o cierpliwosc. Po kilku minutach 
                             w folderze z projektem powinien pojawic sie plik EntrezIDs.gmt, ktory mozna wczytac do analizy GSA"),
                          br(),br(),
                          actionButton("import_z", "Importuj sciezki sygnalowe z kolekcji msigdb"),
                          br(),br(),
                          
                          #DO USUNIECIA?
                          # downloadButton("zapisz_sciezki", "Pobierz sciezki sygnalowe"),
                          # br(),br(),
                          # 
                          
                          tableOutput("gsa")
                          
                          
                          ),
                 tabPanel("Analiza sciezek sygnalowych za pomoca metody GSA",
                          br(),br(),
                          h5("Analiza GSA moze potrwac do kilkanastu minut."),
                          h5("Zaleca sie wczytanie pliku ze sciezkami
                             sygnalowymi (opcja dostepna po lewej stronie) w pierwszej
                             kolejnosci."),
                          h5("Nastepnie, w przypadku, gdy uzytkownik
                             posiada obiekt ExprSet, nalezy wybrac macierz ekspresji i adnotacje
                             z Panelu 'IDENTYFIKACJA SYGNATUR GENOWYCH
                             DLA OBIEKTU EXPRESSIONSET'."),
                          h5("W przypadku, gdy uzytkownik posiada liste sond, nalezy przejsc do Panelu
                          'IDENTYFIKACJA SYGNATUR GENOWYCH DLA LISTY SOND' i wybrac plik."),
                          h5("Ostni etap obejmuje przejscie do Panelu
                             'ANALIZA SCIEZEK SYGNALOWYCH', wybranie zakladki
                             'ANALIZA SCIEZEK SYGNALOWYCH ZA POMOCA METODY GSA'
                             oraz klikniecie przycisku 'PRZEPROWADZ ANALIZE GSA.'"),
                          h5("W przypadku zbyt malej liczby sond roznicujacych 
                             wprowadzanych do analizy GSA (wowczas mniej niz 5% 
                             entrez id pokrywa sie ze sciezkami) moze pojawic sie blad."),
                          br(),br(),
                          actionButton("wykonaj_gsa","Przeprowadz analize GSA")
                           ),
                 tabPanel("WYNIKI ANALIZY GSA",
                          tableOutput("wynik_analizy")
                          )
               )
             )
             ),
    tabPanel("Identyfikacja sygnatur genowych dla listy sond",
             br(),
             sidebarPanel(
               br(),
               fileInput('input3','Wybierz plik z lista sond:',accept='txt')
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Wczytana lista sond",
                          tableOutput("tabela_sondy")
                 ),
                 tabPanel("Wynikowa tabela cech",
                          tableOutput("tabela_cechy_sondy")
                 ),
                 tabPanel("Pobierz wyniki",
                          downloadButton("downloadData4", "Pobierz tabele cech .xlsx")
                          )
               )
             )
    )
    
)))


