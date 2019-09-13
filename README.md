Celem projektu, którego przedstawieniem jest repozytorium było stworzenie aplikacji/programu w celu identyfikacji sygnatur genów różnicujących dwóch grup danych mikromacierzowych oraz analizy ścieżek sygnałowych.

Podział pracy w grupie:
1)	Identyfikacja sygnatur genowych różnicujących dla dwóch grup danych mikromacierzowych 
> NOCUŃ WIKTORIA, KUKLA ARTUR

2)	Analiza ścieżek sygnałowych powiązanych z  listą sygnatur genowych
> WENZ KALINA, PRZYGOŃSKI EMIL
 
3)  Stworzenie interfejsu graficznego 
> BONIEK ANNA, WIECZOREK MARTA 

Uruchomienie aplikacji jest zalecane w oprogramowaniu R w wersji 3.3.3 lub 3.3.0, gdyż w tych wersjach kompatybilne są wymagane pakiety: 

*Biobase
*affy
*tcltk
*gahgu95av2.db
*gahgu95av2cdf
*org.Hs.eg.db
*openxlsx
*limma
*gplots
*msigdbr
*GSEABase
*GSA
*Rtools (Rtools 3.4 dla wersji R 3.3.x)
*plotly
*shiny
*shinyFiles
*shinythemes
*shinyTime
*log4r
*oligo


Opis interfejsu aplikacji:

UWAGA! Prosimy o cierpliwość podczas obsługiwania aplikacji. Wczytywana macierz ekspresji zazwyczaj ma duży rozmiar i jej wczytanie może wymagać chwili czasu. Zalecamy również w pierwszej kolejności wczytać plik z adnotacjami, dopiero później z macierzą ekspresji.

Aplikacja umozliwia selekcje zestawow sond roznicujacych dla dwoch grup danych mikromacierzowych oraz analize sciezek sygnalowych.

ANALIZA ROZNICUJACA
W celu wykonania analizy roznicujacej nalezy wybrac plik z macierza ekspresji oraz plik z adnotacjami w panelu 'Identyfikacja sygnatur genowych dla obiektu ExpressionSet'.

W panelu w poszczegolnych zakladkach uzytkownik ma mozliwosc podgladu wczytanych danych oraz tabeli cech, ktora zawiera wynik konwersji oznaczen sond Affymetrix popularnych mikromacierzy na symbole genow, nazwy genow, EntrezId, opisowe. W zakladce 'Pobierz' jest mozliwosc pobrania wynikowej tabeli cech.

Kolejnym krokiem w analizie roznicujacej jest przejscie do panelu 'Selekcja zestawow sond roznicujacych'. W pierwszej zakladce, po wybraniu ww. panelu, zostaje wyswietlona tabela zawierajaca wynik analizy roznicujacej. Zakladka 'Pobierz wynik' umozliwia pobranie tabeli. Z kolei w zakladce 'Heatmapa' jest mozliwosc pobrania wygenerowanej na podstawie analizy roznicujacej heatmapy do folderu, w ktorym znajduja sie pliki aplikacji, tj. serwer.R oraz ui.R.

ANALIZA SCIEZEK SYGNALOWYCH
Aby wykonac analize sciezek sygnalowych nalezy przejsc do panelu 'Analiza sciezek sygnalowych', a nastepnie, po lewej stronie, wybrac kategorie. Kolejny krok polega na imporcie sciezek sygnalowych z kolekcji msigdb. Aby automatycznie pobrac sciezki sygnalowe, nalezy kliknac przycisk 'Importuj sciezki sygnalowe'. W zakladce 'Analiza sciezek sygnalowych za pomoca metody GSA' nalezy wczytac liste sciezek z pliku "msigdb_entrez_id.gmt" klikajac przycisk 'Wybierz plik z sciezkami', a nastepnie wybrac przycisk 'Wykonaj analize GSA'.


Na całość aplikacji składa się szereg funkcji wymaganych do poprawnego działania:

- Create_ExprSet() --> tworzy obiekt ExpSet na podstawie plików .CEL, użytkownik wskazuje folder, w którym znajdują się pliki oraz plik z adnotacjami zawierający nazwy wcześniej wspomnianych plików. 
UWAGA: FUNKCJA TA ZOSTAŁA WYKONANA NA PLIKACH .CEL POZA APLIKACJĄ, PONIEWAŻ SHINY NIE BYŁO W STANIE POMIEŚCIĆ TAK DUŻYCH PLIKÓW W PAMIĘCI. STWORZONO ZATEM WZORCOWY ExpressionSet Z UDOSTĘPNIONYCH PRZEZ PROWADZĄCEGO PLIKÓW .CEL I ZAPISANO JEGO MACIERZ EKSPRESJI W ZAŁĄCZONYM DO PROJEKTU PLIKU _____________.

- create_ExprSet2() --> tworzy obiekt ExprSet, gdy użytkownik posiada gotową macierz ekspresji i plik z adnotacjami. 
Wzorcowa macierz ekspresji ze wspomnianych wcześniej plików .CEL została załączona do projektu w pliku ________ oraz plik z adnotacjami __________.

- wczytaj_sondy() --> funkcja ma możliwość wczytania listy nazw sond z pliku .txt, gdy użytkownik np.: chce tylko otrzymać ich symbole i entrez_id z funkcji updated_ExprSet(). 

- updated_ExprSet(ExprSet) --> dodaje symbole genów i ich nazwy oraz entrez_id do listy sond (z obiektu ExprSet lub z pliku z listą sond); tabela z danymi zostaje dodana do obiektu ExprSet (featureData) lub jest ona tworzona w przypadku pobrania listy sond z pliku.
Z niewiadomych przyczyn po wczytaniu do tej funkcji wyniku funkcji wczytaj_sondy(), czyli listy nazw sond, aplikacja wyświetla błąd. Dzieje się tak pomimo właściwego działania kodu w środowisku R poza pakietem shiny.

- zapis_xlsx(up_ExprSet) --> eksport informacji o sondach do pliku Excel.

- summary_table(ExprSet, class1, class2, correction_method, FoldChangeCutoff, sort_criterion, number, sep)  --> analiza cech różnicujących (test t lub Foldchange). 

** ExprSet - obiekt ExpressionSet

** class1, class2 - klasy obserwacji w danych, które użytkownik wpisuje "z ręki" (powinien znać nazwy klas które chce porównywać)

***** użytkownik musi znać oznaczenia klas w ExprSet@phenoData@data[["CLASS"]]

** correction_method - metoda korekcji p-wartości (do wyboru holm, hochberg, hommel, bonferroni, BH, BY, fdr, none)

** FoldChangeCutoff - poziom odcięcia FoldChange (użytkownik podaje liczę większą bądź równą 1)

** sort_criterion - poziom ufności dla p-wartości (np.: 0.05)

** number - liczba wyników do zwrócenia (0 zwraca wszystkie)

** sep - separator przy zapisie wyników do pliku .txt

- heatmapa(ExprSet, TAB_geny_roznicujace) --> wyrysowanie heatmapy ekspresji dla sond różnicujących i zapis do pliku .png

- import_zmsigdbr() --> funkcja umożliwia import danych z bazy MSigDB, użytkownik wybiera kategorię, która go interesuje. Następnie tworzone są pliki .gmt potrzebne do analizy sygnałowej GSA.
NIESTETY TWORZONE PLIKI .GMT NIE DO KOŃCA PASUJĄ POD WZGLĘDEM STRUKTURALNYM DO WEJŚCIA WYKORZYSTYWANEJ PÓŹNIEJ FUNKCJI GSA(). Z TEGO WZGLĘDU ZAREJESTROWANO SIĘ DO BAZY MSigDB I POBRANO PLIK _______.GMT ZAWIERAJĄCY WSZYSTKIE KATEGORIE ŚCIEŻEK I WYKORZYSTANO GO W ANALIZIE ŚCIEŻEK SYGNAŁOWYCH W FUNKCJI przeprowadzanie_GSA().

*****
Lista kategorii i oznaczenia (jak po myślniku to subkategoria)

  #H, C1, C2, C2-CGP 
  
  #C2-CP, C2-CP:BIOCARTA, C2-CP:KEGG, C2-CP:PID, C2-CP:REACTOME
  
  #C3, C3-MIR, C3-TFT, C4, C4-CGN, C4-CM, C5, C5-BP, C5-CC
  
  #C5-MF, C6, C7

- przeprowadzenie_GSA() - przeprowadzenie analizy ścieżek sygnałowych metodą GSA.
ZE WZGLĘDU NA NIE DO KOŃCA PASUJĄCĄ STRUKTURĘ TWORZONYCH W FUNKCJI import_zmsigdbr() PLIKÓW .GMT, WYKORZYSTYWANY JEST TU PLIK _______.GMT, POBRANY Z BAZY MSigDB, ZAWIERAJĄCY WSZYSTKIE KATEGORIE ŚCIEŻEK.
