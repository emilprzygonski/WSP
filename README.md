Celem projektu, którego przedstawieniem jest repozytorium było stworzenie aplikacji/programu w celu identyfikacji sygnatur genów różnicujących dwóch grup danych mikromacierzowych oraz analizy ścieżek sygnałowych.

Uruchomienie aplikacji jest możliwe tylko dzięki oprogramowaniu R w wersji 3.3.3 lub 3.3.0 oraz następujących pakietów:

-biobase

-affy

-tcltk

-gahgu95av2.db

-gahgu95av2cdf

-org.Hs.eg.db

-openx1sx

-limma

-gplots

-msigdbr

-GSEABase


Na całość aplikacji składa się szereg funkcji wymaganych do poprawnego działania:

- Create_ExprSet() --> tworzy obiekt ExpSet na podstawie plików .CEL, użytkownik wskazuje folder, w którym znajdują się pliki oraz plik z adnotacjami zawierający nazwy wcześniej wspomnianych plików

- wczytaj_sondy() --> funkcja ma możliwość wczytania listy nazw sond z pliku .txt, gdy użytkownik np.: chce tylko otrzymać ich symbole i entrez_id z funkcji updated_ExprSet

- create_ExprSet2() --> tworzy obiekt ExprSet, gdy użytkownik posiada gotową macierz ekspresji i plik z adnotacjami 

- updated_ExprSet(ExprSet) --> dodaje symbole genów i ich nazwy oraz entrez_id do listy sond (z obiektu ExprSet lub z pliku z listą sond); tabela z danymi zostaje dodana do obiektu ExprSet (featureData) lub jest ona tworzona w przypadku pobrania listy z pliku

- zapis_xlsx(up_ExprSet) --> eksport informacji o sondach do pliku Excel

- summary_table(ExprSet, class1, class2, correction_method, FoldChangeCutoff, sort_critenon, number, sep)  --> analiza cech różnicujących 

* ExprSet - obiekt ExpressionSet

* class1 - 1 klasa obserwacji 

* class2 - 2 klasa obserwacji 

***** użytkownik musi znać oznaczenia klas w ExprSet@phenoData@data[["CLASS"]]

* correction_method - metoda korekcji p-wartości (do wyboru holm, hochberg, hommel, bonferroni, BH, BY, fdr, none)

* FoldChangeCutoff - poziom odcięcia FoldChange (użytkownik podaje liczę większą bądź równą 1)

* sort_critenon - poziom ufności p-wartości (np.: 0.05)

* number - liczba wyników do zwrócenia (0 zwraca wszystkie)

* sep - separator przy zapisie do pliku .txt

- heatmapa(ExprSet, TAB_geny_roznicujace) --> wyrysowanie heatmapy ekspresji dla sond różnicujących i zapis do pliku .png

