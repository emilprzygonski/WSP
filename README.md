# WSP
Repozytorium projektu przedmiotu WSP 

Wytyczne projektu

Stworzyć program dostępny z poziomu strony WWW realizujący następujące zadania:

1) Identyfikacja sygnatur genowych różnicujących dla dwóch grup danych mikromacierzowych

a. Możliwość wczytania standardowego obiektu klasy ExpressionSet (ESet)

2) Analiza ścieżek sygnałowych powiązanych z dowolną listą sygnatur genowych

a. Możliwość wczytania list ścieżek ze strony (przykład: c2.all.v6.2.entrez.gmt): http://software.broadinstitute.org/gsea/msigdb/collections.jsp

b. Możliwość wykorzystania własnej listy ścieżek

Przydatne pakiety:

library(stringr)

library(limma)

library(GSA)

library(piano)

library(openxlsx)

library(shiny)


Projekt powinien zostać umieszczony w dowolnym repozytorium git.

Wymagane elementy:

a) Konwersja oznaczeń sond Affymetrix popularnych mikromacierzy na symbole genów, EntrezId, ewentualnie opisowe, nazwy genów (możliwość eksportu danych do postaci pliku Excel). Funkcja powinna działać dla argumentu w postaci obiektu ExpressionSet, plik tekstowy w postaci listy nazw sond)

b) Selekcja zestawów sond różnicujących wykorzystując proste metody statystyczne (np. test t i krotność zmian)

c) Dla dowolnej listy genów (patrz punkt a,b) program powinien pozwolić na analizę ścieżek sygnałowych

a. Funkcja importu ścieżek sygnałowych z kolekcji msigdb. Wymagana możliwość wczytania ręcznie zdefiniowanej list, zalecana możliwość automatycznego inportu)

b. Funkcja analizy ścieżek sygnałowych za pomocą metody GSA, opcjonalnie GSEA

d) Stworzenie interfejsu graficznego dla wszystkich elementów opisanych w podpunktach a)-c). Wymagania: dostęp przez stronę www, reaktywność działania, aplikacja typu klient-serwer. Zalecane wykorzystanie biblioteki shiny.


