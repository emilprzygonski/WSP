# WSP
Repozytorium projektu przedmiotu WSP 

Wytyczne projektu

Stworzyć program dostępny z poziomu strony WWW realizujący następujące zadania:

Identyfikacja sygnatur genowych różnicujących dla dwóch grup danych mikromacierzowych

Możliwość wczytania standardowego obiektu klasy ExpressionSet (ESet)

Analiza ścieżek sygnałowych powiązanych z dowolną listą sygnatur genowych

Możliwość wczytania list ścieżek ze strony (przykład: c2.all.v6.2.entrez.gmt): http://software.broadinstitute.org/gsea/msigdb/collections.jsp

Możliwość wykorzystania własnej listy ścieżek

Przydatne pakiety:

library(stringr)

library(limma)

library(GSA)

library(piano)

library(openxlsx)

library(shiny)


Projekt powinien zostać umieszczony w dowolnym repozytorium git.

Wymagane elementy:

Konwersja oznaczeń sond Affymetrix popularnych mikromacierzy na symbole genów, EntrezId, ewentualnie opisowe, nazwy genów (możliwość eksportu danych do postaci pliku Excel). Funkcja powinna działać dla argumentu w postaci obiektu ExpressionSet, plik tekstowy w postaci listy nazw sond)

Selekcja zestawów sond różnicujących wykorzystując proste metody statystyczne (np. test t i krotność zmian)

Dla dowolnej listy genów (patrz punkt a,b) program powinien pozwolić na analizę ścieżek sygnałowych

Funkcja importu ścieżek sygnałowych z kolekcji msigdb. Wymagana możliwość wczytania ręcznie zdefiniowanej list, zalecana możliwość automatycznego inportu)

Funkcja analizy ścieżek sygnałowych za pomocą metody GSA, opcjonalnie GSEA

Stworzenie interfejsu graficznego dla wszystkich elementów opisanych w podpunktach a)-c). Wymagania: dostęp przez stronę www, reaktywność działania, aplikacja typu klient-serwer. Zalecane wykorzystanie biblioteki shiny.


