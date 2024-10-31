# Algorytm-Neddlman


## Opis
Program implementuje algorytm Needleman-Wunsch do globalnego dopasowania dwóch sekwencji aminokwasowych. Umożliwia analizę podobieństwa sekwencji oraz oblicza procent identyczności.

## Sposób użycia
Aby uruchomić program, należy przekazać jako argument ścieżkę do pliku FASTA zawierającego dwie sekwencje do dopasowania.

**Składnia:**
```bash
python nw.py <ścieżka_do_pliku_fasta>

##Przykład użycia
>sekwencja1
MGSSHHHHHHSQDPMHKLAIIGYGAAGFAAMIK
>sekwencja2
YKKAREVTGSEVYPPFSSFQEKDGLVQEMRKT
Po uruchomieniu programu poleceniem:

python nw.py para_sekwencji.fasta

Oczekiwane wyjście zapisane w pliku alignment_output.txt:

Dopasowanie:
MGSSHHHHHHS-QDPMHKLAIIGYGAAGFAAMIK-
-YKKAREVTGSEVYPPFSSFQEKDG--LVQEMRKT
Procent identyczno ci: 14.29%
