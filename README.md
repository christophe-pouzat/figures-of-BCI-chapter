# figures-of-BCI-chapter
R code for regenerating the figures of "Extracellular data analysis" in book Brain-Computer Interfaces

## English version

This repository contains the `R` code necessary for re-generating all (except one) of the figures of chapter _Analysis of extracellular data_ in book _Brain-Computer Interfaces 1: Methods and Perspectives_ edited by Maureen Clerc, Laurent Bougrain and Fabien Lotte, [iSTE/Wiley](http://eu.wiley.com/WileyCDA/WileyTitle/productCd-1848218265.html).

The relevant files are:

- `do_sorting.R`, that takes care of getting the (real) data and doing the spike sorting.
- `the_figures.R`, that takes the oputput of the previous file and generates the figures.
- `BCI_Pouzat_Les_Figures.org`, the source of the previous two files.

The (very detailed) code generating the second figure can be found on my [LASCON course page](http://christophe-pouzat.github.io/LASCON2016/OriginOfTheHighFrequencyExtraCellularSignal.html).

## Version française

Ce dépôt contient le code `R` nécessaire à la ré-génération de toutes les figures, sauf une, de mon chapitre « Analyse des données extracellulaires » du livre « Les interfaces cerveau-ordinateur 1. Fondements et méthodes. » édité par Maureen Clerc, Laurent Bougrain et Fabien Lotte chez [iSTE editions](http://iste-editions.fr/products/les-interfaces-cerveau-ordinateur-1).

Les fichiers pertinents sont :

- `do_sorting.R`, qui télécharge les données réelles et effectue le tri des potentiels d'action ;
- `les_figures.R`, qui prend le résultat de l'exécution du code précédent et génère les figures ;
- `BCI_Pouzat_Les_Figures.org`, le fichier source à partir duquel les deux autres sont générés.

Une version très détaillée (mais en anglais) du code générant la deuxième figure se trouve sur [la page de mon cours au LASCON](http://christophe-pouzat.github.io/LASCON2016/OriginOfTheHighFrequencyExtraCellularSignal.html).

