# Projet de MAP431 : _Écoulements visqueux_

Ce dépôt contient le travail relatif au projet _Écoulements visqueux_ du cours MAP431 _Analyse variationnelle des équations aux dérivées partielles_,  proposé par Martin Averseng.

On y trouvera le sujet au format PDF, les fichiers source LaTeX du rapport, le PDF du rapport, les scripts `FreeFem++` et un script pour générer des animations au format GIF.

## Script pour animation

Le code `FreeFem++` exporte les graphes au format `.eps`. Pour avoir une visualisation de la déformation du domaine, on convertit ces fichiers en un `.gif` en utilisant [ImageMagick](www.imagemagick.org) dans la ligne de commande (un script `togif.sh` est donné).
