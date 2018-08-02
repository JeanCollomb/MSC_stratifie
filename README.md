# MSC_stratifie

Lors de la conception et du dimensionnement de structures, il est préférable dans un premier temps de réaliser 
une homogénéisation du composite afin de simplifier l’étude et de réduire les temps de calculs.
L’homogénéisation consiste au passage d’une structure hétérogène (multi-matériaux) 
à une structure homogène (création d’un matériau équivalent).

Ce script python ce propose de faire se travail.


## Pour commencer

Récupérer les fichiers :
- main_mecanique_stratifie.py
- homogeneisation_mecanique_stratifie.py
- empilement.xlsx


```
Pré-requis :
- Posséder un tableur type Excel
- Posséder les packages : numpy, pandas, time
```


## Fonctionnement

Deux utilisations principales sont prévues par le script.

### Empilement connu

Dans ce premier cas, l'empilement du stratifié est connu.

* Renseigner dans "empilement.xlsx" les différents plis (propriétés matériaux, épaisseurs plis, orientations).
* Laisser la colonne proportion vide ou mettre à 0.
* Excécuter "main_mecanique_stratifie.py"

```
Excécution par l'invite de commande dans le dossier de travail :
python main_mecanique_stratifie.py
```

### Proportions connues

Dans ce second cas, l'empilement du stratifié est inconnu, mais les proportions de chacun de angle l'est.

* Renseigner dans "empilement.xlsx" les différentes orientations de plis (propriétés matériau, épaisseurs pli, orientation).
* Renseigner la proportion pour chacune des orientations
* Excécuter "main_mecanique_stratifie.py"

```
Excécution par l'invite de commande dans le dossier de travail :
python main_mecanique_stratifie.py
```

### Résultats

Après excécution du fichier "main_mecanique_stratifie.py", un fichier texte est généré dans le dossier de travail.
Ce fichier contient les valeurs homogénéisées du stratifié.

