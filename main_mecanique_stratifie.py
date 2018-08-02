# -*- coding: utf-8 -*-
"""
Cree par    : Jean COLLOMB
Societe     : CT1 - groupe Compose
Laboratoire : SYMME - Polytech Annecy-Chambéry - Université Savoie Mont Blanc
Date        : janvier 2017

Ce script a pour objectif la determination des proprietes mecaniques
homogeneisees d'un stratifie composite, a l'aide des proprietes 
mecaniques des plis constitutifs, ainsi que des angles de drapage.

FONCTIONNEMENT :
----------------
    Etape 1 :   Renseigner l'empilement des plis et les proprietes mecaniques
                dans le fichier Excel 'empilement.xlsx'
    Etape 2 :   Lecture du fichier Excel par Python et recuperation des 
                informations
    Etape 3 :   Calculs MSC d'homogeneisation mecanique
    Etape 4 :   Stockage des resultats dans un fichier .txt
"""

#------------------------ Importation des Modules
import time
import pandas as pd
import numpy as np
from homogeneisation_mecanique_stratifie import Homogeneisation_Mecanique_Stratifie



#------------------------ Extraction des informations

temps_depart_extraction = time.time()
informations = pd.read_excel('empilement.xlsx')


"""
Stockage de chaque propriete dans des listes.
"""
liste_angles        = informations['angle'].tolist()
liste_epaisseurs    = informations['epaisseur'].tolist()
liste_El            = informations['El'].tolist()
liste_Et            = informations['Et'].tolist()
liste_Glt           = informations['Glt'].tolist()
liste_Nult          = informations['Nult'].tolist()
liste_pli           = informations['Pli'].tolist()
liste_proportion    = informations['proportion (%)'].fillna(0).tolist()
temps_fin_extraction= time.time()


#------------------------ Realisation des calculs

"""
Realisation du calcul a l'aide de la Class Homogeneisation Mecanique
presente dans le fichier homogeneisation_mecanique.py
Extraction des proprietes homogeneisees du stratifie.
"""
temps_depart_calculs = time.time()
Calcul = Homogeneisation_Mecanique_Stratifie(liste_angles, liste_epaisseurs, liste_El, liste_Et, liste_Glt, liste_Nult, liste_proportion)
Ex_stratifie = Calcul.proprietes_mecaniques_stratifie()[0]
Ey_stratifie = Calcul.proprietes_mecaniques_stratifie()[1]
Gxy_stratifie = Calcul.proprietes_mecaniques_stratifie()[2]
Nuxy_stratifie = Calcul.proprietes_mecaniques_stratifie()[3]
temps_fin_calculs = time.time()


#------------------------ Exportation des resultats

#--- Fichier texte
temps_depart_ecriture = time.time()
resultats_txt = open("resultats.txt", "w")

resultats_txt.write("----------------------")
resultats_txt.write("\n")
resultats_txt.write("-------------RESULTATS")
resultats_txt.write("\n")
resultats_txt.write("----------------------")
resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("Ci-après, les resultats obtenus pour l'homogénéisation des propriétés mécaniques du stratifié")
resultats_txt.write("\n")
resultats_txt.write("Les calculs MSC permettent d'estimer les proprietes mecaniques finales du stratifie composite \n")
resultats_txt.write("en se basant sur les proprietes mecaniques de chaque pli. \n")
resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("Le stratifie etudie est compose de : " + str(len(liste_pli)) + " plis")
resultats_txt.write("\n")
resultats_txt.write("L'epaisseur du stratifie est de : " + str(round(Calcul.epaisseur_stratifie,3)) + " mm")

resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("-----Donnees d'entrees\n")
resultats_txt.write("----------------------\n")
pd.set_option('max_rows', len(liste_pli))
informations.set_index('Pli', inplace = True)
resultats_txt.write(str(informations))

resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("-----Donnees intermediaires")
resultats_txt.write("\n")
resultats_txt.write("---------------------------")
resultats_txt.write("\n")
resultats_txt.write("Matrice A         : ")
resultats_txt.write("\n")
resultats_txt.write(str(np.around(sum(Calcul.liste_A()), decimals = 1)))
resultats_txt.write("\n")
resultats_txt.write("Matrice B         : ")
resultats_txt.write("\n")
resultats_txt.write(str(np.around(sum(Calcul.liste_B()), decimals = 1)))
resultats_txt.write("\n")
resultats_txt.write("Matrice D         : ")
resultats_txt.write("\n")
resultats_txt.write(str(np.around(sum(Calcul.liste_D()), decimals = 1)))

resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("-----Donnees de sorties")
resultats_txt.write("\n")
resultats_txt.write("-----------------------")
resultats_txt.write("\n")
resultats_txt.write("Module Ex (MPa)    : " + str(round(Ex_stratifie, 1)))
resultats_txt.write("\n")
resultats_txt.write("Module Ey (MPa)    : " + str(round(Ey_stratifie, 1)))
resultats_txt.write("\n")
resultats_txt.write("Module Ez (MPa)    : " + str(round(np.mean(liste_Et), 1)))
resultats_txt.write("\n")
resultats_txt.write("Module Gxy (MPa)   : " + str(round(Gxy_stratifie, 1)))
resultats_txt.write("\n")
resultats_txt.write("Coefficient Nuxy   : " + str(round(Nuxy_stratifie, 3)))

resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("\n")
resultats_txt.write("-----Information")
resultats_txt.write("\n")
resultats_txt.write("-----------------------")
resultats_txt.write("\n")
temps_fin_ecriture = time.time()
resultats_txt.write("\n")
resultats_txt.write("temps d'extraction des donnees : " + str(round(temps_fin_extraction - temps_depart_extraction, 3)) + " secondes")
resultats_txt.write("\n")
resultats_txt.write("temps de calculs               : " + str(round(temps_fin_calculs - temps_depart_calculs, 3)) + " secondes")
resultats_txt.write("\n")
resultats_txt.write("temps d'ecriture               : " + str(round(temps_fin_ecriture - temps_depart_ecriture, 3)) + " secondes")
resultats_txt.write("\n")
resultats_txt.write("temps total                    : " + str(round(temps_fin_ecriture - temps_depart_extraction, 3)) + " secondes")

resultats_txt.close()







