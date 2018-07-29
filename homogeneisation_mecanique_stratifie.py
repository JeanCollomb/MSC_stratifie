# -*- coding: utf-8 -*-
"""
Cree par    : Jean COLLOMB
Societe     : CT1 - groupe Compose
Laboratoire : SYMME - Polytech Annecy-Chambéry - Université Savoie Mont Blanc
Date        : janvier 2017

Fichier contenant la Class Homogeneisation_Mecanique.
Cette Class contient l'ensemble des fonctions permettant le calcul
des proprietes mecaniques homogeneisees du stratifie.
"""

#------------------------ Importation des Modules
import numpy as np


#------------- Creation de la Class Homogeneisation_Mecanique
class Homogeneisation_Mecanique_Stratifie():
    """
    Class regroupant les fonctions permettant les calculs MSC
    d'homogeneisation mecanique.
    ------
    Donnees d'entrees :
    liste_angles, liste_epaisseurs, liste_El, liste_Et, liste_Glt, liste_Nult
    ------
    liste_angles [list] : Liste comportant les angles de drapage
    liste_epaisseurs [list] : Liste comportant les epaisseurs des plis
    liste_El [list] : Liste comportant les Modules longitudinaux
    liste_Et [list] : Liste comportant les Modules transversaux
    liste_Glt [list] : Liste comportant les Modules de cisaillement
    liste_Nult [list] : Liste comportant les coefficients de poissons
    ------
    Exemple : 
    angle = [0, 90, 90, 0]
    epaisseurs = [0.15, 0.15, 0.15, 0.15]
    El = [150000, 150000, 150000, 150000]
    Et = [15000, 15000, 15000, 15000]
    Glt = [5000, 5000, 5000, 5000]
    Nult = [0.2, 0.2, 0.2, 0.2]
    Cas_1 = Homogeneisation_Mecanique(angles, epaisseurs, El, Et, Glt, Nult)
    """
    
    def __init__ (self, liste_angles, liste_epaisseurs, liste_El, liste_Et, liste_Glt, liste_Nult) :
        """
        Fonction d'initialisation de la Class
        Donnees d'entrees : proprietes materiaux et angles de drapage.
        Les donnees d'entrees sont sous forme de listes.
        """
        self.liste_angles = liste_angles
        self.liste_epaisseurs = liste_epaisseurs
        self.liste_El = liste_El
        self.liste_Et = liste_Et
        self.liste_Glt = liste_Glt
        self.liste_Nult = liste_Nult
        self.nombre_plis = len(liste_angles)
        self.epaisseur_stratifie = sum(liste_epaisseurs)
        self.altitude_inf_pli = self.altitudes()[0]
        self.altitude_sup_pli = self.altitudes()[1]
        self.liste_J = self.liste_J()
        self.liste_Q0  = self.liste_Q0()
        self.liste_Qx  = self.liste_Qx()
    
    def altitudes (self):
        """
        Fonction permettant de calculer les altitudes inferieures et
        superieures des plis du stratifie.
        La fonction retourne une liste comprenant l'ensemble des altitudes inferieures
        ainsi qu'une liste comprenant l'ensemble des altitudes superieures.
        """
        liste_altitude_inf_pli = [-1 * self.epaisseur_stratifie / 2]
        liste_altitude_sup_pli = [liste_altitude_inf_pli[0] + self.liste_epaisseurs[0]]
        for pli in range(len(self.liste_epaisseurs) - 1) :
            altitude_inf_pli = liste_altitude_sup_pli[pli]
            altitude_sup_pli = altitude_inf_pli + self.liste_epaisseurs[pli + 1]
            liste_altitude_inf_pli.append(altitude_inf_pli)
            liste_altitude_sup_pli.append(altitude_sup_pli)
        return liste_altitude_inf_pli, liste_altitude_sup_pli
        
    
    def J (self, angle):
        """
        Fonction permettant le calcul de  :
            la matrice de changement de repere
            l'inverse de la matrice de changement de repere
            la transposee de la matrice de changement de repere
            la transposee de l'inverse de la matrice de changement de repere
        ------
        Entree : angle [float]
        """
        angle = np.radians(float(angle))
        J11 = np.cos(angle)**2
        J12 = np.sin(angle)**2
        J13 = 2 * np.cos(angle) * np.sin(angle)
        J21 = np.sin(angle)**2
        J22 = np.cos(angle)**2
        J23 = -2 * np.cos(angle) * np.sin(angle)
        J31 = -1 * np.cos(angle) * np.sin(angle)
        J32 = np.cos(angle) * np.sin(angle)
        J33 = np.cos(angle)**2 - np.sin(angle)**2
        matrice_J = np.array([[J11, J12, J13],[J21, J22, J23],[J31, J32, J33]])
        matrice_J_inv = np.linalg.inv(matrice_J)
        matrice_J_trans = matrice_J.transpose()
        matrice_J_trans_inv = np.linalg.inv(matrice_J_trans)
        return matrice_J, matrice_J_inv, matrice_J_trans, matrice_J_trans_inv
    
    def Q0 (self, El, Et, Glt, Nult):
        """
        Fonction permettant le calcul de la matrice de rigidite du pli
        dans le repere local L, T.
        -------
        Entrees : 
            El [float]
            Et [float]
            Glt [float]
            Nult [float]
        """
        El = float(El)
        Et = float(Et)
        Glt = float(Glt)
        Nult = float(Nult)
        Nutl = float(Nult * Et / El)
        Q11 = El / (1 - Nult * Nutl)
        Q12 = (Nult * Et) / (1 - Nult * Nutl)
        Q13 = 0
        Q21 = (Nutl * El) / (1 - Nult * Nutl)
        Q22 = Et / (1 - Nult * Nutl)
        Q23 = 0
        Q31 = 0
        Q32 = 0
        Q33 = Glt
        matrice_Q0 = np.array([[Q11, Q12, Q13],[Q21, Q22, Q23],[Q31, Q32, Q33]])
        return matrice_Q0
    
    def Qx (self, matrice_J_inv, matrice_J_trans_inv, matrice_Q0):
        """
        Fonction permettant le calcul de la matrice de rigidite du pli
        dans le repere global x, y.
        -------
        Entrees :
            matrice_J_trans [array]
            matrice_J [array]
            matrice_Q0 [array]
        """
        matrice_Qx = np.dot(np.dot(matrice_J_inv, matrice_Q0), matrice_J_trans_inv)
        return matrice_Qx
    
    def A (self, Qx, altitude_inf_pli, altitude_sup_pli):
        """
        Fonction permettant le calcul de la matrice A du pli.
        -------
        Entrees :
            Qx [array] 
            altitude_inf_pli [float]
            altitude_sup_pli [float]
        """
        A11 = Qx[0][0] * (altitude_sup_pli - altitude_inf_pli)
        A12 = Qx[0][1] * (altitude_sup_pli - altitude_inf_pli)
        A13 = Qx[0][2] * (altitude_sup_pli - altitude_inf_pli)
        A21 = Qx[1][0] * (altitude_sup_pli - altitude_inf_pli)
        A22 = Qx[1][1] * (altitude_sup_pli - altitude_inf_pli)
        A23 = Qx[1][2] * (altitude_sup_pli - altitude_inf_pli)
        A31 = Qx[2][0] * (altitude_sup_pli - altitude_inf_pli)
        A32 = Qx[2][1] * (altitude_sup_pli - altitude_inf_pli)
        A33 = Qx[2][2] * (altitude_sup_pli - altitude_inf_pli)
        matrice_A = np.array([[A11, A12, A13], [A21, A22, A23], [A31, A32, A33]])
        return matrice_A
    
    def B (self, Qx, altitude_inf_pli, altitude_sup_pli):
        """
        Fonction permettant le calcul de la matrice B.
        -------
        Entrees :
            matrice_J_trans [array]
            matrice_J [array]
            matrice_Q0 [array]
        """
        B11 = 0.5 * Qx[0][0] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        B12 = 0.5 * Qx[0][1] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        B13 = 0.5 * Qx[0][2] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        B21 = 0.5 * Qx[1][0] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        B22 = 0.5 * Qx[1][1] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        B23 = 0.5 * Qx[1][2] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        B31 = 0.5 * Qx[2][0] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        B32 = 0.5 * Qx[2][1] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        B33 = 0.5 * Qx[2][2] * (altitude_sup_pli**2 - altitude_inf_pli**2)
        matrice_B = np.array([[B11, B12, B13], [B21, B22, B23], [B31, B32, B33]])
        return matrice_B
    
    def D (self, Qx, altitude_inf_pli, altitude_sup_pli):
        """
        Fonction permettant le calcul de la matrice D.
        -------
        Entrees :
            matrice_J_trans [array]
            matrice_J [array]
            matrice_Q0 [array]
        """
        D11 = (1/3) * Qx[0][0] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        D12 = (1/3) * Qx[0][1] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        D13 = (1/3) * Qx[0][2] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        D21 = (1/3) * Qx[1][0] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        D22 = (1/3) * Qx[1][1] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        D23 = (1/3) * Qx[1][2] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        D31 = (1/3) * Qx[2][0] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        D32 = (1/3) * Qx[2][1] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        D33 = (1/3) * Qx[2][2] * (altitude_sup_pli**3 - altitude_inf_pli**3)
        matrice_D = np.array([[D11, D12, D13], [D21, D22, D23], [D31, D32, D33]])
        return matrice_D
    
    def liste_Q0 (self):
        """
        Fonction permettant de stocker l'ensemble des matrices Q0 dans une liste.
        """
        liste_matrices_Q0 = [self.Q0(self.liste_El[pli], self.liste_Et[pli],
                                        self.liste_Glt[pli], self.liste_Nult[pli]) 
                                        for pli in range(len(self.liste_angles))]
        return liste_matrices_Q0
    
    def liste_J (self):
        """
        Fonction permettant de stocker l'ensemble des matrices 
        J, J_inv, J_trans_inv dans une liste.
        """
        liste_matrices_J = {}
        for angle in range(-90,91,1):
            if angle in self.liste_angles:
                liste_matrices_J[angle] = self.J(angle)
            else:
                liste_matrices_J[angle] = 0
        return liste_matrices_J
    
    def liste_Qx (self):
        """
        Fonction permettant de stocker l'ensemble des matrices Qx dans une liste.
        """
        liste_matrices_Qx = [self.Qx(self.liste_J[self.liste_angles[pli]][1],
                                     self.liste_J[self.liste_angles[pli]][3],
                                     self.liste_Q0[pli]) 
                                    for pli in range(len(self.liste_angles))]
        return liste_matrices_Qx
    
    def liste_A (self):
        """
        Fonction permettant de stocker l'ensemble des matrices A dans une liste.
        """
        liste_matrices_A = [self.A(self.liste_Qx[pli], 
                                   self.altitude_inf_pli[pli], 
                                   self.altitude_sup_pli[pli])
                                    for pli in range(len(self.liste_angles))]
        return liste_matrices_A
    
    def liste_B (self):
        """
        Fonction permettant de stocker l'ensemble des matrices B dans une liste.
        """
        liste_matrices_B = [self.B(self.liste_Qx[pli], 
                                   self.altitude_inf_pli[pli], 
                                   self.altitude_sup_pli[pli])
                                    for pli in range(len(self.liste_angles))]
        return liste_matrices_B
    
    def liste_D (self):
        """
        Fonction permettant de stocker l'ensemble des matrices D dans une liste.
        """
        liste_matrices_D = [self.D(self.liste_Qx[pli], 
                                   self.altitude_inf_pli[pli], 
                                   self.altitude_sup_pli[pli])
                                    for pli in range(len(self.liste_angles))]
        return liste_matrices_D
    
    def ABD (self):
        """
        Fonction permettant de calculer la matrice ABD correspondant
        a l'empilement du stratifie etudie.
        """
        A_globale = sum(self.liste_A())
        B_globale = sum(self.liste_B())
        D_globale = sum(self.liste_D())
        AB = np.concatenate((A_globale, B_globale), axis = 0)
        BD = np.concatenate((B_globale, D_globale), axis = 0)
        matrice_ABD = np.concatenate((AB, BD), axis = 1)
        return matrice_ABD
    
    def proprietes_mecaniques_stratifie (self):
        """
        Fonction permettant de calculer les proprietes homogeneises du
        stratifie composite.
        """
        matrice_A_globale = sum(self.liste_A()) / self.epaisseur_stratifie
        matrice_A_globale_inv = np.linalg.inv(matrice_A_globale)
        Ex_stratifie = 1 / (matrice_A_globale_inv[0][0])
        Ey_stratifie = 1 / (matrice_A_globale_inv[1][1])
        Gxy_stratifie = 1 / (matrice_A_globale_inv[2][2])
        Nuxy_stratifie = -Ex_stratifie * matrice_A_globale_inv[1][0]
        return Ex_stratifie, Ey_stratifie, Gxy_stratifie, Nuxy_stratifie


#
#
#liste_angles = [0, 90, 90, 0]
#liste_epaisseurs = [0.25, 0.25, 0.25, 0.25]
#liste_El = [45600, 45600, 45600, 45600]
#liste_Et = [16200, 16200, 16200, 16200]
#liste_Glt = [5830, 5830, 5830, 5830]
#liste_Nult = [0.278, 0.278, 0.278, 0.278]
#
#test_1 = Homogeneisation_Mecanique_Stratifie(liste_angles, liste_epaisseurs, liste_El, liste_Et, liste_Glt, liste_Nult)