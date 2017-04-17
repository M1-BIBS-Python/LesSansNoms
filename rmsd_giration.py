#!/usr/bin/python
# @author : Alexandra Benamar, Robin Droit
import sys
import string
import math
import matplotlib.pyplot as plt


def RMSD_local(dico,model1,model2): #on donne soit un dico specifique soit le dico general et les numeros des modeles et on renvoi une liste de tous les rmsd de tous les residus
    """Cette fonction nous permet de calculer l'ensemble des distances entre les differents carbonnes alpha entre deux modeles
       en entree on donne le dictionnaire obtenu par le parseur et les numeros des modeles a comparer"""

    liste_rmsd = []
    for i in dico[model1]["listChains"]:
        liste_rmsd.append(sqrt(dico[(model2][i]["CA"]["x"] - dico[model1][i]["CA"]["x"])**2 + (dico[model2][i]["CA"]["y"] - dico[model1][i]["CA"]["y"])**2 + (dico[model2][i]["CA"]["z"] - dico[model1][i]["CA"]["z"])**2))

    return liste_rmsd


def Giration_redidu(centre,residu): #calcul de la distance entre deux residus on donne les deux listes des coordonnees presentes dans le dico du parseur
    """Cette fonction nous donne la distance entre deux points dans des reperes en trois dimensions
       applique au centre et a un autre residu nous obtenons un rayon de giration
       on donne en entree deux listes correspondantes aux coordonnees du centre et du residu considere"""

    return sqrt((residu[0]-centre[0])**2 + (residu[1]-centre[1])**2 + (residu[2]-centre[2])**2)
