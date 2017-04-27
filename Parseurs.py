#!/usr/bin/python
# @author : Alexandra Benamar, Robin Droit
# date: 06/03/2017

import os                      # gestion de fichiers et de dossiers
import sys                     # gestion des erreurs et des arguments
import string


def pars_RMSD(dPDB):
    """ Cette fonction prend le dictionnaire obtenu en parsant
    le fichier PDB d'origine et renvoie un nouveau dictionnaire
    qui prend uniquement en compte les carbones alpha

    Input : dictionnaire PDB.
    Output : dictionnaire carbones alpha.

    """

    dRMSD = {}                         # initialisation du dictionnaire RMSD
    compteur=0
    list_CA=[]

    for i in dPDB.keys():              # On parcourt le dico dans sa longueur
        list_CA = dPDB[i]              # on recupere le dic par residus
        new_list = []
        for j in list_CA:
            if dicPDB[modelnumber][chain][number]["atomlist"] == "CA"
                new_list.append(list_CA[j*6 + 1])
                new_list.append(list_CA[j*6 + 2])
                new_list.append(list_CA[j*6 + 3])
                new_list.append(list_CA[j*6 + 4])
                new_list.append(list_CA[j*6 + 5])
        dRMSD[compteur] = new_list
        compteur += 1

    return dRMSD
