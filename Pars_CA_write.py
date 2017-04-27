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



def write2PDB(dPDB, filout = "PDB2_out.pdb") :
    """ Ecriture des coordonnees de chaque atome dans un fichier .pdb
        Input : PDB dictionnary
        Output : file containing the coordinates
    """

    fout = open(filout, "w")

    for confo in dPDB["listChains"]:
        for res in dPDB[confo]["reslist"]:
            fout.write("RMSD        %8.3f\n"%(dPDB[confo][res]["rmsd"]))
            for atom in dPDB[confo][res]["atomlist"] :
                fout.write("ATOM  %4s      %4s     %8.3f  %8.3f  %8.3f\n"%(confo, res,dPDB[confo][res][atom]["x"], dPDB[confo][res][atom]["y"],dPDB[confo][res][atom]["z"]))
                    #      "ATOM" confo     RES       X     Y      Z
    fout.close()
