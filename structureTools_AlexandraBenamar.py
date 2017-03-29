#!/usr/bin/python3
# @author : Alexandra Benamar, Robin Droit
# date: 06/03/2017

import os                      # gestion de fichiers et de dossiers
import sys                     # gestion des erreurs et des arguments
import string


def parser_pdb():
    """ Cette fonction a pour but de parser un fichier de type pdb afin
    d'en recuperer les informations sur les differents atomes qui composent
    les acides amines.

    Input : fichier pdb a parser.
    Output : dictionnaire contenant les informations sur le fichier.

    Usage : ./parser.py <fichier.pdb>

    """

    # Le nom du fichier est passe en argument
    if len(sys.argv) != 2:
        sys.exit("ERREUR : il faut exactement un argument.")

    # Test d'ouverture du fichier
    try:
        f = open(sys.argv[1],'r')
    except:
        print "Erreur, le fichier n'a pas pu s'ouvrir"
        sys.exit(1)

    lines  = f.readlines()
    number_of_lines=len(lines)


    # Test fichier vide
    if number_of_lines==0:
		print "Erreur, le fichier est vide"
		sys.exit(1)

    # Initialisation d'un dictionnaire
    dicPDB = {}

    for line in lines:

		if line[0:5] == "MODEL":
			modelnumber = line[10:14]
			dicPDB[modelnumber] = {}
			dicPDB[modelnumber]["listChains"] = []     # le dictionnaire a la cle "listChains" qui prend une liste

												   # Pour toutes les lignes qui commencent par ATOM (celles qui ont des atomes)
		elif line[0:4] == "ATOM":
			chain = line[21]
												   # on ne selectionne que les lignes qui contiennent des ATOM
			if not chain in dicPDB[modelnumber]["listChains"]:
				dicPDB[modelnumber]["listChains"].append(chain)
				dicPDB[modelnumber][chain] = {}                  # Creation d'une cle chain dans le dictionnaire

			dicPDB[modelnumber][chain]["reslist"] = []           # qui a une cle reslist prenant une liste comme valeur
			number = "%s"%(line[22:26]).strip()     # numero du residu

			if not number in dicPDB[modelnumber][chain]["reslist"]:
				dicPDB[modelnumber][chain]["reslist"].append(number)
				dicPDB[modelnumber][chain][number] = {}          # Creation d'un dictionnaire dans dPBD[chain]
                                                    # pour la cle number ayant pour cle "resname"
				dicPDB[modelnumber][chain][number]["resname"] = string.strip(line[17:20])

			dicPDB[modelnumber][chain][number]["atomlist"] = []  # a pour cle atomlist et prend une liste

			atomtype = string.strip(line[13:16])

			dicPDB[modelnumber][chain][number]["atomlist"].append(atomtype) # ajout de l'atome a la liste

			dicPDB[modelnumber][chain][number][atomtype] = {}    # cree un dictionnaire dans dicPBD[chain][number]
                                                       # pour la cle "atomtype"

			dicPDB[modelnumber][chain][number][atomtype]["x"] = float(line[30:38])
			dicPDB[modelnumber][chain][number][atomtype]["y"] = float(line[38:46])
			dicPDB[modelnumber][chain][number][atomtype]["z"] = float(line[46:54])
			dicPDB[modelnumber][chain][number][atomtype]["id"] = line[6:11].strip()

	# Test presence d'ATOM
    if dicPDB==0:
		print "Le fichier ne contient pas d'ATOM"
		sys.exit(1)

    # Fermeture du fichier
    f.close()

    return dicPDB



def pars_RMSD(dPDB):
    """ Cette fonction prend le dictionnaire obtenu en parsant
    le fichier PDB d'origine et renvoie un nouveau dictionnaire
    qui prend uniquement en compte les carbones alpha

    Input : dictionnaire PDB.
    Output : dictionnaire carbones alpha.

    """

    dRMSD = {}               # initialisation du dictionnaire RMSD
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

    print dRMSD
    return dRMSD

dictionnaire = {}
dictionnaire = parser_pdb()
A=pars_RMSD(dictionnaire)
print A
