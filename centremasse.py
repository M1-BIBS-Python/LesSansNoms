#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Ce programme calcule le centre de masse d'une protéine
# Le résultat sera utilisé pour caculer le rayon de gyration de la protéine


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

    dicPDB["listChains"] = []     # le dictionnaire a la cle "listChains" qui prend une liste


    for line in lines:
        # Pour toutes les lignes qui commencent par ATOM (celles qui ont des atomes)
        if (line[0:4] == "ATOM"):                  # on ne selectionne que les lignes qui contiennent des ATOM

            chain = line[21]

            if not chain in dicPDB["listChains"]:
                dicPDB["listChains"].append(chain)
                dicPDB[chain] = {}                  # Creation d'une cle chain dans le dictionnaire

            dicPDB[chain]["reslist"] = []           # qui a une cle reslist prenant une liste comme valeur
            number = "%s"%(line[22:26]).strip()     # numero du residu

            if not number in dicPDB[chain]["reslist"] :
                dicPDB[chain]["reslist"].append(number)
                dicPDB[chain][number] = {}          # Creation d'un dictionnaire dans dPBD[chain]
                                                    # pour la cle number ayant pour cle "resname"
                dicPDB[chain][number]["resname"] = string.strip(line[17:20])

            dicPDB[chain][number]["atomlist"] = []  # a pour cle atomlist et prend une liste

            atomtype = string.strip(line[13:16])

            dicPDB[chain][number]["atomlist"].append(atomtype) # ajout de l'atome a la liste

            dicPDB[chain][number][atomtype] = {}    # cree un dictionnaire dans dicPBD[chain][number]
                                                       # pour la cle "atomtype"

            dicPDB[chain][number][atomtype]["x"] = float(line[30:38])
            dicPDB[chain][number][atomtype]["y"] = float(line[38:46])
            dicPDB[chain][number][atomtype]["z"] = float(line[46:54])
            dicPDB[chain][number][atomtype]["id"] = line[6:11].strip()

    # Test presence d'ATOM
    if dicPDB==0:
        print "Le fichier ne contient pas d'ATOM"
        sys.exit(1)

    # Fermeture du fichier
    f.close()

    # Affichage
    for i in dicPDB.keys():
        print dicPDB[i]

    return dicPDB


def centredemasse(dicopdb):
    dico_atomes = dict()
    dico_atomes["C"]=12
    dico_atomes["H"]=1
    dico_atomes["0"]=16
    dico_atomes["N"]=14
    dico_atomes["S"]=32
    dict_coord = {}
    m = 0
    x=0
    y=0
    z=0

    for atomtype in dicopdb.keys():
        if atomtype != 'resname' and atomtype != 'atomlist':
            for key in dico_atomes.keys():
                mass = dico_atomes[key]
            x += mass * dico_atomes[atomtype]['x']
            y += mass * dico_atomes[atomtype]['y']
            z += mass * dico_atomes[atomtype]['z']
            m += mass

    dict_coord['x']= x/m
    dict_coord['y']= y/m
    dict_coord['z']= z/m

    return dict_coord


dic = parser_pdb()
x = centredemasse(dic)
print x
