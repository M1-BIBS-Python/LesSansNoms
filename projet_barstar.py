#!/usr/bin/python
# Author : Alexandra Benamar, Robin Droit

import os                      # gestion de fichiers et de dossiers
import sys                     # gestion des erreurs et des arguments
import string
import math
import numpy as np

#------------------------------------------------------------------
# Fonction permettant de parser le fichier

def parserPDB(infile):
    """ Cette fonction a pour but de parser un fichier de type pdb afin
    d'en recuperer les informations sur les differents atomes qui composent
    les acides amines.

    Input : fichier pdb a parser.
    Output : dictionnaire contenant les informations sur le fichier.

    """

    # Test d'ouverture du fichier
    try:
		f = open(infile,'r')
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
    # Initialisation d'un flag pour remplir la liste des atomes
    Flag = False

    for line in lines:

		if line[0:5] == "MODEL":
			modelnumber = int(string.strip(line[10:14]))
			dicPDB[modelnumber] = {}
			modelnumber = int(modelnumber)
			dicPDB[modelnumber]["listChains"] = []  # numero du residu
			dicPDB[modelnumber]["listRes"] = []     # noms des residus
			# le dictionnaire a la cle "listChains" qui prend une liste

												   # Pour toutes les lignes qui commencent par ATOM (celles qui ont des atomes)
		elif line[0:4] == "ATOM":
			chain = line[24:27]

												   # on ne selectionne que les lignes qui contiennent des ATOM
			if chain not in dicPDB[modelnumber]["listChains"]:
				dicPDB[modelnumber]["listChains"].append(chain)
				dicPDB[modelnumber]["listRes"].append(string.strip(line[17:20]))
				dicPDB[modelnumber][chain] = {}
                                                    # pour la cle number ayant pour cle "resname"
			if dicPDB[modelnumber][chain].has_key("resname") == False:
				dicPDB[modelnumber][chain]["resname"] = string.strip(line[17:20])

				dicPDB[modelnumber][chain]["atomlist"] = []  # a pour cle atomlist et prend une liste

			atomtype = string.strip(line[13:16])
			#print atomtype

			dicPDB[modelnumber][chain]["atomlist"].append(atomtype) # ajout de l'atome a la liste

			dicPDB[modelnumber][chain][atomtype] = {}    # cree un dictionnaire dans dicPBD[chain][number]


			dicPDB[modelnumber][chain][atomtype]["x"] = float(line[30:38])
			dicPDB[modelnumber][chain][atomtype]["y"] = float(line[38:46])
			dicPDB[modelnumber][chain][atomtype]["z"] = float(line[46:54])
			dicPDB[modelnumber][chain][atomtype]["id"] = line[6:11].strip()
			dicPDB[modelnumber][chain][atomtype]["lyst"] = [float(line[30:38]),float(line[38:46]),float(line[46:54]),line[6:11].strip()]



	# Test presence d'ATOM
    if dicPDB==0:
		print "Le fichier ne contient pas d'ATOM"
		sys.exit(1)

    # Fermeture du fichier
    f.close()
    return dicPDB


#------------------------------------------------------------------
# Fontions concernant un parseur avec carbones alpha uniquement

def parserPDB_CA(infile) :
    """ Cette fonction a pour but de parser un fichier de type pdb afin
    d'en recuperer les informations sur les carbones alpha.

    Input : fichier pdb a parser.
    Output : dictionnaire contenant les informations sur les carbones alpha.

    """

    # Test d'ouverture du fichier
    try:
		f = open(infile,'r')
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
    # Initialisation d'un flag pour remplir la liste des atomes
    Flag = False

    for line in lines:

		if line[0:5] == "MODEL":
			modelnumber = string.strip(line[10:14])
			dicPDB[modelnumber] = {}
			dicPDB[modelnumber]["listChains"] = []  # numero du residu
			dicPDB[modelnumber]["listRes"] = []     # noms des residus
			# le dictionnaire a la cle "listChains" qui prend une liste

												   # Pour toutes les lignes qui commencent par ATOM (celles qui ont des atomes)
		elif line[0:4] == "ATOM":
			chain = line[24:27]
												   # on ne selectionne que les lignes qui contiennent des ATOM
			if chain not in dicPDB[modelnumber]["listChains"]:
				dicPDB[modelnumber]["listChains"].append(chain)
				dicPDB[modelnumber]["listRes"].append(string.strip(line[17:20]))
				dicPDB[modelnumber][chain] = {}
                                                    # pour la cle number ayant pour cle "resname"
			if dicPDB[modelnumber][chain].has_key("resname") == False:
				dicPDB[modelnumber][chain]["resname"] = string.strip(line[17:20])
				dicPDB[modelnumber][chain]["atomlist"] = []  # a pour cle atomlist et prend une liste

			atomtype = string.strip(line[13:16])
			dicPDB[modelnumber][chain]["atomlist"].append(atomtype) # ajout de l'atome a la liste
			dicPDB[modelnumber][chain][atomtype] = {}    # cree un dictionnaire dans dicPBD[chain][number]

			if atomtype == "CA": #on conserve dans une autre partie du dictionnaire les carbones alpha
				dicPDB[modelnumber][chain]["CA"] = {}
				dicPDB[modelnumber][chain]["CA"]["x"] = float(line[30:38])
				dicPDB[modelnumber][chain]["CA"]["y"] = float(line[38:46])
				dicPDB[modelnumber][chain]["CA"]["z"] = float(line[46:54])
				dicPDB[modelnumber][chain]["CA"]["id"] = line[6:11].strip()
				dicPDB[modelnumber][chain]["CA"]["lyst"] = [float(line[30:38]),float(line[38:46]),float(line[46:54]),line[6:11].strip()]

	# Test presence d'ATOM
    if dicPDB==0:
		print "Le fichier ne contient pas d'ATOM"
		sys.exit(1)

    # Fermeture du fichier
    f.close()

    return dicPDB



def writePDB_CA(dPDB, filout = "PDB_CA_out.pdb") :
    """ Ecriture des coordonnees de chaque atome dans un fichier .pdb
        Input : dictionnaire PDB
        Output : fichier contenant les coordonnees
    """

    fout = open(filout, "w")

    for confo in dPDB["listChains"]:
        for res in dPDB[confo]["listRes"]:
            fout.write("RMSD        %8.3f\n"%(dPDB[confo][res]["rmsd"]))
            for atom in dPDB[confo][res]["atomlist"] :
                fout.write("ATOM  %4s      %4s     %8.3f  %8.3f  %8.3f\n"%(confo, res,dPDB[confo][res][atom]["x"], dPDB[confo][res][atom]["y"],dPDB[confo][res][atom]["z"]))
                    #      "ATOM" confo     RES       X     Y      Z
    fout.close()


#------------------------------------------------------------------
# RMSD global et RMSD local

def rmsd_global(dPDB):
    """ Calcul du RMSD global de chaque modele avec uniquement les carbones alpha
        Output : RMSD entre le modele pris en compte et le modele de reference
    """
    for model in range(0,len(dPDB)):
        n=somme=0
        conflist = dPDB[model]["listChains"]
        for confo in conflist:
            atomlist = dPDB[model][confo]["atomlist"]
            for atom in atomlist:
                distance=(dPDB[model][confo][atom]["x"]-dPDB[0][confo][atom]["x"])**2
                +(dPDB[model][confo][atom]["y"]-dPDB[0][confo][atom]["y"])**2
                +(dPDB[model][confo][atom]["z"]-dPDB[0][confo][atom]["z"])**2
                somme+=distance
                n+=1

        dPDB[model]["rmsd"] = math.sqrt(somme/n)



def rmsd_local(dPDB): #on donne soit un dico specifique soit le dico general et les numeros des modeles et on renvoi une liste de tous les rmsd de tous les residus
    """
        Permet de calculer l'ensemble des distances entre les differents carbones alpha
        entre deux modeles.
        Input : numeros des modeles a comparer et un dictionnaire. Le dictionnaire peut etre specifique
        ou un dictionnaire general.
        Output : liste des RMSD de tous les residus.
    """

    for model in range(0,len(dPDB)):
        conflist = dPDB[model]["listChains"]
        for confo in conflist:
            n=somme=0
            atomlist = dPDB[model][confo]["atomlist"]
            for atom in atomlist:
                distance=(dPDB[model][confo][atom]["x"]-dPDB[0][confo][atom]["x"])**2
                +(dPDB[model][confo][atom]["y"]-dPDB[0][confo][atom]["y"])**2
                +(dPDB[model][confo][atom]["z"]-dPDB[0][confo][atom]["z"])**2
                somme+=distance
                n+=1

            dPDB[model][confo]["rmsd"] = math.sqrt(somme/n)

#------------------------------------------------------------------
# Rayon de giration

def rayon_giration(dPDB) :
    """ Calcule le rayon de giration de chaque conformation du dictionnaire donne en entree
        Rayon de giration : distance, pour chaque conformation, entre son centre de masse
        et l'atome le plus eloigne de ce centre de masse.
    """
    for model in range(0,len(dPDB)):
        conflist = dPDB[model]["listChains"]
        raymax=0
        for conf in conflist:
            atomlist = dPDB[model][conf]["atomlist"]
            for atom in atomlist:
                dist = math.sqrt((dPDB[model][conf][atom]["x"]-dPDB[model]["XCM"])**2
                +(dPDB[model][conf][atom]["y"]-dPDB[model]["YCM"])**2
                +(dPDB[model][conf][atom]["z"]-dPDB[model]["ZCM"])**2)
                if (dist>=raymax):
                    raymax = dist
        dPDB[model]["giration"] = raymax

#------------------------------------------------------------------
# Distance entre le centre de masse d'un residu et la conformation
def distance(dPDB) :
    """
        Calcul de la distance entre le centre de masse d'un residu et la conformation
    """
    for model in range(0,len(dPDB)):
        confo_list = dPDB[model]["listChains"]
        for confo in confo_list:
            distance = math.sqrt((dPDB[model][confo]["XCM"]-dPDB[model]["XCM"])**2
            +(dPDB[model][confo]["YCM"]-dPDB[model]["YCM"])**2
            +(dPDB[model][confo]["ZCM"]-dPDB[model]["ZCM"])**2)
            dPDB[model][confo]["dist_CM"] = distance

#------------------------------------------------------------------
# Centre de masse

def centerMassOfConf(dPDB):
    """ Calcul du centre de masse de chaque conformation.
        Tient compte de tous les atomes de chaque residus.

        Input : dictionnaire
        Output : centre de masse.
    """

    for model in range(0,len(dPDB)):
        x=y=z=cpt=0
        confo_list = dPDB[model]["listChains"]
        for confo in confo_list :
            atomlt = dPDB[model][confo]["atomlist"]
            for atom in atomlt:
                x += dPDB[model][confo][atom]['x']
                y += dPDB[model][confo][atom]['y']
                z += dPDB[model][confo][atom]['z']
                cpt += 1
        Xcm = float(x)/cpt
        Ycm = float(y)/cpt
        Zcm = float(z)/cpt
        dPDB[model]["XCM"] = Xcm
        dPDB[model]["YCM"] = Ycm
        dPDB[model]["ZCM"] = Zcm


def centerMassOfRes(dPDB):
    """ Calcule le centre de masse d'un residu en tenant compte
        de tous ses atomes.
    """

    for model in range(0,len(dPDB)):
        confo_list = dPDB[model]["listChains"]
        for confo in confo_list :
            cpt=x=y=z=0
            atomlt = dPDB[model][confo]["atomlist"]
            for atom in atomlt:
                x += dPDB[model][confo][atom]['x']
                y += dPDB[model][confo][atom]['y']
                z += dPDB[model][confo][atom]['z']
                cpt += 1
            Xcm = float(x)/cpt
            Ycm = float(y)/cpt
            Zcm = float(z)/cpt
            dPDB[model][confo]["XCM"] = Xcm
            dPDB[model][confo]["YCM"] = Ycm
            dPDB[model][confo]["ZCM"] = Zcm

#------------------------------------------------------------------
# Calcul des distances moyennes

def distance_moyenne(dPDB):
    """ Calcul de la distance moyenne au cours du temps entre le
        centre de masse d'un residu et le centre de masse de la proteine
    """
    dist_moy = {}
    dist_moy["reslist"] = []
    n=0

    for model in range(0,len(dPDB)):
        confo_list = dPDB[model]["listChains"]
        for confo in confo_list:
            if not confo in dist_moy["reslist"]:
                dist_moy["reslist"].append(confo)
                dist_moy[confo]=dPDB[model][confo]["dist_CM"]
            else:
                dist_moy[confo]+=dPDB[model][confo]["dist_CM"]
        n+=1
    reslist=dist_moy["reslist"]
    for res in reslist:
        dist_moy[res]=dist_moy[res]/n

    return dist_moy


def rmsd_moyen(dPDB):
    """
        Calcul du RMSD local moyen au cours du temps pour chaque acide amine.
    """
    rmsd_moy={}
    rmsd_moy["reslist"]=[]
    n=0

    for model in range(0,len(dPDB)):
        confo_list = dPDB[model]["listChains"]
        for conf in confo_list:
            if not conf in rmsd_moy["reslist"]:
                rmsd_moy["reslist"].append(conf)
                rmsd_moy[conf]=dPDB[model][conf]["rmsd"]
            else:
                rmsd_moy[conf]+=dPDB[model][conf]["rmsd"]
        n+=1
    reslist=rmsd_moy["reslist"]
    for res in reslist:
        rmsd_moy[res]=rmsd_moy[res]/n

    return rmsd_moy

#------------------------------------------------------------------
# Ecriture des resultats dans des fichiers de sortie

def print_PDB(dPDB, filout="output_pdb.txt"):
    """
        Ecriture des coordonnees de chaque atome dans un fichier .txt
    """
    fout = open(filout, "w")

    for model in range(0,len(dPDB)):
        fout.write("CEMA  %4s                          %8.3f%8.3f%8.3f\n"
        %(model, dPDB[model]["XCM"],dPDB[model]["YCM"],dPDB[model]["ZCM"]))
        fout.write("RAGI  %4s                                                   %8.3f\n"
        %(model, dPDB[model]["giration"]))
        for conf in dPDB[model]["listChains"]:
            fout.write("CMRE   %4s                  %3s     %8.3f%8.3f%8.3f\n"
            %(model, conf, dPDB[model][conf]["XCM"], dPDB[model][conf]["YCM"], dPDB[model][conf]["ZCM"]))
            fout.write("D_CM   %4s                  %3s                     %8.3f\n"
            %(model, conf, dPDB[model][conf]["dist_CM"]))
            fout.write("RMSD   %4s                  %3s                     %8.3f\n"
            %(model, conf, dPDB[model][conf]["rmsd"]))
            for atom in dPDB[model][conf]["atomlist"]:
                fout.write("ATOM  %4s       %4s      %4s     %8.3f%8.3f%8.3f\n"
                %(model, conf, atom, dPDB[model][conf][atom]["x"], dPDB[model][conf][atom]["y"], dPDB[model][conf][atom]["z"]))
    fout.close()

def print_PDB_CA(dPDB, filout="output_pdb.txt"):
    """
        Ecriture des coordonnees de chaque CA dans un fichier .txt
    """
    fout = open(filout, "w")

    for model in range(0,len(dPDB)):
        for conf in dPDB[model]["listChains"]:
            fout.write("RMSD   %4s                  %3s                     %8.3f\n"
            %(model, conf, dPDB[model][conf]["rmsd"]))
            for atom in dPDB[model][conf]["atomlist"]:
                fout.write("ATOM  %4s       %4s      %4s     %8.3f%8.3f%8.3f\n"
                %(model, conf, atom, dPDB[model][conf][atom]["x"], dPDB[model][conf][atom]["y"], dPDB[model][conf][atom]["z"]))
    fout.close()


def print_Gira(dPDB, fileout = "giration.pdb"):

    fout = open(fileout, "w")
    for model in range(0,len(dPDB)):
        fout.write("CONFO   %4s   RAYON_GIRATION   %8.3f\n"
        %(model, dPDB[model]["giration"]))
    fout.close()


def print_DMoy(dist_moy, fileout="dist_moy.pdb"):

    fout=open(fileout,"w")
    reslist=dist_moy["reslist"]
    for res in reslist:
        fout.write("RES   %4s   DIST_MOY   %8.3f\n"
        %(res, dist_moy[res]))
    fout.close()

def print_rmsd_moy(dRMSD, fileout="dist_moy.pdb"):

    fout=open(fileout,"w")
    reslist=dRMSD["reslist"]
    for res in reslist:
        fout.write("RES   %4s   RMSD_MOY   %8.3f\n"
        %(res, dRMSD[res]))
    fout.close()


#------------------------------------------------------------------
# Informations sur le projet

def usage() :
    """
        Affichage des informations concernant le projet et son utilisation
    """
    print("\n\nProjet d'Alexandra BENAMAR et Robin DROIT\nCe programme permet de parser un fichier .pdb\net calcule les distances globales et locales de la proteine.\nUsage : ./main.py -pdb <nom du fichier PDB>\n")

#------------------------------------------------------------------
# MAIN

###############################################
###        OUVERTURE DU FICHIER .PDB        ###
###############################################

try:
    # on veut que le fichier .pdb soit en argument juste apres -pdb
    pdb_file = sys.argv[1]
    # affichage
    print "Fichier pdb: ", pdb_file

except:
    usage()
    print "ERROR: please, enter the name of the pdb input"
    sys.exit()


dico_1 = parserPDB(pdb_file)
centerMassOfConf(dico_1)
centerMassOfRes(dico_1)
distance(dico_1)
rayon_giration(dico_1)
rmsd_local(dico_1)
dist_moy = distance_moyenne(dico_1)
rmsd_moy = rmsd_moyen(dico_1)
print_PDB(dico_1)
print_Gira(dico_1)
print_DMoy(dist_moy)
print_rmsd_moy(rmsd_moy)


dico_2 = parserPDB_CA(pdb_file)
rmsd_global(dico_2)
rmsd_local(dico_2)
print_PDB_CA(dico_2)
print_PDB(dico_2)
#print CenterOfMassConf(dico_1)
#print Giration_local()
