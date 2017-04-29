#!/usr/bin/python2.7
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
    # Ce parseur (qui fonctionne de la meme maniere que le fichier parse.pdb) permet de parser le fichier pdb en ne retenant que les atomes du
    # squelette peptidique de chaque residu (CA, C, O, N)

        # lecture du fichier PDB
        f = open(infile, "r")
        lines = f.readlines()
        f.close()

        # var init
        conf= True
        firstline = True
        prevres = None
        dPDB = {}
        dPDB["conformation"] = []
        dPDB["liste_residus"] = []

        # parcoure le PDB
        for line in lines :
                if line[0:5] == "MODEL" :
                        conf=int(line[6:14])
                        if not conf in dPDB["conformation"] :
                                dPDB["conformation"].append(conf)
                                dPDB[conf] = {}
                                dPDB[conf]["liste_residus"] = []
                elif line[0:4] == "ATOM" :
                    if line[13:15] == "CA":
                        # On ne retient que les atomes du squelette peptidique
                        curres = "%s"%(line[22:26]).strip()
                        if not curres in dPDB[conf]["liste_residus"] :
                                dPDB[conf]["liste_residus"].append(curres)
                                dPDB[conf][curres] = {}
                                dPDB[conf][curres]["nom_resid"] = string.strip(line[17:20])
                                dPDB[conf][curres]["liste_atomes"] = []

                        atomtype = string.strip(line[12:16])
                        dPDB[conf][curres]["liste_atomes"].append(atomtype)
                        dPDB[conf][curres][atomtype] = {}
                        dPDB[conf][curres][atomtype]["x"] = float(line[30:38])
                        dPDB[conf][curres][atomtype]["y"] = float(line[38:46])
                        dPDB[conf][curres][atomtype]["z"] = float(line[46:54])


        return dPDB




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


#------------------------------------------------------------------
# RMSD global et RMSD local

def RMSD_global(model1, model2):
	longueur = (len(model1))/5
	distance = []
	somme=0.0
	for i in range(longueur):
		distance.append((model1[i*5+2] - model2[i*5+2])**2 + (model1[i*5+3] - model2[i*5+3])**2 + (model1[i*5+4]-model2[i*5+4])**2)
	for i in distance:
		somme += i
	RMSD = math.sqrt(((1/float(longueur))*somme))
	return RMSD

def RMSD_local(dico,model1,model2): #on donne soit un dico specifique soit le dico general et les numeros des modeles et on renvoi une liste de tous les rmsd de tous les residus
    """
        Permet de calculer l'ensemble des distances entre les differents carbones alpha
        entre deux modeles.
        Input : numeros des modeles a comparer et un dictionnaire. Le dictionnaire peut etre specifique
        ou un dictionnaire general.
        Output : liste des RMSD de tous les residus.
    """

    liste_rmsd = []
    for i in dico[model1]["listChains"]:
        rmsd = sqrt((dico[model2][i]["CA"]["x"] - dico[model1][i]["CA"]["x"])**2
        + (dico[model2][i]["CA"]["y"] - dico[model1][i]["CA"]["y"])**2
        + (dico[model2][i]["CA"]["z"] - dico[model1][i]["CA"]["z"])**2)
        liste_rmsd.append(rmsd)

    return liste_rmsd

#------------------------------------------------------------------
# Rayons de giration global et local

def Giration_redidu(centre,residu):
    """
        Cette fonction nous donne la distance entre deux points dans des reperes en trois dimensions
        applique au centre et a un autre residu nous obtenons un rayon de giration
        on donne en entree deux listes correspondantes aux coordonnees du centre et du residu considere
    """
    return sqrt((residu[0]-centre[0])**2 + (residu[1]-centre[1])**2 + (residu[2]-centre[2])**2)

#------------------------------------------------------------------

def distanceMoy(dPDB) :
    # La fonction distanceMoy permet de calculer la distance moyenne au cours du temps entre le centre de masse d'un residu et le centre de masse de la proteine

    dMEAN = {}
    dMEAN["liste_residus"] = []
    n = 0
    conflist = dPDB["conformation"]
    for conf in conflist :
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist :
            if not res in dMEAN["liste_residus"] :
                dMEAN["liste_residus"].append(res)
                dMEAN[res] = dPDB[conf][res]["dist_CM"]
            else :
                dMEAN[res] += dPDB[conf][res]["dist_CM"]
        n += 1
    reslist = dMEAN["liste_residus"]
    for res in reslist :
        dMEAN[res] = dMEAN[res]/n

    return dMEAN

def distance(dPDB) :
    """
        Calcul de la distance entre le centre de masse d'un residu et la conformation
    """
    confo_list = dPDB["listChains"]
    for confo in confo_list :
        reslist = dPDB[confo]["reslist"]
        for res in reslist :
            dist = sqrt((dPDB[conf][res]["XCM"]-dPDB[conf]["XCM"])**2 +(dPDB[conf][res]["YCM"]-dPDB[conf]["YCM"])**2 +(dPDB[conf][res]["ZCM"]-dPDB[conf]["ZCM"])**2)
            dPDB[conf][res]["dist_CM"] = dist


def print_Rayon_Giration(dPDB, fileout = "giration.pdb") :

    fout = open(fileout, "w")

    conflist = dPDB["listChains"]
    for conf in conflist :
        fout.write("CONFORMATION  %4s    GIRATION  %8.3f\n" %(conf, dPDB[conf]["rayon_giration"]))
        #          "CONFORMATION" conf  "GIRATION" rayon

    fout.close()

#------------------------------------------------------------------
# Centre de masse

def CenterOfMassConf(dPDB):
    """ Calcul du centre de masse de chaque conformation

    Input : dictionnaire
    Output : centre de masse.
    """

    dico_atomes = {"C":12, "H":1, "O":16, "N":14, "S":32}

    dict_coord = {}
    m = 0
    x=y=z=0

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

#------------------------------------------------------------------


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
    print "pdb a ouvrir: ", pdb_file

except:
    usage()
    print "ERROR: please, enter the name of the pdb input"
    sys.exit()


dico_1 = parserPDB(pdb_file)
dico_2 = parserPDB_CA(pdb_file)
