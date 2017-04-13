#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Ce programme calcule le centre de masse d'une protéine
# Le résultat sera utilisé pour caculer le rayon de gyration de la protéine


import os                      # gestion de fichiers et de dossiers
import sys                     # gestion des erreurs et des arguments
import string

def CenterOfMass(dPDB):
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
