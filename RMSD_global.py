
import numpy as np

def RMSD_global(model1, model2):
    """ Calcul du RMSD global entre deux modeles : le modele souhaite
        contre le modele de reference.

        Input : dictionnaire des carbones alpha.
        Output : RMSD (float).
	"""
    D = len(model1[0])
    N = len(model1)
    rmsd = 0.0
    for m1, m2 in zip(model1, model2)     # zip retourne une liste de tuples
        rmsd += sum([(m1[i]-m2[i])**2.0 for i in range(D)])

    return np.sqrt(rmsd/N)    # racine carree sur chaque element de la liste
