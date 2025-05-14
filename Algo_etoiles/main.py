from PIL import Image
from math import *
import csv
import random

import config
import display
import data
import geometry
import identification
import cmlcm
import kdtree

#ALGORITHMES: 
# https://www.mdpi.com/1424-8220/20/11/3027 MVDT-SI
# https://www.mdpi.com/2304-6732/9/1/13 CMLCM

#Traitement image:

img_original = Image.open(config.IMAGE_PATH)
img_grayscale = img_original.convert('L')

img_size = img_width, img_height = img_original.size

img_nb =  (cmlcm.cmlcmtotal(img_grayscale, img_width, img_height, config.BLOCKSIZE, config.S, config.L) if config.CMLCM_OR_NOT
           else img_grayscale.point(lambda px: 1 if px>config.BLACK_WHITE_THRESHOLD else 0, mode='1'))

img_nb.show()

LISTE_ETOILES_IMAGE = cmlcm.star_list(img_nb, img_width, img_height)

#Import BDD

DATA_BASE = data.parse_database_file(config.DATA_BASE_TREATED_PATH)
CON_DIC = data.parse_constellation_file(config.CON_DIC_PATH)

#Boucles principales

#LISTE_ETOILES_REF = data.choose_random(LISTE_ETOILES_IMAGE, LAMBD)
LISTE_ETOILES_REF = LISTE_ETOILES_IMAGE

TREE = kdtree.build_tree(LISTE_ETOILES_IMAGE, dim=2, dir=0)
for etoile in LISTE_ETOILES_REF:
    geometry.calcul_map_dtf_tfl_2d(etoile, LISTE_ETOILES_IMAGE, TREE)

bestr_score, bestmatchlist = -1, []
for etoile in LISTE_ETOILES_REF:

    D01 = identification.closest_dtf(etoile.F1, DATA_BASE)
    D01.load_gnomic()
    r_score, matchlist = identification.match_maps(etoile.normalized_map, D01.gnomic_projection_map, config.ID_THRESHOLD)
    if r_score > bestr_score:
        bestr_score, bestmatchlist, bestcentral_star, bestmatch = r_score, matchlist, D01, etoile 

    D02 = identification.closest_dtf(etoile.F2, DATA_BASE)
    D02.load_gnomic()
    r_score, matchlist = identification.match_maps(etoile.normalized_map, D02.gnomic_projection_map, config.ID_THRESHOLD)
    if r_score > bestr_score:
        bestr_score, bestmatchlist, bestcentral_star, bestmatch = r_score, matchlist, D02, etoile 

#chargement des etoiles depuis la bdd
bestmatchlist = [(data.get_by_attribute(DATA_BASE, "id", starid), etoile) for (starid, etoile) in bestmatchlist]
print(f"R_score avant iterations : {bestr_score}")

for i in range(config.N_ITE):
    (starA, etoileA), (starB, etoileB) = (random.choice(bestmatchlist), random.choice(bestmatchlist))
    if starA!=starB:
        imgmap = geometry.changement_image_vers_normalise(etoileA, etoileB, LISTE_ETOILES_IMAGE)
        mapbdd = geometry.calcul_gnomic(starA, starB, DATA_BASE)
        newr_score, newmatchlist = identification.match_maps(imgmap, mapbdd, config.ID_THRESHOLD2)
        if newr_score > bestr_score:
            bestr_score = newr_score
            bestmatchlist = newmatchlist 
            bestcentral_star = starA
            bestmatch = etoileA

print(f"R_score apres iterations: {bestr_score}")

for (star, etoile) in bestmatchlist:
    star.imagematch = etoile
    etoile.starmatch = star

results = display.affiche_resultat_pillow(LISTE_ETOILES_IMAGE, img_original, img_size, CON_DIC)
results.show() #Correspondances trouvées

#Sauvegarde
if config.SAVE_CENTROIDS:
    img_withcentroids = display.display_centroids(LISTE_ETOILES_IMAGE, img_original, img_size)
    img_withcentroids.show("Centroïdes des étoiles repérées")
    img_withcentroids.save(config.CENTROIDS_SAVE_PATH)
if config.SAVE_IMAGE:
    results.save(config.RESULTS_SAVE_PATH)