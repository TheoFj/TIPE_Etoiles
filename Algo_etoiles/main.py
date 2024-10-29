from PIL import Image
from math import *
import config
import types_perso
import display
import data
import geometry
import identification
import image_treatment

#ALGORITHME: https://www.mdpi.com/1424-8220/20/11/3027


def get_by_attribute(data_base, attribute, val):
    if attribute == "bayer":
        return next((star for star in data_base if star.bayer == val), [None])
    elif attribute == "id":
        return next((star for star in data_base if star.id == val), [None])
    elif attribute == "hip":
        return next((star for star in data_base if star.hip == val), [None])


#LISTE_ETOILES_REF = data.choose_random(LISTE_ETOILES_IMAGE, LAMBD)
LISTE_ETOILES_REF = image_treatment.LISTE_ETOILES_IMAGE

for etoile in LISTE_ETOILES_REF:
    geometry.calcul_map_dtf_tfl_2d(etoile, image_treatment.LISTE_ETOILES_IMAGE)


bestr_score, bestmatchlist = -1, []
for etoile in LISTE_ETOILES_REF:

    D01 = identification.closest_dtf(etoile.F1, data.DATA_BASE)
    r_score, matchlist, match = identification.identify(D01, etoile)
    if r_score > bestr_score:
        bestr_score, bestmatchlist, bestcentral_star, bestmatch = r_score, matchlist, D01, match 

    D02 = identification.closest_dtf(etoile.F2, data.DATA_BASE)
    r_score, matchlist, match = identification.identify(D02, etoile)
    if r_score > bestr_score:
        bestr_score, bestmatchlist, bestcentral_star, bestmatch = r_score, matchlist, D02, match 

for (starid, etoile) in bestmatchlist:
    star = get_by_attribute(data.DATA_BASE, "id", starid)
    star.imagematch = etoile
    etoile.starmatch = star

back_projection = display.affiche_etoiles(geometry.changement_normalise_vers_image(bestmatch, bestmatch.closest_star, bestcentral_star.gnomic_projection_map), bestmatch, image_treatment.SIZE, image_treatment.IMAGE_ORIGINAL)
back_projection.show()

results = display.affiche_resultat_pillow(image_treatment.LISTE_ETOILES_IMAGE, image_treatment.IMAGE_ORIGINAL, image_treatment.SIZE, data.CON_DIC)
results.show("Étoiles reconnues")