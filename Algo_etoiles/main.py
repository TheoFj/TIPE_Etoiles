from PIL import Image
from math import *
import csv
import config
import types_perso
import display
import data
import geometry
import identification
import img_fnc as img_fnc

#ALGORITHME: https://www.mdpi.com/1424-8220/20/11/3027

#Traitement image:

img_original = Image.open(config.IMAGE_PATH)
img_grayscale = img_original.convert('L')

img_size = img_width, img_height = img_original.size

if config.CONVOLUTION_OR_NOT == True:
    img_treated = img_fnc.convolve(img_grayscale, img_width, img_height, img_fnc.MAT)
    img_treated.show()
else:
    img_treated = img_grayscale

img_nb = img_treated.point(img_fnc.pixel_to_NB, mode='1')

img_nb.show()

image_temp = img_nb.copy() #utile car la fonction "new_star" a besoin de modifier l'image
LISTE_ETOILES_IMAGE = []
for y in range(img_height): #remplit la liste des coordonnéees des étoiles sur l'image
    for x in range(img_width):
        if image_temp.getpixel((x,y)) == 1:
            img_fnc.new_star(image_temp, (x,y), LISTE_ETOILES_IMAGE)

img_withcentroids = display.display_centroids(LISTE_ETOILES_IMAGE, img_original, img_size)
#img_withcentroids.show("Centroïdes des étoiles repérées")

#Import BDD

DATA_BASE = data.parse_database_file(config.DATA_BASE_TREATED_PATH)
CON_DIC = data.parse_constellation_file(config.CON_DIC_PATH)

#Accesseurs

def get_by_attribute(data_base, attribute, val):
    if attribute == "bayer":
        return next((star for star in data_base if star.bayer == val), [None])
    elif attribute == "id":
        return next((star for star in data_base if star.id == val), [None])
    elif attribute == "hip":
        return next((star for star in data_base if star.hip == val), [None])

#Boucles principales

#LISTE_ETOILES_REF = data.choose_random(LISTE_ETOILES_IMAGE, LAMBD)
LISTE_ETOILES_REF = LISTE_ETOILES_IMAGE

for etoile in LISTE_ETOILES_REF:
    geometry.calcul_map_dtf_tfl_2d(etoile, LISTE_ETOILES_IMAGE)

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

(idA, etoileA), (idB, etoileB) = geometry.furthest_stars(bestmatchlist)
starA = get_by_attribute(DATA_BASE, "id", idA)
starB = get_by_attribute(DATA_BASE, "id", idB)
imgmap = geometry.changement_image_vers_normalise(etoileA, etoileB, LISTE_ETOILES_IMAGE)
mapbdd = geometry.calcul_gnomic(starA, starB, DATA_BASE)
newr_score, newmatchlist = identification.match_maps(imgmap, mapbdd, config.ID_THRESHOLD2)

print(f"r_score étoiles eloignées: {newr_score}\n Ancien r_score: {bestr_score}")

#Affichage des résultats:

if newr_score >= bestr_score:
    for (star, etoile) in newmatchlist:
        star.imagematch = etoile
        etoile.starmatch = star
    back_projection = display.affiche_etoiles(geometry.changement_normalise_vers_image(etoileA, etoileB, mapbdd), etoileA, img_size, img_original)

else:
    for (starid, etoile) in bestmatchlist:
        star = get_by_attribute(DATA_BASE, "id", starid)
        star.imagematch = etoile
        etoile.starmatch = star
    back_projection = display.affiche_etoiles(geometry.changement_normalise_vers_image(bestmatch, bestmatch.closest_star, bestcentral_star.gnomic_projection_map), bestmatch, img_size, img_original)


results = display.affiche_resultat_pillow(LISTE_ETOILES_IMAGE, img_original, img_size, CON_DIC)


back_projection.show() #Projection de la base de données sur l'image: permet de
results.show() #Correspondances trouvées


#Sauvegarde

if config.SAVE_CENTROIDS:
    img_withcentroids.save(config.CENTROIDS_SAVE_PATH)
if config.SAVE_IMAGE:
    results.save(config.RESULTS_SAVE_PATH)
if config.SAVE_PDF:
    pass