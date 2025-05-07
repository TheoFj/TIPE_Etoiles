from math import *
import csv
import geometry
import types_perso
import kdtree
import heapq

#LIRE https://www.mdpi.com/1424-8220/20/11/3027

"""
CONSTANTES
"""
FOV = 50
L = FOV*pi/180 #fov en radians 
L2 = (L**2)/2 #rayon du disque des etoiles prises en compte dans "calcul_gnomic_dtf_tfl" (cf. 4.1 https://www.mdpi.com/1424-8220/20/11/3027) (ce calcul est super chelou)
DATA_BASE_NAME = "athyg_modified_vmagmax6.csv"
DATA_BASE_PATH = "Algo_etoiles/databasecsv/"+DATA_BASE_NAME
DATA_BASE_TREATED_PATH = "Algo_etoiles/databasecsv/treated_"+DATA_BASE_NAME

def intbis(s):
    return None if (s=='' or s=='#N/A') else int(s)
def strbis(s):
    return None if (s=='' or s=='#N/A') else str(s)
def fltbis(s):
    return None if (s=='' or s=='#N/A') else float(s)

'''
FONCTIONS IMPORTANTES
'''

def proj(star, D0):
    #PROJECTION SUR LE PLAN TANGEANT AU CERCLE UNITE
    #permet de simuler une photo prise par un capteur plan
    p = geometry.dot_prod3(D0.pos, star.pos)
    return (star.x/p - D0.x), (star.y/p - D0.y), (star.z/p - D0.z)

def calcul_gnomic_dtf_tfl(D0, star_list, tree):
    '''
    Calculs à faire sur chaque étoile de la base donnée, 
    à exécuter avant pour ne pas refaire ces calculs (indépendants de l'image) à chaque image.
    '''
    # (cf. https://www.mdpi.com/1424-8220/20/11/3027 section 4.2)

    nearest_3 = kdtree.nearest_nstars(tree, D0.pos, n=3, dim=3, dir=0, withtarget=False)
    D3 = heapq.heappop(nearest_3)[2]
    D2 = heapq.heappop(nearest_3)[2]
    D1 = heapq.heappop(nearest_3)[2]

    # (E1, E2) est la nouvelle base dans laquelle sont exprimees les coordonnes des etoiles projetees
    # permet d'aligner les systemes de coordonnees basee de donnee/image afin d'avoir des donnees comparables
    E1 = proj(D1, D0)
    E2 = geometry.cross_prod3(D0.pos, E1) #rotation de pi/2 de E1
    k = geometry.norm3_sqr(E1)

    L = []
    for star in star_list:
        p = geometry.dot_prod3(D0.pos, star.pos)
        proj_vect = proj(star, D0)
        if p>0 and geometry.norm3_sqr(proj_vect) < L2: #condition pour ne pas prendre des etoiles trop eloignees ou bien derriere la Terre
            L.append((star, (geometry.dot_prod3(proj_vect, E1)/k, geometry.dot_prod3(proj_vect, E2)/k)))
    
    M0 = (0., 0.)
    M1 = (1., 0.)
    p2 = proj(D2, D0)
    p3 = proj(D3, D0)
    M2 = (geometry.dot_prod3(p2, E1)/k, geometry.dot_prod3(p2, E2)/k)
    M3 = (geometry.dot_prod3(p3, E1)/k, geometry.dot_prod3(p3, E2)/k)


    dtf = geometry.calcul_double_triangle_feature(M0, M1, M2, M3)
    
    #Enregistrement dans l'objet Star D0
    D0.gnomic_projection_map = L
    D0.double_triangle_feature = dtf
    D0.total_feature_length = geometry.calcul_total_feature_length(dtf)
    return

'''
LECTURE, CALCULS ET ECRITURE SUR LA BDD
'''

csv_file = open(DATA_BASE_PATH, "r", encoding="utf-8")
next(csv_file)

reader = csv.reader(csv_file, delimiter=',')
DATA_BASE = [types_perso.Star(                      #construction de DATA_BASE: liste d'objets de type Star, conversions un peu inutiles mais copiees
                    id=int(row[0]),
                    hip=strbis(row[1]),
                    bayer=strbis(row[2]),
                    flam=strbis(row[3]),
                    con=str(row[4]),
                    proper=strbis(row[5]),
                    ra=float(row[6]),
                    dec=float(row[7]),
                    mag=float(row[8]),
                    full=strbis(row[9]),
                    gen=strbis(row[10]),
                    greek_bay=strbis(row[11]),
                    x=None, y=None, z=None, 
                    tfl=None, 
                    dtf=None, 
                    ) for row in reader]
csv_file.close()


for star in DATA_BASE:
    star.calcxyz()
TREE = kdtree.build_tree(DATA_BASE.copy(), dim=3, dir=0)

csv_file2 = open(DATA_BASE_PATH, "r", encoding="utf-8")
csv_save = open(DATA_BASE_TREATED_PATH, 'w', newline='', encoding="utf-8")
reader2 = csv.reader(csv_file2, delimiter=',')
writer = csv.writer(csv_save, delimiter=',')

row1 = next(reader2)
writer.writerow(row1+["x","y","z","tfl","th01","th021","th12","th023","th03","th23","l01","l02","l12","l02","l03","l23"])

for star in DATA_BASE:
    row = next(reader2)
    assert(star.id == int(row[0]))
    calcul_gnomic_dtf_tfl(star, DATA_BASE, TREE)
    writer.writerow(row+[star.x,star.y,star.z,star.total_feature_length]+star.double_triangle_feature)
    star.save_gnomic()

csv_save.close()
csv_file2.close()