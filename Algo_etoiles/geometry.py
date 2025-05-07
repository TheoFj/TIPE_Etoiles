from math import*
import types_perso
import config
import kdtree
import heapq

#certaines fonctions ne sont pas utilisees

#fonctions sur vecteurs de dimension 2
def dot_prod2(u, v):
    return u[0]*v[0] + u[1]*v[1]
def norm2(v):
    return sqrt(v[0]**2 + v[1]**2)
def norm2_sqr(v):
    return v[0]**2 + v[1]**2
def round_xy(coord):
    return (round(coord[0]), round(coord[1]))
def vectpp(p1, p2):
    return (p2[0] - p1[0], p2[1] - p1[1])

#fonctions sur vecteurs de dimension 3
def dot_prod3(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def cross_prod3(u,v):
    return (u[1]*v[2] - u[2]*v[1], 
            u[2]*v[0] - u[0]*v[2], 
            u[0]*v[1] - u[1]*v[0])
def norm3(v):
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)
def norm3_sqr(v):
    return v[0]**2 + v[1]**2 + v[2]**2

#angles
def ang_dist(Star1, Star2): #https://en.wikipedia.org/wiki/Angular_distance
    p = dot_prod3(Star1.pos, Star2.pos)
    if p>1: #gère les imprecisions de calcul genre 1.0000000000000002
        return 0.0
    elif p<-1:
        return pi
    return acos(p)

def angle_2d(v1,v2, norm_prod):
    costheta = dot_prod2(v1,v2)/norm_prod
    if costheta>1: #gere les imprecisions de calcul genre 1.0000000000000002
        return 0.0
    elif costheta<-1:
        return pi
    return acos(costheta)


def calcul_double_triangle_feature(M0, M1, M2, M3):
    '''
    Calcule les donnees longueurs et angles à comparer (cf. section 4.3 https://www.mdpi.com/1424-8220/20/11/3027#sec4dot3-sensors-20-03027 )
    '''
    v01 = vectpp(M0, M1); v10 = (-v01[0],-v01[1])
    v02 = vectpp(M0, M2); v20 = (-v02[0],-v02[1])
    v03 = vectpp(M0, M3); v30 = (-v03[0],-v03[1])
    v21 = vectpp(M2, M1); v12 = (-v21[0],-v21[1])
    v23 = vectpp(M2, M3); v32 = (-v23[0],-v23[1])
    n01 = norm2(v01)
    n02 = norm2(v02)
    n03 = norm2(v03)
    n21 = norm2(v21)
    n32 = norm2(v32)
    return [angle_2d(v21, v20, n21*n02), angle_2d(v10, v12, n21*n01), 
            angle_2d(v01, v02, n01*n02), angle_2d(v30, v32, n03*n32), 
            angle_2d(v23, v20, n32*n02), angle_2d(v03, v02, n03*n02), 
            n01, n02, n21, n02, n03, n32]

def calcul_total_feature_length(dtf):
    return sum([dtf[i] for i in range(6,12)])

def calcul_gnomic(starA, starB, star_list):

    #PROJECTION SUR LE PLAN TANGEANT AU CERCLE UNITE
    #permet de simuler une photo prise par un capteur plan
    # (cf. https://www.mdpi.com/1424-8220/20/11/3027 section 4.1)
    C = []
    for star in star_list:
        p = dot_prod3(starA.pos, star.pos)
        if p>0:
            proj_vect = (star.x/p - starA.x), (star.y/p - starA.y), (star.z/p - starA.z)
            C.append((star, proj_vect))

    #CHANGEMENT DE BASE
    pAB = dot_prod3(starA.pos, starB.pos)
    E1 = (starB.x/pAB - starA.x), (starB.y/pAB - starA.y), (starB.z/pAB - starA.z)
    E2 = cross_prod3(starA.pos, E1) #rotation de pi/2 de E1
    k = norm3_sqr(E1)

    L = [(starandvect[0], (dot_prod3(starandvect[1], E1)/k, dot_prod3(starandvect[1], E2)/k)) for starandvect in C]

    return L


def changement_image_vers_normalise(M0, M1, image_star_list):
    #Normalisation de la liste des coordonnees des etoiles en fonction de deux points/etoiles de l'image
    E1 = vectpp(M0.pos, M1.pos)
    E2 = (-E1[1], E1[0])
    k = norm2_sqr(E1)
    return [(etoile, (dot_prod2(vectpp(M0.pos, etoile.pos), E1)/k, dot_prod2(vectpp(M0.pos, etoile.pos), E2)/k)) for etoile in image_star_list]

def changement_normalise_vers_image(M0, M1, map):
    #Calcul des coordonnées sur l'image d'une liste de coordonnées normalsiées
    E1 = vectpp(M0.pos, M1.pos)
    E2 = (-E1[1], E1[0])

    L=[]
    for starandvect in map:
        x = E1[0]*starandvect[1][0] + E2[0]*starandvect[1][1] + M0.x
        y = E1[1]*starandvect[1][0] + E2[1]*starandvect[1][1] + M0.y
        L.append((x,y))
    return L


def calcul_map_dtf_tfl_2d(M0, liste_etoiles_image, tree):
    '''
    Calculs à faire sur chaque étoile (ou une partie des étoiles) de l'image
    '''
    nearest_4 = kdtree.nearest_nstars(tree, M0.pos, n=4, dim=2, dir=0, withtarget=False)
    M4 = heapq.heappop(nearest_4)[2]
    M3 = heapq.heappop(nearest_4)[2]
    M2 = heapq.heappop(nearest_4)[2]
    M1 = heapq.heappop(nearest_4)[2]
    M0.closest_star = M1

    E1 = vectpp(M0.pos, M1.pos)
    E2 = (-E1[1], E1[0])
    k = norm2_sqr(E1)

    #PROJECTIONS DANS LE NOUVEAU REPERE
    L = [(M, (dot_prod2(vectpp(M0.pos, M.pos), E1)/k, dot_prod2(vectpp(M0.pos, M.pos), E2)/k)) for M in liste_etoiles_image]

    N0 = (0., 0.)       
    N1 = (1., 0.)
    N2 = (dot_prod2(vectpp(M0.pos, M2.pos), E1)/k, dot_prod2(vectpp(M0.pos, M2.pos), E2)/k)
    N3 = (dot_prod2(vectpp(M0.pos, M3.pos), E1)/k, dot_prod2(vectpp(M0.pos, M3.pos), E2)/k)
    N4 = (dot_prod2(vectpp(M0.pos, M4.pos), E1)/k, dot_prod2(vectpp(M0.pos, M4.pos), E2)/k)

    M0.F1 = F1 = calcul_double_triangle_feature(N0, N1, N2, N3)
    M0.F2 = F2 = calcul_double_triangle_feature(N0, N2, N3, N4)
    M0.F1_length = calcul_total_feature_length(F1)
    M0.F2_length = calcul_total_feature_length(F2)
    M0.normalized_map = L
    return
