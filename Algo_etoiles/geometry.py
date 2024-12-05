from math import*
import types_perso
import config

#certaines fonctions ne sont pas utilisees

#fonctions sur vecteurs de dimension 2
def flat_dist2(pos1, pos2):
    return sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
def flat_dist2_sqr(pos1, pos2): #permet d'eviter de calculer une racine
    return (pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2
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
def flat_dist3(pos1, pos2):
    return sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2 + (pos1[2] - pos2[2])**2)
def flat_dist3_sqr(pos1, pos2): #permet d'eviter de calculer une racine
    return (pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1] )**2 + (pos1[2] - pos2[2])**2
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
    p = dot_prod3(Star1.xyz, Star2.xyz)
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
    return [angle_2d(v21, v20, n21*n02),
            angle_2d(v10, v12, n21*n01), 
            angle_2d(v01, v02, n01*n02), 
            angle_2d(v30, v32, n03*n32), 
            angle_2d(v23, v20, n32*n02), 
            angle_2d(v03, v02, n03*n02), 
            n01, 
            n02, 
            n21, 
            n02,
            n03, 
            n32, 
            ]

def calcul_total_feature_length(dtf):
    return sum([dtf[i] for i in range(6,12)])

def closest_nstars(starandvect_list, n, with_zero, dim):
    '''
    Renvoie les n vecteurs les plus proches de 0, en dimension 2 ou 3, en incluant le 0 ou non.
    O(len(starandvect_list)*n) ATTENTION: compare par rapport à (0,0)/(0,0,0), pour comparer à un autre point, mettre en entree les vecteurs "P0->P1"
    '''
    if dim==2: L_dist = [(norm2_sqr(starandvect[1]), starandvect) for starandvect in starandvect_list]
    elif dim==3: L_dist = [(norm3_sqr(starandvect[1]), starandvect) for starandvect in starandvect_list]
    else: print("????????????")
    best_n = []

    for (dist, starandvect) in L_dist:
        i=0
        while i<(len(best_n)):
            if dist <= best_n[i][0]:
                break
            i+=1
        best_n.insert(i, (dist, starandvect))
        best_n = best_n[:(n+1)]

    if with_zero:
        return [best_n[i][1] for i in range(0,n)] #on ne renvoie pas le dernier, indicé n, car on ne renvoie que les n meilleures
    else:
        return [best_n[i][1] for i in range(1,n+1)] #on ne renvoie pas le numero 0 car il s'agit de l'etoile "pos"

def furthest_stars(matchlist):

    a = matchlist[0]
    b = matchlist[0]
    maxdsqr = 0
    n = len(matchlist)
    for i in range(n):
        for j in range(i+1,n):
            c = matchlist[i][1]
            d = matchlist[j][1]
            dsqr = (d.x-c.x)**2 + (d.y-c.y)**2
            if dsqr>maxdsqr:
                a = matchlist[i]
                b = matchlist[j]
                maxdsqr = dsqr
    return (a, b)

def calcul_gnomic(starA, starB, star_list):

    #PROJECTION SUR LE PLAN TANGEANT AU CERCLE UNITE
    #permet de simuler une photo prise par un capteur plan
    # (cf. https://www.mdpi.com/1424-8220/20/11/3027 section 4.1)
    C = []
    for star in star_list:
        p = dot_prod3(starA.xyz, star.xyz)
        if p>0:
            proj_vect = (star.x/p - starA.x), (star.y/p - starA.y), (star.z/p - starA.z)
            C.append((star, proj_vect))

    #CHANGEMENT DE BASE
    pAB = dot_prod3(starA.xyz, starB.xyz)
    E1 = (starB.x/pAB - starA.x), (starB.y/pAB - starA.y), (starB.z/pAB - starA.z)
    E2 = cross_prod3(starA.xyz, E1) #rotation de pi/2 de E1
    k = norm3_sqr(E1)

    L = [(starandvect[0], (dot_prod3(starandvect[1], E1)/k, dot_prod3(starandvect[1], E2)/k)) for starandvect in C]

    return L


def changement_image_vers_normalise(M0, M1, image_star_list):
    #Normalisation de la liste des coordonnees des etoiles en fonction de deux points/etoiles de l'image
    E1 = vectpp(M0.xy, M1.xy)
    E2 = (-E1[1], E1[0])
    k = norm2_sqr(E1)
    return [(etoile, (dot_prod2(vectpp(M0.xy, etoile.xy), E1)/k, dot_prod2(vectpp(M0.xy, etoile.xy), E2)/k)) for etoile in image_star_list]

def changement_normalise_vers_image(M0, M1, map):
    #Calcul des coordonnées sur l'image d'une liste de coordonnées normalsiées
    E1 = vectpp(M0.xy, M1.xy)
    E2 = (-E1[1], E1[0])

    L=[]
    for starandvect in map:
        x = E1[0]*starandvect[1][0] + E2[0]*starandvect[1][1] + M0.x
        y = E1[1]*starandvect[1][0] + E2[1]*starandvect[1][1] + M0.y
        L.append((x,y))
    return L


def calcul_map_dtf_tfl_2d(M0, image_star_list):
    '''
    Calculs à faire sur chaque étoile (ou une partie des étoiles) de l'image
    '''

    starandvect_list = [(M, vectpp(M0.xy, M.xy)) for M in image_star_list]
    best_4 = closest_nstars(starandvect_list, 4, False, 2)
    M0.closest_star = best_4[0][0]
    E1 = best_4[0][1]
    E2 = (-E1[1], E1[0])
    k = norm2_sqr(E1)

    #PROJECTIONS DANS LA NOUVELLE BASE
    L = [(starandvect[0], (dot_prod2(starandvect[1], E1)/k, dot_prod2(starandvect[1], E2)/k)) for starandvect in starandvect_list]

    M0xy = (0., 0.)       
    M1xy = (1., 0.)
    M2xy = (dot_prod2(best_4[1][1], E1)/k, dot_prod2(best_4[1][1], E2)/k)
    M3xy = (dot_prod2(best_4[2][1], E1)/k, dot_prod2(best_4[2][1], E2)/k)
    M4xy = (dot_prod2(best_4[3][1], E1)/k, dot_prod2(best_4[3][1], E2)/k)

    M0.F1 = F1 = calcul_double_triangle_feature(M0xy, M1xy, M2xy, M3xy)
    M0.F2 = F2 = calcul_double_triangle_feature(M0xy, M2xy, M3xy, M4xy)
    M0.F1_length = calcul_total_feature_length(F1)
    M0.F2_length = calcul_total_feature_length(F2)
    M0.normalized_map = L
    return
