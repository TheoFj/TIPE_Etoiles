from math import*
import kdtree
import heapq
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
def angle_2d(v1, v2, norm_prod):
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
    n01 = norm2(v01); n02 = norm2(v02); n03 = norm2(v03); n21 = norm2(v21); n32 = norm2(v32)
    return [angle_2d(v21, v20, n21*n02), angle_2d(v10, v12, n21*n01), 
            angle_2d(v01, v02, n01*n02), angle_2d(v30, v32, n03*n32), 
            angle_2d(v23, v20, n32*n02), angle_2d(v03, v02, n03*n02), 
            n01, n02, n21, n02, n03, n32]

def calcul_total_feature_length(dtf):
    return sum([dtf[i] for i in range(6,12)])

def proj(star, D0):
    #PROJECTION SUR LE PLAN TANGEANT AU CERCLE UNITE
    #permet de simuler une photo prise par un capteur plan
    p = dot_prod3(D0.pos, star.pos)
    return (star.x/p - D0.x), (star.y/p - D0.y), (star.z/p - D0.z)

def calcul_gnomic_dtf_tfl(D0, star_list, tree, L2):
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
    E2 = cross_prod3(D0.pos, E1) #rotation de pi/2 de E1
    k = norm3_sqr(E1)

    L = []
    for star in star_list:
        p = dot_prod3(D0.pos, star.pos)
        proj_vect = proj(star, D0)
        if p>0 and norm3_sqr(proj_vect) < L2: #condition pour ne pas prendre des etoiles trop eloignees ou bien derriere la Terre
            L.append((star, (dot_prod3(proj_vect, E1)/k, dot_prod3(proj_vect, E2)/k)))
    
    M0 = (0., 0.)
    M1 = (1., 0.)
    p2 = proj(D2, D0)
    p3 = proj(D3, D0)
    M2 = (dot_prod3(p2, E1)/k, dot_prod3(p2, E2)/k)
    M3 = (dot_prod3(p3, E1)/k, dot_prod3(p3, E2)/k)


    dtf = calcul_double_triangle_feature(M0, M1, M2, M3)
    
    #Enregistrement dans l'objet Star D0
    D0.gnomic_projection_map = L
    D0.double_triangle_feature = dtf
    D0.total_feature_length = calcul_total_feature_length(dtf)
    return

def calcul_gnomic(starA, starB, star_list):

    #PROJECTION SUR LE PLAN TANGEANT AU CERCLE UNITE
    #permet de simuler une photo prise par un capteur plan
    # (cf. https://www.mdpi.com/1424-8220/20/11/3027 section 4.1)
    C = []
    for star in star_list:
        p = dot_prod3(starA.pos, star.pos)
        if p>0:
            proj_vect = proj(star, starA)
            C.append((star, proj_vect))

    #CHANGEMENT DE BASE
    E1 = proj(starB,starA)
    E2 = cross_prod3(starA.pos, E1) #rotation de pi/2 de E1
    k = norm3_sqr(E1)

    return [(starandvect[0], (dot_prod3(starandvect[1], E1)/k, dot_prod3(starandvect[1], E2)/k)) for starandvect in C]

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

def changement_image_vers_normalise(M0, M1, image_star_list):
    #Normalisation de la liste des coordonnees des etoiles en fonction de deux points/etoiles de l'image
    E1 = vectpp(M0.pos, M1.pos)
    E2 = (-E1[1], E1[0])
    k = norm2_sqr(E1)
    return [(M, (dot_prod2(vectpp(M0.pos, M.pos), E1)/k, dot_prod2(vectpp(M0.pos, M.pos), E2)/k)) for M in image_star_list]

def changement_normalise_vers_image(M0, M1, map):
    #Calcul des coordonnées sur l'image d'une liste de coordonnées normalsiées
    E1 = vectpp(M0.pos, M1.pos)
    E2 = (-E1[1], E1[0])

    return [(E1[0]*starandvect[1][0] + E2[0]*starandvect[1][1] + M0.x, E1[1]*starandvect[1][0] + E2[1]*starandvect[1][1] + M0.y) for starandvect in map]
