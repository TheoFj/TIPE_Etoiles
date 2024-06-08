from math import *
from PIL import Image, ImageDraw, ImageFont
import csv
import random
#import matplotlib.pyplot as plt

#LIRE https://www.mdpi.com/1424-8220/20/11/3027

"""
CONSTANTES
"""
FOV = 90
L = FOV*pi/180 #fov en radians 
L2 = (L**2)/2 #rayon du disque des etoiles prises en compte dans "calcul_gnomic_dtf_tfl" (cf. 4.1 https://www.mdpi.com/1424-8220/20/11/3027)
DATA_BASE_NAME = "athyg_modified_vmagmax6.csv"
DATA_BASE_PATH = "Algo_etoiles/databasecsv/"+DATA_BASE_NAME
DATA_BASE_TREATED_PATH = "Algo_etoiles/databasecsv/treated_"+DATA_BASE_NAME
DISPLAY_MODE = "mixt" #bayerflam / hip / mixt

"""
TYPES
"""
class Star: #Type enregistrement pour les étoiles de la base de données

    def __init__(self, id, hip, bayer, flam, con, proper, ra, dec, mag, full, gen, unicode): #ra et dec en heures et degres dans le csv
        self.id = id            #Identifiant dans la base de donnéee
        self.hip = hip          #Identifiant hipparcos (si il existe)
        self.bayer = bayer      #Désignation de bayer  (si elle existe)
        self.flam = flam        #Désignation de flamsteed (si elle existe)
        self.con = con          #Abréviation de la constellation du système international
        self.proper = proper    #Nom commun (si il existe)
        self.ra = ra*pi/12      #Ascention droite avec conversion heures(15deg)->radians
        self.dec = dec*pi/180   #Déclinaison avec conversion degrés->radians
        self.mag = mag          #Magnitude (correspond en gros à du -log(luminosité))
        self.full_letter = full #Désignation de bayer en toutes lettres
        self.genitive = gen     #Génitif latin de la constellation en toutes lettres
        self.uni = unicode      #Caractère unicode grec de la désignation de bayer

        self.simbad = ""
        self.wiki = ""
        if self.hip != None: 
            self.simbad = "http://cdsportal.u-strasbg.fr/?target=hip"+self.hip #sinon utiliser les coordonnees mais flm là
        if self.full_letter != None: #ameliorer pour adapter aux Sig-1, Sig-2 etc.
            self.wiki = "https://en.wikipedia.org/wiki/"+self.full_letter+"_"+self.genitive


        cdec, sdec, cra, sra = cos(self.dec), sin(self.dec), cos(self.ra), sin(self.ra)
        self.xyz = (self.x, self.y, self.z) = (cdec*cra, cdec*sra, sdec) #Position sur la sphere celeste unité

        self.double_triangle_feature = None     #Longueurs et angles des deux triangles formé par les 4 etoiles les plus proches
        self.total_feature_length = None        #Longueur totale des segments de double_triangle_feature
        self.gnomic_projection_map = None       #Etoiles proches et leurs coordonnees dans la base adaptée à cet étoile
    
    def save_gnomic(self):
        path = "Algo_etoiles/treated_data/"+str(self.id)+".csv"
        with open(path, 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            #spamwriter.writerow([self.x, self.y, self.z]) #LIGNE 1 : X, Y, Z
            #spamwriter.writerow(self.double_triangle_feature) #LIGNE 2 : DTF
            #spamwriter.writerow([self.total_feature_length]) #LIGNE 3 : TFL
            for (star, (x,y)) in self.gnomic_projection_map: #LIGNE 4 et + : projection
                spamwriter.writerow([star.id, x, y])
            csvfile.close()
    

"""
FONCTIONS PRATIQUES
""" 
#certaines fonctions ne sont pas utilisees

def intbis(s):
    return None if (s=='' or s=='#N/A') else int(s)
def strbis(s):
    return None if (s=='' or s=='#N/A') else str(s)
def fltbis(s):
    return None if (s=='' or s=='#N/A') else float(s)

def ang_dist(Star1, Star2): #https://en.wikipedia.org/wiki/Angular_distance
    p = dot_prod3(Star1.xyz, Star2.xyz)
    if p>1: #gere les imprecisions de calcul genre 1.0000000000000002
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

'''
FONCTIONS IMPORTANTES
'''

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
    
def calcul_double_triangle_feature(M0, M1, M2, M3): #appelé sur M0,M2,M3,M4 pour construire F2
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
    return [angle_2d(v21, v20, n21*n02), #attention: liste ici et tuple dans le main
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

def calcul_gnomic_dtf_tfl(D0, star_list):
    '''
    Calculs à faire sur chaque étoile de la base donnée, 
    à faire idéalement dans un code séparé pour ne pas le recalculer à chaque fois
    car ces calculs sont indépendants de l'imaage.
    '''
    
    #PROJECTION SUR LE PLAN TANGEANT AU CERCLE UNITE
    #permet de simuler une photo prise par un capteur plan
    # (cf. https://www.mdpi.com/1424-8220/20/11/3027 section 4.1)
    # (cf. __ a completer __)
    C = []
    for star in star_list:
        p = dot_prod3(D0.xyz, star.xyz)
        proj_vect = (star.x/p - D0.x), (star.y/p - D0.y), (star.z/p - D0.z)
        if p>0 and norm3_sqr(proj_vect) < L2:
            C.append((star, proj_vect))

    #CHANGEMENT DE BASE EN PRENANT LETOILE LA PLUS PROCHE
    # (cf. https://www.mdpi.com/1424-8220/20/11/3027 section 4.2)
    # permet d'aligner les systemes de coordonnees basee de donnee/image afin d'avoir des donnees comparables
    best_3 = closest_nstars(C, 3, False, 3)
    E1 = best_3[0][1]
    E2 = cross_prod3(D0.xyz, E1) #rotation de pi/2 de E1
    k = norm3_sqr(best_3[0][1])

    L = [(starandvect[0], (dot_prod3(starandvect[1], E1)/k, dot_prod3(starandvect[1], E2)/k)) for starandvect in C]

    C0 = (D0, (0., 0.))
    C1 = (best_3[0][0], (1., 0.))
    C2 = (best_3[1][0], (dot_prod3(best_3[1][1], E1)/k, dot_prod3(best_3[1][1], E2)/k))
    C3 = (best_3[2][0], (dot_prod3(best_3[2][1], E1)/k, dot_prod3(best_3[2][1], E2)/k))

    DTF = calcul_double_triangle_feature(C0[1], C1[1], C2[1], C3[1])
    
    #Enregistrement dans l'objet Star D0
    D0.gnomic_projection_map = L
    D0.double_triangle_feature = DTF
    D0.total_feature_length = calcul_total_feature_length(DTF)
    return


'''
IMPORTATION ET CALCULS SUR LA BASE DE DONNEES
'''

csv_file = open(DATA_BASE_PATH)
next(csv_file)

reader = csv.reader(csv_file, delimiter=',')
DATA_BASE = [Star(                      #construction de DATA_BASE: liste d'objets de type Star, conversions un peu inutiles mais copiees
                int(row[0]),            #id
                strbis(row[1]),         #hip (ici en string mais en int dans main)
                strbis(row[2]),         #bayer
                strbis(row[3]),         #flam
                str(row[4]),            #const
                strbis(row[5]),         #proper
                float(row[6]),          #ra
                float(row[7]),          #dec
                float(row[8]),          #mag
                strbis(row[9]),         #full
                strbis(row[10]),        #gen
                strbis(row[11]),        #greek_unicode
                ) for row in reader]

csv_file.close()

csv_file2 = open(DATA_BASE_PATH)
csv_save = open(DATA_BASE_TREATED_PATH, 'w', newline='')
reader2 = csv.reader(csv_file2, delimiter=',')
writer = csv.writer(csv_save, delimiter=',')

row1 = next(reader2)
writer.writerow(row1+["x","y","z","tfl","th01","th021","th12","th023","th03","th23","l01","l02","l12","l02","l03","l23","wikipedia","simbad"])

for star in DATA_BASE:
    row = next(reader2)
    assert(star.id == int(row[0]))
    calcul_gnomic_dtf_tfl(star,  DATA_BASE)
    writer.writerow(row+[star.x,star.y,star.z,star.total_feature_length]+star.double_triangle_feature+[star.wiki, star.simbad])
    star.save_gnomic()


csv_save.close()
csv_file2.close()