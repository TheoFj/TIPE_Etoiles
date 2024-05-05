from math import *
from PIL import Image
import csv
#import matplotlib.pyplot as plt
import pygame
import sys
import random
import time

#LIRE https://www.mdpi.com/1424-8220/20/11/3027#sec4dot3-sensors-20-03027 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

"""
CONSTANTES
"""
FOV = 90
L = FOV*pi/180 #fov en radians
L2 = (L**2)/2
LAMBD = 20 #nombre maximum d'etoiles de "reference"
IMAGE_PATH = "Algo_etoiles/images/UMa2.png"
DATA_BASE_PATH = "Algo_etoiles/database/UMa_vmagmax6.csv"
TEMP_IMAGE_PATH = "temp.png"


"""
TYPES
"""
class Star:

    def __init__(self, id, hip, bayer, flam, con, proper, ra, dec, mag, wikipedia): #ra et dec en heures et degres dans le csv
        self.id = id
        self.hip = hip
        self.bayer = bayer
        self.flam = flam
        self.con = con
        self.proper = proper
        self.ra = ra*pi/12      #conversion heures(15deg) ->radians
        self.dec = dec*pi/180   #conversion degres ->radians
        self.mag = mag
        self.wiki = wikipedia

        cdec, sdec, cra, sra = cos(self.dec), sin(self.dec), cos(self.ra), sin(self.ra)
        self.xyz = (self.x, self.y, self.z) = (cdec*cra, cdec*sra, sdec) #position sur la sphere celeste unite
        #self.etheta = (-sra, cra, 0)
        #self.ephi = (-sdec*cra, -sdec*sra, cdec)

        self.imagematch = None
        self.double_triangle_feature = None
        self.total_feature_length = None
        self.gnomic_projection_map = None

        self.display_name = None
        if self.bayer != None:
            self.display_name = self.bayer + " " + self.con
        elif self.flam != None: 
            self.display_name = self.flam + " " + self.con

    def affiche(self, window):
        if self.imagematchxy != None:
            print(self.imagematchxy, self.display_name)
            pygame.draw.circle(window, (0,255,0), self.imagematchxy, 8, 1)
            window.blit(FONT.render(self.display_name, False, (0, 255, 0)), (self.imagematchxy[0]-25, self.imagematchxy[1]-25))
        
class Etoile_image:

    def __init__(self, pos):
        self.xy = self.x, self.y = pos

        self.F1 = None
        self.F1_length = None
        self.F2 = None
        self.F2_length = None
        self.normalized_map = None

        self.starmatch = None

    def affiche(self, window):
        if self.starmatch != None:
            print(self.xy, self.starmatch)
            pygame.draw.circle(window, (0,255,0), self.xy, 8, 1)
            window.blit(FONT.render(self.starmatch.display_name, False, (0, 255, 0)), (self.xy[0]-25, self.xy[1]-25))
        else:
            print(self.xy, None)
            pygame.draw.circle(window, (255,0,0), self.xy, 8, 1)



"""
FONCTIONS PRATIQUES
""" 

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


def flat_dist2(pos1, pos2):
    return sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
def flat_dist2_sqr(pos1, pos2):
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


def flat_dist3(pos1, pos2):
    return sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2 + (pos1[2] - pos2[2])**2)
def flat_dist3_sqr(pos1, pos2):
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

def pygameblitfrompillow(image_pillow, window):
    image_pillow.save(TEMP_IMAGE_PATH)
    image_pygame = pygame.image.load(TEMP_IMAGE_PATH)
    window.blit(image_pygame,(0,0))

'''
FONCTIONS IMPORTANTES
'''

def closest_nstars(starandvect_list, n, with_zero, dim): #O(len(starandvect_list)*n) ATTENTION: compare par rapport à (0,0)/(0,0,0), pour comparer à un autre point, mettre en entree les vecteurs "P0->P1"
    if dim==2: L_dist = [(norm2_sqr(starandvect[1]), starandvect) for starandvect in starandvect_list]
    elif dim==3: L_dist = [(norm3_sqr(starandvect[1]), starandvect) for starandvect in starandvect_list]
    else: print("????????????")
    best_n = []

    for couple in L_dist:
        i=0
        while i<(len(best_n)):
            if couple[0]<=best_n[i][0]:
                break
            i+=1
        best_n.insert(i, couple)
        best_n = best_n[:(n+1)]

    if with_zero:
        return [best_n[i][1] for i in range(0,n)] #on ne renvoie pas le dernier, indicé n, car on ne renvoie que les n meilleures
    else:
        return [best_n[i][1] for i in range(1,n+1)] #on ne renvoie pas le numero 0 car il s'agit de l'etoile "pos"
    
def calcul_double_triangle_feature(M0, M1, M2, M3): #appelé sur M0,M2,M3,M4 pour construire F2
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
    return (angle_2d(v21, v20, n21*n02),
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
            )

def calcul_total_feature_length(dtf):
    return sum([dtf[i] for i in range(6,12)])

def calcul_gnomic_dtf_tfl(D0, star_list):
    
    #PROJECTION SUR LE PLAN TANGEANT AU CERCLE UNITE
    C = []
    for star in star_list:
        p = dot_prod3(D0.xyz, star.xyz)
        proj_vect = (star.x/p - D0.x), (star.y/p - D0.y), (star.z/p - D0.z)
        if p>0 and norm3_sqr(proj_vect) < L2:
            C.append((star, proj_vect))

    #CHANGEMENT DE BASE EN PRENANT LETOILE LA PLUS PROCHE
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
    
    D0.gnomic_projection_map = L
    D0.double_triangle_feature = DTF
    D0.total_feature_length = calcul_total_feature_length(DTF)
    return L

def changement_base_2d(M0, M1, image_star_list):
    E1 = vectpp(M0.xy, M1.xy)
    E2 = (-E1[1], E1[0])
    k = norm2_sqr(E1)
    return [(etoile, (dot_prod2(vectpp(M0.xy, etoile.xy), E1)/k, dot_prod2(vectpp(M0.xy, etoile.xy), E2)/k)) for etoile in image_star_list]

def calcul_map_dtf_tfl_2d(M0, image_star_list):
    starandvect_list = [(M, vectpp(M0.xy, M.xy)) for M in image_star_list]
    best_4 = closest_nstars(starandvect_list, 4, False, 2)
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

    return F1, F2


def dtf_diff(dtf1, dtf2): #eq (6)
    return sum([abs(dtf1[i] - dtf2[i]) for i in range(12)])

def select_Ks_stars(data_base, etoile): #eq(7)

    l1 = etoile.F1_length
    l2 = etoile.F2_length
    A1 = min([(star.total_feature_length - l1) for star in data_base])*200
    A2 = min([(star.total_feature_length - l2) for star in data_base])*200

    data_base_reduced1 = []
    data_base_reduced2 = []
    for star in data_base:
        diff1 = (star.total_feature_length - l1)
        diff2 = (star.total_feature_length - l2)
        if diff1 < A1 and diff1<50:
            data_base_reduced1.append(star)
        if diff2 < A2 and diff2<50:
            data_base_reduced2.append(star)
    
    return data_base_reduced1, data_base_reduced2

def identify(D0, M0):
    r_score = 0
    matchlist = []
    for (catalog_star, (x1,y1)) in D0.gnomic_projection_map:
        for (etoile_image, (x2,y2)) in M0.normalized_map:
            if abs(x1-x2)<=0.05 and abs(y1-y2)<=0.05:
                r_score += 1
                matchlist.append((catalog_star, etoile_image))

    return r_score, matchlist

'''
CREATION BASE DE DONNEES
'''

csv_file = open(DATA_BASE_PATH)
next(csv_file) #skip la premiere ligne
reader = csv.reader(csv_file, delimiter=',')

DATA_BASE = [Star(intbis(row[0]),       #construction data_base, liste d'objets de type Star
                    intbis(row[1]),
                    strbis(row[2]),
                    strbis(row[3]),
                    strbis(row[4]),
                    strbis(row[5]),
                    fltbis(row[6]),
                    fltbis(row[7]),
                    fltbis(row[8]),
                    strbis(row[11])) for row in reader]
csv_file.close()

for star in DATA_BASE:
    calcul_gnomic_dtf_tfl(star,  DATA_BASE)


"""
TRAITEMENT IMAGE + TRUCS PYGAME
"""
image_original = Image.open(IMAGE_PATH).convert('L')
fps_pygame = 30
pygame.init()
pygame.font.init()
WINDOW = pygame.display.set_mode(image_original.size)
FONT = pygame.font.SysFont('Arial', 15)

def high_contrast(image, threshold):
    for x in range(image.width):
        for y in range(image.height):
            if image.getpixel((x,y))> threshold:
                image.putpixel((x,y), 255)
            else:
                image.putpixel((x,y), 0)

image=Image.open(IMAGE_PATH).convert('L')
high_contrast(image, 230)

def drawmap(map):
    img = Image.new('1',(300,300),0)
    for starandvect in map:
        x,y = round(starandvect[1][0]*30+150), round(starandvect[1][1]*30+150)
        if 0<=x<300 and 0<=y<300:
            img.putpixel((x,y), 1)
    img.show()

"""
CHOIX DU CONTRASTE
"""
'''
def high_contrast(image, threshold):
    for x in range(image.width):
        for y in range(image.height):
            if image.getpixel((x,y))> threshold:
                image.putpixel((x,y), 255)
            else:
                image.putpixel((x,y), 0)

choosing_contrast = True
threshold = 70
image = Image.open(IMAGE_PATH).convert('L')
high_contrast(image, threshold)

pygameblitfrompillow(image, window)
pygame.display.update()

while choosing_contrast:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            choosing_contrast = False
    if pygame.key.get_pressed()[pygame.K_UP]:
        threshold -= 5
    elif pygame.key.get_pressed()[pygame.K_DOWN]:
        threshold += 5
    if pygame.key.get_pressed()[pygame.K_UP] or pygame.key.get_pressed()[pygame.K_DOWN]:
        image = image_original.copy()
        high_contrast(image, threshold)
        pygameblitfrompillow(image, window)
        pygame.display.update()
    pygame.time.Clock().tick(fps_pygame)

image.convert("1")
#image.show()
'''

"""
DETECTION ETOILES
"""

def average_pix(pix_list): #renvoie le centroide d'une liste de pixels pris par etoile
    return (sum([pix[0] for pix in pix_list])/len(pix_list), sum([pix[1] for pix in pix_list])/len(pix_list))

def pix_list_of_star(image, pos, pixel_list): #détecte tous les pixels blancs "collés" à pos, les met en noir, et renvoie la liste de ces pixels
    pixel_list.append(pos)
    image.putpixel(pos, 0)
    x,y = pos
    width,height = image.size
    if x+1<width and image.getpixel((x+1,y)) == 255:
        pix_list_of_star(image, (x+1,y), pixel_list)
    if y+1<height and image.getpixel((x,y+1)) == 255:
        pix_list_of_star(image, (x,y+1), pixel_list)
    if 0<=x-1 and image.getpixel((x-1,y)) == 255:
        pix_list_of_star(image, (x-1,y), pixel_list)
    if 0<=y-1 and image.getpixel((x,y-1)) == 255:
        pix_list_of_star(image, (x,y-1), pixel_list)
    return pixel_list

def new_star(image, pos, star_list):
    star_list.append(Etoile_image(average_pix(pix_list_of_star(image, pos, [])))) #ajoute à star_list le centre de la liste des pixels d'une étoile détectée en pos


image_temp = image.copy() #utile car la fonction "new_star" a besoin de modifier l'image
LISTE_ETOILES_IMAGE = []
for x in range(image.width): #remplit la liste des coordonnéees des étoiles sur l'image
    for y in range(image.height):
        if image.getpixel((x,y)) == 255:
            new_star(image, (x,y), LISTE_ETOILES_IMAGE)

pygameblitfrompillow(image_original,WINDOW)

for etoile in LISTE_ETOILES_IMAGE:
    pygame.draw.circle(WINDOW, (0,0,255), round_xy(etoile.xy), 2)
pygame.display.update()

'''
IDENTIFICATION
'''

def choose_random(liste_etoiles, lambd):
    n = len(liste_etoiles)
    if(n<4):
        print("PAS ASSEZ D'ETOILES")
        return None
    elif (4<=n<lambd):
        return liste_etoiles
    elif (lambd<n):
        L=[]
        for i in range(n, n-lambd):
            L.append(liste_etoiles.pop(random.randint(0,i-1)))
        return L


def affiche_noms(liste_etoiles_image):
    for etoile in liste_etoiles_image:
        etoile.affiche(WINDOW)
    pygame.display.update()

'''
for star in DATA_BASE:
    drawmap(star.gnomic_projection_map)
'''

for etoile in LISTE_ETOILES_IMAGE:
    calcul_map_dtf_tfl_2d(etoile, LISTE_ETOILES_IMAGE)
for etoile in LISTE_ETOILES_IMAGE:
    print(len(select_Ks_stars(DATA_BASE, etoile)[0]))

while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()