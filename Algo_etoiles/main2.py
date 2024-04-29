from math import *
from PIL import Image
import csv
#import matplotlib.pyplot as plt
import pygame
import sys

"""
CONSTANTES
"""
FOV = 90
L = FOV*pi/180 #fov en radians
L2 = (L**2)/2
lamdda = 50 #nombre maximum d'etoiles de "reference"

#LIRE https://www.mdpi.com/1424-8220/20/11/3027#sec4dot3-sensors-20-03027 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

        self.imagexy = None
        self.double_triangle_feature = None
        self.total_feature_length = None
        self.gnomic = None

        self.display_name = None
        if self.bayer != None:
            self.display_name = self.bayer + " " + self.con
        elif self.flam != None: 
            self.display_name = self.flam + " " + self.con

    hip = 123

    def get_hip(self):
        return self.hip
        

"""
FONCTIONS PRATIQUES
""" 

def intbis(s):
    return None if (s=='' or s=='#N/A') else int(s)
def strbis(s):
    return None if (s=='' or s=='#N/A') else str(s)
def fltbis(s):
    return None if (s=='' or s=='#N/A') else float(s)

def fit(pat_d, pat_s, grid_size): #pat_d et pat_s sont des matrices de 0 et de 1
    score = 0
    for i in range(grid_size):
        for j in range(grid_size):
            score += pat_d[i][j]*pat_s[i][j]
    return score

def ang_dist(Star1, Star2): #https://en.wikipedia.org/wiki/Angular_distance
    p = dot_prod3(Star1.xyz, Star2.xyz)
    if p>1: #gere les imprecisions de calcul genre 1.0000000000000002
        return 0.0
    elif p<-1:
        return pi
    return acos(p)

def angle(v1,v2, norm_prod):
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


'''
CREATION BASE DE DONNEES
'''

with open("Algo_etoiles/database/UMa_vmagmax5.csv") as csv_file:
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

def closest_nstars(star_list, n, with_zero, dim): #O(len(star_list)*n²) ATTENTION: compare par rapport à (0,0)/(0,0,0) pour comparer à un autre point, mettre en entree les vecteurs "P0->P1"
    if dim==2: L_dist = [(norm2_sqr(star), star) for star in star_list]
    elif dim==3: L_dist = [(norm3_sqr(star), star) for star in star_list]
    else: print("????????????")
    best_n = []

    for couple in L_dist:
        best_n.append(couple)
        best_n.sort(key = tuple[0])
        best_n = best_n[:(n+1)]

    if with_zero:
        return [best_n[i][1] for i in range(0,n)] #on ne renvoie pas le dernier, indicé n, car on ne renvoie que les n meilleures
    else:
        return [best_n[i][1] for i in range(1,n+1)] #on ne renvoie pas le numero 0 car il s'agit de l'etoile "pos"
    
def calcul_double_triangle_feature(M0,M1,M2,M3): #appelé sur M0,M2,M3,M4 pour construire F2
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
    return (angle(v21, v20, n21*n02),
            angle(v10, v12, n21*n01), 
            angle(v01, v02, n01*n02), 
            angle(v30, v32, n03*n32), 
            angle(v23, v20, n32*n02), 
            angle(v03, v02, n03*n02), 
            n01, 
            n02, 
            n21, 
            n02,
            n03, 
            n32, 
            )

def gnomic_projection(D0, star_list):
    
    #PROJECTION SUR LE PLAN TANGEANT AU CERCLE UNITE
    C = []
    for star in star_list:
        p = dot_prod3(D0.xyz, star.xyz)
        proj_vect = (star.x/p - D0.x), (star.y/p - D0.y), (star.z/p - D0.z)
        if p>0 and norm3_sqr(proj_vect) < L2:
            C.append(proj_vect)

    #CHANGEMENT DE BASE EN PRENANT LETOILE LA PLUS PROCHE
    best_3 = closest_nstars(C, 3, False, 3)
    A = cross_prod3(D0.xyz, best_3[0]) #rotation de pi/2 de best_3[0]
    K = norm3_sqr(best_3[0])

    L = [(dot_prod3(coord, best_3[0])/K, dot_prod3(coord, best_3[0])/K) for coord in C]
    C0 = (0., 0.)
    C1 = (1., 0.) #projection de best_3[0]
    C2 = (dot_prod3(best_3[1], best_3[0])/K, dot_prod3(best_3[1], A)/K)
    C3 = (dot_prod3(best_3[2], best_3[0])/K, dot_prod3(best_3[2], A)/K)
    DTF = calcul_double_triangle_feature(C0, C1, C2, C3)

    D0.gnomic = L
    D0.double_triangle_feature = DTF
    D0.total_feature_length = calcul_total_feature_length(DTF)
    return L

def changement_base_2d(M0, M1, liste_etoiles):
    E1 = vectpp(M0, M1)
    E2 = (-E1[1], E1[0])
    K = norm2_sqr(E1)
    return [(dot_prod2(vectpp(M0, P), E1)/K, dot_prod2(vectpp(M0, P), E2)/K) for P in liste_etoiles]

"""
TRAITEMENT IMAGE + TRUCS PYGAME
"""
'''
image_path = "Algo_etoiles/images/UMa.png"
image_original = Image.open(image_path).convert('L')
fps_pygame = 30
pygame.init()
pygame.font.init()
window = pygame.display.set_mode(image_original.size)
my_font = pygame.font.SysFont('Arial', 15)

"""
CHOIX DU CONTRASTE
"""

def high_contrast(image, threshold):
    for x in range(image.width):
        for y in range(image.height):
            if image.getpixel((x,y))> threshold:
                image.putpixel((x,y), 255)
            else:
                image.putpixel((x,y), 0)

def pygameblitfrompillow(image, window):
    image.save("temp.png")
    im_pygame = pygame.image.load("temp.png")
    window.blit(im_pygame,(0,0))


choosing_contrast = True
threshold = 70
image = Image.open(image_path).convert('L')
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

"""
DETECTION ETOILES
"""

def average_pix(pix_list):
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
    star_list.append(average_pix(pix_list_of_star(image, pos, []))) #ajoute à star_list le centre de la liste des pixels d'une étoile détectée en pos
    
image_temp = image.copy()
liste_etoiles_image = []
for x in range(image.width): #remplit la liste des coordonnéees des étoiles sur l'image
    for y in range(image.height):
        if image.getpixel((x,y)) == 255:
            new_star(image, (x,y), liste_etoiles_image)

pygameblitfrompillow(image_original,window)
for etoile in liste_etoiles_image:
    pygame.draw.circle(window, (0,0,255), round_xy(etoile), 2)
pygame.display.update()

"""
MATCHING
"""
'''

def dtf_diff(dtf1, dtf2): #eq (6)
    return sum([abs(dtf1[i] - dtf2[i]) for i in range(12)])

def calcul_total_feature_length(feature):
    return sum([feature[i] for i in range(6,12)])

def select_Ks_stars(data_base, feature_etoile): #eq(7)
    lFr = calcul_total_feature_length(feature_etoile)
    A = 200*min([(star.total_feature_length - lFr) for star in data_base])

    data_base_reduced = []
    for star in data_base:
        diff = (star.total_feature_length - lFr)
        if diff < A and diff<50:
            data_base_reduced.append(star)
    
    return data_base_reduced



'''
def get_by_bayer(bayer_l, data_base):
    pass
    return data_base.index(bayer_l, key = Star.bayer)

while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()
'''