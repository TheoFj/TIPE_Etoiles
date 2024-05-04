from math import *
from PIL import Image
import csv
#import matplotlib.pyplot as plt
import pygame
import sys

"""
CONSTANTES
"""
GRID_SIZE = 50
FOV = pi/2

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
        self.etheta = (-sra, cra, 0)
        self.ephi = (-sdec*cra, -sdec*sra, cdec)

        self.grid = None
        self.imagexy = None
        self.best_fit = 0

        self.bf = None
        if self.bayer != None:
            self.bf = self.bayer + " " + self.con
        elif self.flam != None: 
            self.bf = self.flam + " " + self.con
        

"""
FONCTIONS PRATIQUES
""" 

def intbis(s):
    return None if (s=='' or s=='#N/A') else int(s)
def strbis(s):
    return None if (s=='' or s=='#N/A') else str(s)
def fltbis(s):
    return None if (s=='' or s=='#N/A') else float(s)
def round_xy(coord):
    return (round(coord[0]), round(coord[1]))
def vectpp(p1, p2):
    return (p2[0] - p1[0], p2[1] - p1[1])

def fit(pat_d, pat_s, grid_size): #pat_d et pat_s sont des matrices de 0 et de 1
    score = 0
    for i in range(grid_size):
        for j in range(grid_size):
            score += pat_d[i][j]*pat_s[i][j]
    return score

def ang_dist(Star1, Star2): #https://en.wikipedia.org/wiki/Angular_distance
    p = dot_prod3(Star1.xyz, Star2.xyz)
    if p>1: #gere les imprecisions de calcul genre 1.0000000000000002
        return 0
    elif p<-1:
        return pi
    return acos(dot_prod3(Star1.xyz, Star2.xyz))

def flat_dist(pos1, pos2):
    return sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
def flat_dist_sqr(pos1, pos2):
    return (pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2

def dot_prod3(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
def norm3(v):
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

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

def gnomic_projection(cs, star_list):#cs: center_star
    list = []
    for star in star_list:
        if ang_dist(cs, star) <= FOV:
            p = dot_prod3(cs.xyz, star.xyz)
            #(xp, yp, zp) = (star.x/p, star.y/p, star.z/p) #projection sur le plan perpendiculaire à l'etoile cs sur la sphere unite
            X = dot_prod3(star.xyz, cs.etheta)/p #coordonnees de la projection dans la base etheta ephi du plan
            Y = dot_prod3(star.xyz, cs.ephi)/p    
            list.append((X,Y))
    return list

def closest_star(pos, list, with_zero): #renvoie l'indice de l'etoile/pixel dans la liste la plus proche de pos
    if flat_dist_sqr(list[0], pos) != 0:
        closest_index = 0
    else: #gere le cas où pos correspond aux coordonnees de l'etoile 0 dans list
        closest_index = 1
    smallestdist = flat_dist_sqr(list[closest_index], pos)
    for i in range(len(list)):
        dist = flat_dist_sqr(list[i], pos)
        if (dist < smallestdist) and (with_zero or dist>0): #with_zero sert à demander si on compte une etoile pile sur pos ou non
            closest_index = i
            smallestdist = dist
    return closest_index

def generate_grid_image(m, n, list, gsize): #m: etoile milieu, n: etoile de calibration (normalement la plus proche de m)
    distsqr = flat_dist_sqr(list[m], list[n])
    xn, yn = vectpp(list[m], list[n]) #vecteur etoile centre->etoile plus proche
    grid = Image.new("1", (gsize, gsize), 0)
    for i in range(len(list)):
        if (i!=m):
            x, y = vectpp(list[m], list[i])
            xgrid = gsize//2 + round(4*(1/distsqr)*(xn*x + yn*y)) #changement de base en multipliant (x,y) par la matrice inverse pour exprimer leurs coord en fonction de x
            ygrid = gsize//2 + round(4*(1/distsqr)*(-yn*x + xn*y))
            if (0 <= xgrid <= gsize-1 and 0 <= ygrid <= gsize-1):
                grid.putpixel((xgrid, ygrid), 1)
    return grid

def generate_grid_matrix(m, n, list, gsize): #m: etoile milieu, n: etoile de calibration (normalement la plus proche de m)
    distsqr = flat_dist_sqr(list[m], list[n])
    xn, yn = vectpp(list[m], list[n]) #vecteur etoile centre->etoile plus proche
    grid = [[0 for j in range(gsize)] for j in range(gsize)]
    for i in range(len(list)):
        if (i!=m):
            x, y = vectpp(list[m], list[i])
            xgrid = gsize//2 + round(4*(1/distsqr)*(xn*x + yn*y)) #changement de base en multipliant (x,y) par la matrice inverse pour exprimer leurs coord en fonction de x
            ygrid = gsize//2 + round(4*(1/distsqr)*(-yn*x + xn*y))
            if (0 <= xgrid <= gsize-1 and 0 <= ygrid <= gsize-1):
                grid[xgrid][ygrid] = 1
    return grid

for star in DATA_BASE:  #CALCUL DES GRILLES DE LA BASE DE DONNEES
    L = gnomic_projection(star, DATA_BASE)
    m = closest_star((0,0), L, True)
    n = closest_star(L[m], L, False)
    star.grid = generate_grid_matrix(m, n, L, GRID_SIZE)


"""
TRAITEMENT IMAGE + TRUCS PYGAME
"""
image_path = "Algo_etoiles/images/UMa2.png"
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
threshold = 240
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
def get_by_bayer(bayer_l, data_base):
    output = []
    for bayer in bayer_l:
        for star in data_base:
            if star.bayer == bayer:
                output.append(star)
                break
    return output

#LISTE_MATCHS = []
N = len(liste_etoiles_image)
for m in range(N):
    n = closest_star(liste_etoiles_image[m], liste_etoiles_image, False)
    grille = generate_grid_matrix(m, n, liste_etoiles_image, GRID_SIZE)
    matching = None
    imagestar_best_fit = -1
    for star in DATA_BASE:
        f = fit(star.grid, grille, GRID_SIZE)
        if f>imagestar_best_fit and f>star.best_fit:
            matching = star
            star.imagexy = liste_etoiles_image[m]
            star.best_fit = f
            imagestar_best_fit = f
    
    #LISTE_MATCHS.append((liste_etoiles_image[m], matching))

'''
def affiche_noms(liste_matchs):
    for match in liste_matchs:
        if match[1] != None:
            print(LISTE_MATCHS.index(match), match[0], match[1].bf)
            pygame.draw.circle(window, (0,255,0), match[0], 8, 1)
            window.blit(my_font.render(match[1].bf, False, (0, 255, 0)), (match[0][0]-25, match[0][1]-25))
        else:
            print(LISTE_MATCHS.index(match), match[0], match[1])
            pygame.draw.circle(window, (255,0,0), match[0], 8, 1)
    pygame.display.update()
'''

def affiche_noms2(data_base):
    for star in data_base:
        if star.imagexy != None:
            print(star.imagexy, star.bf)
            pygame.draw.circle(window, (0,255,0), star.imagexy, 8, 1)
            window.blit(my_font.render(star.bf, False, (0, 255, 0)), (star.imagexy[0]-25, star.imagexy[1]-25))
    
    pygame.display.update()

affiche_noms2(DATA_BASE)
print(len(DATA_BASE))
print(len(liste_etoiles_image))

while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()