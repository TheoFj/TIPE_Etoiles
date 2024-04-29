from math import *
from PIL import Image
import csv
#import matplotlib.pyplot as plt
import pygame
import sys

"""
CONSTANTES
"""
GRID_SIZE = 30
n=2

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

"""
FONCTIONS PRATIQUES
"""
def intbis(s):
    return None if (s=='') else int(s)
def strbis(s):
    return None if (s=='') else str(s)
def fltbis(s):
    return None if (s=='') else float(s)
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
    return acos(sin(Star1.dec)*sin(Star2.dec) + cos(Star1.dec)*cos(Star2.dec)*cos(Star1.ra - Star2.ra))

def flat_dist(pos1, pos2):
    return sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
def flat_dist_sqr(pos1, pos2):
    return (pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2

def dot_prod3(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
def norm3(v):
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)


with open("Algo_etoiles/database/UMa_reduced.csv") as csv_file:
    next(csv_file) #skip la premiere ligne
    reader = csv.reader(csv_file, delimiter=',')

    data_base = [Star(intbis(row[0]),       #construction data_base, liste d'objets de type Star
                      intbis(row[1]),
                      strbis(row[2]),
                      intbis(row[3]),
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
        #p = dot_prod3(cs.xyz, star.xyz)
        #(xp, yp, zp) = (star.x/p, star.y/p, star.z/p) #projection sur le plan perpendiculaire à l'etoile cs sur la sphere unite
        X = dot_prod3(star.xyz, cs.etheta) #coordonnees de la projection dans la base etheta ephi du plan
        Y = dot_prod3(star.xyz, cs.ephi)
        list.append((X,Y))
    return list

L = gnomic_projection(data_base[n], data_base)
print(L)

fps_pygame = 30
pygame.init()
window = pygame.display.set_mode((1000,1000))
pygame.Surface.fill(window, 0)

for pos in L:
    pygame.draw.circle(window, (255,255,255), (500+ round(1000*pos[0]), 500+round(1000*pos[1])), 5)

pygame.display.update()


def closest_star(pos,list): #renvoie l'indice de l'etoile/pixel dans la liste la plus proche de pos
    closest_index = 0
    smallestdist = flat_dist_sqr(list[0], pos)
    for i in range(len(list)):
        dist = flat_dist_sqr(list[i], pos)
        if (0 < dist < smallestdist):
            closest_index = i
            smallestdist = dist
    return closest_index

def generate_grid(m, n, list, gsize): #m: etoile milieu, n: etoile de calibration (normalement la plus proche de m)
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

G = generate_grid(n, closest_star(L[n], L), L, GRID_SIZE)
G.show()



while True:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()