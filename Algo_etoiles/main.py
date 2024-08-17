from math import *
from PIL import Image, ImageDraw, ImageFont
import csv
import random
from fpdf import FPDF

#import matplotlib.pyplot as plt

#LIRE https://www.mdpi.com/1424-8220/20/11/3027
#UTILISER https://pyfpdf.readthedocs.io/en/latest/index.html pour renvoyer un fichier pdf avec les hyperliens wikipedia/simbad

"""
CONSTANTES
"""

#LAMBD = 50 #nombre maximum d'etoiles de "reference"
IMAGE_PATH = "Algo_etoiles/images/test2.png"
DATA_BASE_PATH = "Algo_etoiles/databasecsv/treated_athyg_modified_vmagmax6.csv"
CENTROID_SAVE_PATH = "resultats/centroids.png"
RESULTS_SAVE_PATH = "resultats/results7.png"
RESULTSPDF_SAVE_PATH = "resultatspdf/results3.pdf"
FONT_PATH = "Algo_etoiles/arial.ttf"
ID_THRESHOLD = 0.1
DISPLAY_MODE = "mixt" #bayerflam / hip / mixt


"""
IMPORTATION ET TRAITEMENT D'IMAGE
"""
FONT_SIZE = 8
DISPLAY_CIRCLE_RADIUS = 4 #rayon des cercles lors de l'affichage des etoiles trouvees
FONT = ImageFont.truetype(font = FONT_PATH, size = FONT_SIZE)

BLACK_WHITE_THRESHOLD = 100
def pixel_to_NB(pixel): #permet de choisir le contraste de la conversion en noir et blanc
    return 1 if pixel > BLACK_WHITE_THRESHOLD else 0

IMAGE_ORIGINAL = Image.open(IMAGE_PATH)
IMAGE_GRAYSCALE = IMAGE_ORIGINAL.convert('L')
IMAGE_NB = IMAGE_GRAYSCALE.point(pixel_to_NB, mode='1')

SIZE = WIDTH, HEIGHT = IMAGE_ORIGINAL.size

"""
TYPES
"""
class Star: #Type enregistrement pour les étoiles de la base de données

    def __init__(self, id, hip, bayer, flam, con, proper, ra, dec, mag, unicode, x, y, z, tfl, dtf, wikipedia, simbad): #ra et dec en heures et degres dans le csv
        self.id = id            #Identifiant dans la base de donnéee
        self.hip = hip          #Identifiant hipparcos (si il existe)
        self.bayer = bayer      #Désignation de bayer  (si elle existe)
        self.flam = flam        #Désignation de flamsteed (si elle existe)
        self.con = con          #Abréviation de la constellation du système international
        self.proper = proper    #Nom commun (si il existe)
        self.ra = ra*pi/12      #Ascention droite avec conversion heures(15deg)->radians
        self.dec = dec*pi/180   #Déclinaison avec conversion degrés->radians
        self.mag = mag          #Magnitude (correspond en gros à du -log(luminosité))
        self.xyz = (self.x, self.y, self.z) = x, y, z #(cosdec*cosra, cosdec*sinra, sindec) #Position sur la sphere celeste unité
        self.total_feature_length = tfl        #Longueur totale des segments de double_triangle_feature
        self.double_triangle_feature = dtf     #Longueurs et angles des deux triangles formé par les 4 etoiles les plus proches
        self.uni = unicode      #Caractère unicode grec de la désignation de bayer
        self.wiki = wikipedia   #Lien vers la page wikipedia (si elle existe)
        self.simbad = simbad    #Lien vers la page simbad

        #cdec, sdec, cra, sra = cos(self.dec), sin(self.dec), cos(self.ra), sin(self.ra)

        self.imagematch = None                  #Objet Etoile_Image correspondant si trouvé
        self.gnomic_projection_map = None      #Etoiles proches et leurs coordonnees dans la base adaptée à cet étoile

        self.display_name = None    #Nom à afficher si trouvée (bayer, flamsteed ou id)
        if DISPLAY_MODE == "bayerflam" and self.bayer != None:
            self.display_name = self.bayer + " " + self.con
        elif DISPLAY_MODE == "bayerflam" and self.flam != None:
            self.display_name = self.flam + " " + self.con
        elif DISPLAY_MODE == "hip" and self.hip != None:
            self.display_name = "HIP"+str(self.hip)
        elif DISPLAY_MODE == "mixt":
            if self.proper != None:
                self.display_name = self.proper
            elif self.bayer != None:
                self.display_name = self.bayer + " " + self.con
            elif self.flam != None: 
                self.display_name = self.flam + " " + self.con
            elif self.hip != None:
                self.display_name = "HIP"+str(self.hip)
        else:
            self.display_name = str(self.id)

    def load_gnomic(self):
        path = "Algo_etoiles/treated_data/"+str(self.id)+".csv"
        self.gnomic_projection_map = []
        with open(path, newline='') as file:
            reader = csv.reader(file, delimiter=',')
            for row in reader:
                starid, x, y = int(row[0]), float(row[1]), float(row[2])
                self.gnomic_projection_map.append((starid,(x,y)))


    def draw_name_pillow(self, drawing_instance):
        if self.imagematch.xy != None:
            x, y = self.imagematch.xy
            r = DISPLAY_CIRCLE_RADIUS
            if y-r-4-FONT_SIZE>0:
                textxy = (x, y-r-4)
                mode = "ms" #middle baseline
            else:
                textxy = (x, y+r+4)
                mode = "mt" #middle top
            circle_x0y0x1y1 = (x-r, y-r, x+r, y+r)

            drawing_instance.text(textxy, self.display_name, font = FONT, fill = (0, 255, 0, 255), anchor = mode)
            drawing_instance.ellipse(circle_x0y0x1y1, outline = (0, 255, 0, 255), width = 1)
        
class Etoile_image: #Type enregistrement pour les étoiles de l'image

    def __init__(self, pos):
        self.xy = self.x, self.y = pos

        self.F1 = None                  #Première feature
        self.F1_length = None           #Longueur totale des segments de F1
        self.F2 = None                  #Deuxième feature
        self.F2_length = None           #Longueur totale des segments de F2
        self.normalized_map = None      #Etoiles proches et leurs coordonnees dans la base adaptée à cet étoile

        self.starmatch = None           #Objet Star correspondant si trouvé

    def draw_name_pillow(self, drawing_instance):
        x, y = self.xy
        r = DISPLAY_CIRCLE_RADIUS
        if y-r-4-FONT_SIZE>0:
            textxy = (x, y-r-4)
            mode = "ms" #middle baseline
        else:
            textxy = (x, y+r+4)
            mode = "mt" #middle top
        circle_x0y0x1y1 = (x-r, y-r, x+r, y+r)

        if self.starmatch != None:
            #print(self.xy, self.starmatch)
            drawing_instance.text(textxy, self.starmatch.display_name, font = FONT, fill = (0, 255, 0, 255), anchor = mode)
            drawing_instance.ellipse(circle_x0y0x1y1, outline = (0, 255, 0, 255), width = 1)
        else:
            #print(self.xy, None)
            drawing_instance.ellipse(circle_x0y0x1y1, outline = (255, 0, 0, 255), width = 1)
            
    def draw_name_pdf(self, pdf):
        x, y = self.xy
        r = DISPLAY_CIRCLE_RADIUS
        textx, texty = x-2.5*FONT_SIZE, y-r-4-FONT_SIZE
        if texty<0:
            texty = y+r+4
        circlex, circley = x-r, y-r

        if self.starmatch != None:
            #print(self.xy, self.starmatch)
            pdf.set_xy(textx, texty)
            pdf.set_text_color(0, 255, 0)
            pdf.set_draw_color(0, 255, 0)
            pdf.write(5, self.starmatch.display_name, self.starmatch.simbad)
        else:
            pdf.set_draw_color(255, 0, 0)
        pdf.ellipse(circlex, circley, 2*r, 2*r)



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
def unidecode(s):
    return None if (s=='' or s=="#N/A") else bytes(s).decode("utf-8")

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

def changement_base_2d(M0, M1, image_star_list):
    #Normalisation de la liste des coordonnees des etoiles en fonction de deux points/etoiles de l'image
    E1 = vectpp(M0.xy, M1.xy)
    E2 = (-E1[1], E1[0])
    k = norm2_sqr(E1)
    return [(etoile, (dot_prod2(vectpp(M0.xy, etoile.xy), E1)/k, dot_prod2(vectpp(M0.xy, etoile.xy), E2)/k)) for etoile in image_star_list]

def calcul_map_dtf_tfl_2d(M0, image_star_list):
    '''
    Calculs à faire sur chaque étoile (ou une partie des étoiles) de l'image
    '''

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

    return


def dtf_diff(dtf1, dtf2): #cf. eq(6) de https://www.mdpi.com/1424-8220/20/11/3027
    return sum([abs(dtf1[i] - dtf2[i]) for i in range(12)])
'''
def select_Ks_stars(data_base, etoile): #eq(7) MARCHE TRES TRES MAL

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
'''

def identify(D0, M0):
    '''
        Une fois que des etoiles D0 et M0 semblent correspondre via la methode des doubles triangles,
        on compare leurs projections normalisees afin d'établir des correspondances,
        le résuultat final est l'identification avec le meilleur "r_score": le nombre d'étoiles identifiées.
    '''
    D0.load_gnomic() #besoin d'importer ces valeurs seulement pour les etoiles choisies par la "premiere selction"
    r_score = 0
    matchlist = []
    for (etoile_image, (x2,y2)) in M0.normalized_map:
        for (catalog_star_id, (x1,y1)) in D0.gnomic_projection_map:
            if abs(x1-x2)<=ID_THRESHOLD and abs(y1-y2)<=ID_THRESHOLD:
                r_score += 1
                matchlist.append((catalog_star_id, etoile_image))
                break #a ameliorer

    return r_score, matchlist

'''
IMPORTATION ET CALCULS SUR LA BASE DE DONNEES
'''

csv_file = open(DATA_BASE_PATH)
next(csv_file) #skip la premiere ligne
reader = csv.reader(csv_file, delimiter=',')

reader = csv.reader(csv_file, delimiter=',')
DATA_BASE = [Star(                      #construction de DATA_BASE: liste d'objets de type Star, conversions un peu inutiles mais copiees
                int(row[0]),            #id
                intbis(row[1]),         #hip
                strbis(row[2]),         #bayer
                strbis(row[3]),         #flam
                str(row[4]),            #const
                strbis(row[5]),         #proper
                float(row[6]),          #ra
                float(row[7]),          #dec
                float(row[8]),          #mag
                strbis(row[11]),        #greek unicode
                str(row[12]),           #x
                str(row[13]),           #y
                str(row[14]),           #z
                str(row[15]),           #tfl
                [float(row[i]) for i in range(16,28)],     #dtf
                strbis(row[28]),        #wiki
                strbis(row[29]),        #simbad
                ) for row in reader]

csv_file.close()

"""
TRAITEMENT IMAGE
"""

def drawmap(map):
    #affiche un apperçu des etoiles dans le repere normalisé autour d'une etoile
    img = Image.new('1',(300,300),0)
    for starandvect in map:
        x,y = round(starandvect[1][0]*30+150), round(starandvect[1][1]*30+150)
        if 0<=x<300 and 0<=y<300:
            img.putpixel((x,y), 1)
    img.show()

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
    if x+1<width and image.getpixel((x+1,y)) == 1:
        pix_list_of_star(image, (x+1,y), pixel_list)
    if y+1<height and image.getpixel((x,y+1)) == 1:
        pix_list_of_star(image, (x,y+1), pixel_list)
    if 0<=x-1 and image.getpixel((x-1,y)) == 1:
        pix_list_of_star(image, (x-1,y), pixel_list)
    if 0<=y-1 and image.getpixel((x,y-1)) == 1:
        pix_list_of_star(image, (x,y-1), pixel_list)
    return pixel_list

def new_star(image, pos, star_list):
    star_list.append(Etoile_image(average_pix(pix_list_of_star(image, pos, [])))) #ajoute à star_list le centre de la liste des pixels d'une étoile détectée en pos


image_temp = IMAGE_NB.copy() #utile car la fonction "new_star" a besoin de modifier l'image
LISTE_ETOILES_IMAGE = []
for y in range(HEIGHT): #remplit la liste des coordonnéees des étoiles sur l'image
    for x in range(WIDTH):
        if image_temp.getpixel((x,y)) == 1:
            new_star(image_temp, (x,y), LISTE_ETOILES_IMAGE)

#AFFICHAGE TEMPORAIRE

new_layer = Image.new("RGBA", SIZE, (255, 255, 255, 0)) #crée une image transparente pour afficher le texte sur fond transparent avant de le combiner à image_original
drawing_instance = ImageDraw.Draw(new_layer)
image_copy = IMAGE_ORIGINAL.copy().convert('RGBA')

for etoile in LISTE_ETOILES_IMAGE: #affichage
    drawing_instance.rectangle((round(etoile.x)-1, round(etoile.y)-1, round(etoile.x)+1, round(etoile.y)+1), (0,0,255,255), (0,0,255,255), width=1)

centroids = Image.alpha_composite(image_copy, new_layer)
centroids.show("Centroïdes des étoiles repérées")
centroids.save(CENTROID_SAVE_PATH)


'''
IDENTIFICATION
'''

def choose_random(liste_etoiles, lambd): #choisit "lambd" etoiles au hasard sur l'image
    n = len(liste_etoiles)
    if(n<4):
        print("PAS ASSEZ D'ETOILES")
        return None
    elif (4<=n<lambd):
        return liste_etoiles
    elif (lambd<n):
        L=[]
        for i in range(n, n-lambd): #???? je fous quoi
            L.append(liste_etoiles.pop(random.randint(0,i-1))) #detruit liste_etoiles ??? a refaire
        return L

def affiche_noms_pillow(liste_etoiles_image): #affiche le nom de l'etoile identifiee pour chaque etoile de l'image
    text_layer = Image.new("RGBA", SIZE, (255, 255, 255, 0)) #cree une image transparente pour afficher le texte sur fond transparent avant de le combiner à image_original
    drawing_instance = ImageDraw.Draw(text_layer)
    image_copy = IMAGE_ORIGINAL.copy().convert('RGBA')
    
    for etoile in liste_etoiles_image:
        etoile.draw_name_pillow(drawing_instance)
    
    return Image.alpha_composite(image_copy, text_layer)

def affiche_noms_pdf(liste_etoiles_image): #affiche le nom de l'etoile identifiee pour chaque etoile de l'image
    pdf = FPDF("P", "pt", (WIDTH, HEIGHT))
    pdf.add_font("police", "", FONT_PATH, uni=True)
    pdf.set_font("police", "", FONT_SIZE)
    pdf.set_margins(0, 0, 0)
    pdf.set_auto_page_break(False)
    pdf.add_page()
    pdf.image(IMAGE_PATH, 0, 0, WIDTH, HEIGHT)
    for etoile in liste_etoiles_image:
        etoile.draw_name_pdf(pdf)
    
    return pdf

def closest_dtf(dtf, data_base): #renvoie l'objet Star de data_base avec la feature la plus proche de dtf
    best = data_base[0]
    mindiff = dtf_diff(best.double_triangle_feature, dtf)
    for star in data_base:
        diff = dtf_diff(star.double_triangle_feature, dtf)
        if diff < mindiff:
            best, mindiff = star, diff
    return best

def get_by_attribute(data_base, attribute, val): #sûrement possible de faire plus proprement mais pas grave
    if attribute == "bayer":
        return next((star for star in data_base if star.bayer == val), [None])
    elif attribute == "id":
        return next((star for star in data_base if star.id == val), [None])
    elif attribute == "hip":
        return next((star for star in data_base if star.hip == val), [None])


'''
for star in DATA_BASE:
    drawmap(star.gnomic_projection_map)
'''

#LISTE_ETOILES_REF = choose_random(LISTE_ETOILES_IMAGE, LAMBD)
LISTE_ETOILES_REF = LISTE_ETOILES_IMAGE

for etoile in LISTE_ETOILES_REF:
    calcul_map_dtf_tfl_2d(etoile, LISTE_ETOILES_IMAGE)
    #drawmap(etoile.normalized_map)


bestr_score, bestmatchlist = -1, []
for etoile in LISTE_ETOILES_REF:
    D0 = closest_dtf(etoile.F1, DATA_BASE)
    r_score, matchlist = identify(D0, etoile)
    if r_score > bestr_score:
        bestr_score, bestmatchlist = r_score, matchlist

for (starid, etoile) in bestmatchlist:
    star = get_by_attribute(DATA_BASE, "id", starid)
    star.imagematch = etoile
    etoile.starmatch = star

results = affiche_noms_pillow(LISTE_ETOILES_IMAGE)
results.show("Étoiles reconnues")
#results.save(RESULTS_SAVE_PATH)


resultspdf = affiche_noms_pdf(LISTE_ETOILES_IMAGE)
resultspdf.output(RESULTSPDF_SAVE_PATH, 'F')


#Del = get_by_attribute(DATA_BASE, "bayer", "Del")
#drawmap(Del.gnomic_projection_map)