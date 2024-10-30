from math import *
import csv
import config

def display_name_calc(proper, con, bayer, flam, hip, id):
    display_name = ""
    if config.DISPLAY_MODE == "bayerflam" and bayer != None:
        display_name = bayer + " " + con
    elif config.DISPLAY_MODE  == "bayerflam" and flam != None:
        display_name = flam + " " + con
    elif config.DISPLAY_MODE  == "hip" and hip != None:
        display_name = "HIP" + str(hip)
    elif config.DISPLAY_MODE  == "mixt":
        if proper != None:
            display_name = proper
        elif bayer != None:
            display_name = bayer + " " + con
        elif flam != None: 
            display_name = flam + " " + con
        elif hip != None:
            display_name = "HIP"+str(hip)
    else:
        display_name = str(id)
    return display_name


class Star: #Type enregistrement pour les étoiles de la base de données

    def __init__(self, id, hip, bayer, flam, con, proper, ra, dec, mag, greek_bay, x, y, z, tfl, dtf, wikipedia, simbad): #ra et dec en heures et degres dans le csv
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
        self.greek_bay = greek_bay      #Caractère unicode grec de la désignation de bayer
        self.wiki = wikipedia   #Lien vers la page wikipedia (si elle existe)
        self.simbad = simbad    #Lien vers la page simbad

        #cdec, sdec, cra, sra = cos(self.dec), sin(self.dec), cos(self.ra), sin(self.ra)

        self.imagematch = None                  #Objet Etoile_Image correspondant si trouvé
        self.gnomic_projection_map = None      #Etoiles proches et leurs coordonnees dans la base adaptée à cet étoile

        if config.GREEK_DISPLAY:
            self.display_name = display_name_calc(proper, con, greek_bay, flam, hip, id)    #Nom à afficher si trouvée (bayer, flamsteed ou id)
        else:
            self.display_name = display_name_calc(proper, con, bayer, flam, hip, id)    #Nom à afficher si trouvée (bayer, flamsteed ou id)

        
    def load_gnomic(self):
        path = config.GNOMIC_PATH + str(self.id) + ".csv"
        self.gnomic_projection_map = []
        with open(path, newline='') as file:
            reader = csv.reader(file, delimiter=',')
            for row in reader:
                starid, x, y = int(row[0]), float(row[1]), float(row[2])
                self.gnomic_projection_map.append((starid,(x,y)))


class Etoile_image: #Type enregistrement pour les étoiles de l'image

    def __init__(self, pos):
        self.xy = self.x, self.y = pos

        self.F1 = None                  #Première feature
        self.F1_length = None           #Longueur totale des segments de F1
        self.F2 = None                  #Deuxième feature
        self.F2_length = None           #Longueur totale des segments de F2
        self.normalized_map = None      #Etoiles proches et leurs coordonnees dans la base adaptée à cet étoile
        self.closest_star = None

        self.starmatch = None           #Objet Star correspondant si trouvé
'''
    def draw_name_pillow(self, drawing_instance):
        x, y = self.xy
        r = config.DISPLAY_CIRCLE_RADIUS
        
        if y-r-4-config.FONT_SIZE>0:
            textxy = (x, y-r-4)
            mode = "ms" #middle baseline
        else:
            textxy = (x, y+r+4)
            mode = "mt" #middle top
        circle_x0y0x1y1 = (x-r, y-r, x+r, y+r)

        if self.starmatch != None:
            #print(self.xy, self.starmatch)
            drawing_instance.text(textxy, self.starmatch.display_name, font = config.FONT, fill = (0, 255, 0, 255), anchor = mode)
            drawing_instance.ellipse(circle_x0y0x1y1, outline = (0, 255, 0, 255), width = 1)
        else:
            #print(self.xy, None)
            drawing_instance.ellipse(circle_x0y0x1y1, outline = (255, 0, 0, 255), width = 1)
            
'''