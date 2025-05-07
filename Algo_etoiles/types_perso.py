from math import *
import csv
import config

class Star: #Type enregistrement pour les étoiles de la base de données

    def __init__(self, id, hip, bayer, flam, con, proper, ra, dec, mag, full, gen, greek_bay, x, y, z, tfl, dtf): #ra et dec en heures et degres dans le csv
        self.id = id                    #Identifiant dans la base de donnéee
        self.hip = hip                  #Identifiant hipparcos (si il existe)
        self.bayer = bayer              #Désignation de bayer  (si elle existe)
        self.flam = flam                #Désignation de flamsteed (si elle existe)
        self.con = con                  #Abréviation de la constellation du système international
        self.proper = proper            #Nom commun (si il existe)
        self.ra = ra*pi/12              #Ascention droite avec conversion heures(15deg)->radians
        self.dec = dec*pi/180           #Déclinaison avec conversion degrés->radians
        self.mag = mag                  #Magnitude (correspond en gros à du -log(luminosité))
        self.full_letter = full         #Désignation de bayer en toutes lettres
        self.genitive = gen             #Génitif latin de la constellation en toutes lettres
        self.greek_bay = greek_bay      #Caractère unicode grec de la désignation de bayer
        self.pos = (self.x, self.y, self.z) = x, y, z #(cosdec*cosra, cosdec*sinra, sindec) #Position sur la sphere celeste unité
        self.total_feature_length = tfl        #Longueur totale des segments de double_triangle_feature
        self.double_triangle_feature = dtf     #Longueurs et angles des deux triangles formé par les 4 etoiles les plus proches

        self.imagematch = None                  #Objet Etoile_Image correspondant si trouvé
        self.gnomic_projection_map = None       #Etoiles proches et leurs coordonnees dans la base adaptée à cet étoile

    def calcxyz(self):
        cdec, sdec, cra, sra = cos(self.dec), sin(self.dec), cos(self.ra), sin(self.ra)
        self.pos = (self.x, self.y, self.z) = (cdec*cra, cdec*sra, sdec) #Position sur la sphere celeste unité

    def displayname(self):  #Nom à afficher si trouvée (bayer, flamsteed ou id)
        display_name = ""
        bay = self.greek_bay if config.GREEK_DISPLAY else self.bayer
        if config.DISPLAY_MODE == "bayerflam" and bay != None:
            display_name = bay + " " + self.con
        elif config.DISPLAY_MODE  == "bayerflam" and self.flam != None:
            display_name = self.flam + " " + self.con
        elif config.DISPLAY_MODE  == "hip" and self.hip != None:
            display_name = "HIP" + str(self.hip)
        elif config.DISPLAY_MODE  == "mixt":
            if self.proper != None:
                display_name = self.proper
            elif bay != None:
                display_name = bay + " " + self.con
            elif self.flam != None: 
                display_name = self.flam + " " + self.con
            elif self.hip != None:
                display_name = "HIP"+str(self.hip)
        else:
            display_name = str(self.id)
        return display_name
        
    def simbad(self): #Lien vers la page Simbad de l'etoile
        if self.hip != None: 
            return "http://cdsportal.u-strasbg.fr/?target=hip"+self.hip #sinon utiliser les coordonnees mais flm là
        else: return False

    def wikipedia(self): #Lien vers la page Wikipedia de l'etoile
        if (self.full_letter != None and self.genitive != None): #a ameliorer pour adapter aux Sig-1, Sig-2 etc.
            return "https://en.wikipedia.org/wiki/"+self.full_letter+"_"+self.genitive
        else: return False

    def load_gnomic(self):
        if self.gnomic_projection_map == None:
            path = config.GNOMIC_PATH + str(self.id) + ".csv"
            self.gnomic_projection_map = []
            with open(path, newline='') as file:
                reader = csv.reader(file, delimiter=',')
                for row in reader:
                    starid, x, y = int(row[0]), float(row[1]), float(row[2])
                    self.gnomic_projection_map.append((starid,(x,y)))

    def save_gnomic(self):
        path = "Algo_etoiles/treated_data/"+str(self.id)+".csv"
        with open(path, 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            for (star, (x,y)) in self.gnomic_projection_map:
                spamwriter.writerow([star.id, x, y])
            csvfile.close()

class Etoile_image: #Type enregistrement pour les étoiles de l'image

    def __init__(self, pos):
        self.pos = self.x, self.y = pos

        self.F1 = None                  #Première feature
        self.F1_length = None           #Longueur totale des segments de F1
        self.F2 = None                  #Deuxième feature
        self.F2_length = None           #Longueur totale des segments de F2
        self.normalized_map = None      #Etoiles proches et leurs coordonnees dans la base adaptée à cet étoile
        self.closest_star = None

        self.starmatch = None           #Objet Star correspondant si trouvé
