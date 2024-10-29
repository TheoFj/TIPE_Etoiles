#from fpdf import FPDF
from PIL import Image, ImageDraw, ImageFont
import config

def draw_name_pillow(etoile, drawing_instance):
    x, y = etoile.xy
    r = config.DISPLAY_CIRCLE_RADIUS
    
    if y-r-4-config.FONT_SIZE>0:
        textxy = (x, y-r-4)
        mode = "ms" #middle baseline
    else:
        textxy = (x, y+r+4)
        mode = "mt" #middle top
    circle_x0y0x1y1 = (x-r, y-r, x+r, y+r)

    if etoile.starmatch != None:
        #print(self.xy, self.starmatch)
        drawing_instance.text(textxy, etoile.starmatch.display_name, font = config.FONT, fill = (0, 255, 0, 255), anchor = mode)
        drawing_instance.ellipse(circle_x0y0x1y1, outline = (0, 255, 0, 255), width = 1)
    else:
        #print(self.xy, None)
        drawing_instance.ellipse(circle_x0y0x1y1, outline = (255, 0, 0, 255), width = 1)
        

def draw_name_pdf(etoile, pdf):
    x, y = etoile.xy
    r = config.DISPLAY_CIRCLE_RADIUS
    textx, texty = x-2.5*config.FONT_SIZE, y-r-4-config.FONT_SIZE
    if texty<0:
        texty = y+r+4
    circlex, circley = x-r, y-r

    if etoile.starmatch != None:
        pdf.set_xy(textx, texty)
        pdf.set_text_color(0, 255, 0)
        pdf.set_draw_color(0, 255, 0)
        pdf.write(5, etoile.starmatch.display_name, etoile.starmatch.simbad)
    else:
        pdf.set_draw_color(255, 0, 0)
    pdf.ellipse(circlex, circley, 2*r, 2*r)

def affiche_resultat_pdf(liste_etoiles_image, width, height): #affiche le nom de l'etoile identifiee pour chaque etoile de l'image
    pdf = FPDF("P", "pt", (width, height))
    pdf.add_font("police", "", config.FONT_PATH, uni=True)
    pdf.set_font("police", "", config.FONT_SIZE)
    pdf.set_margins(0, 0, 0)
    pdf.set_auto_page_break(False)
    pdf.add_page()
    pdf.image(config.IMAGE_PATH, 0, 0, width, height)
    for etoile in liste_etoiles_image:
        draw_name_pdf(etoile, pdf)
    
    return pdf

def affiche_etoiles(L, central, size, image_original):
    layer = Image.new("RGBA", size, (255, 255, 255, 0)) #cree une image transparente pour afficher le texte sur fond transparent avant de le combiner à image_original
    drawing_instance = ImageDraw.Draw(layer)
    image_copy = image_original.copy().convert('RGBA')
    r = 2
    for (x,y) in L:
        drawing_instance.ellipse((x-r,y-r,x+r,y+r), outline = (0, 255, 255, 255), width = 1)
    r=2*r
    (xc, yc) = central.xy
    drawing_instance.ellipse((xc-r,yc-r,xc+r,yc+r), outline = (255, 0, 255, 255), width = 2)
    return Image.alpha_composite(image_copy, layer)


def drawmap(map):
    #affiche un aperçu des etoiles dans le repere normalisé autour d'une etoile
    img = Image.new('1',(300,300),0)
    for starandvect in map:
        x,y = round(starandvect[1][0]*30+150), round(starandvect[1][1]*30+150)
        if 0<=x<300 and 0<=y<300:
            img.putpixel((x,y), 1)
    img.show()

def display_centroids(liste_etoiles_image, image_original, size):
    new_layer = Image.new("RGBA", size, (255, 255, 255, 0)) #crée une image transparente pour afficher le texte sur fond transparent avant de le combiner à image_original
    drawing_instance = ImageDraw.Draw(new_layer)
    image_copy = image_original.copy().convert('RGBA')

    for etoile in liste_etoiles_image: #affichage
        drawing_instance.rectangle((round(etoile.x)-1, round(etoile.y)-1, round(etoile.x)+1, round(etoile.y)+1), (0,0,255,255), (0,0,255,255), width=1)

    centroids = Image.alpha_composite(image_copy, new_layer)
    centroids.show("Centroïdes des étoiles repérées")
    centroids.save(config.CENTROIDS_SAVE_PATH)


def affiche_resultat_pillow(liste_etoiles_image, image_original, size, constellations): #affiche le nom de l'etoile identifiee pour chaque etoile de l'image
    text_layer = Image.new("RGBA", size, (255, 255, 255, 0)) #cree une image transparente pour afficher le texte sur fond transparent avant de le combiner à image_original
    drawing_instance = ImageDraw.Draw(text_layer)
    image_copy = image_original.copy().convert('RGBA')
    
    for et1 in liste_etoiles_image:
        draw_name_pillow(et1, drawing_instance)
        if et1.starmatch != None:
            con = et1.starmatch.con
            hip1 = et1.starmatch.hip
            n = constellations[con][0]
            l = constellations[con][1]
            for et2 in liste_etoiles_image:
                if et2.starmatch != None:
                    hip2 = et2.starmatch.hip
                    for i in range(n):
                        if l[i] == (hip1,hip2):
                            drawing_instance.line((et1.x, et1.y, et2.x, et2.y), (0, 255, 0, 255), width = 1)

    return Image.alpha_composite(image_copy, text_layer)
