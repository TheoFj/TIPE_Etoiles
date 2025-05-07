from PIL import Image, ImageDraw
import config

#APERCUS

def affiche_etoiles(L, central, size, image_original):
    layer = Image.new("RGBA", size, (255, 255, 255, 0)) #cree une image transparente pour afficher le texte sur fond transparent avant de le combiner à image_original
    drawing_instance = ImageDraw.Draw(layer)
    image_copy = image_original.copy().convert('RGBA')
    r = 2
    for (x,y) in L:
        drawing_instance.ellipse((x-r,y-r,x+r,y+r), outline = (0, 255, 255, 255), width = 1)
    r=2*r
    (xc, yc) = central.pos
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


#CENTROIDES

def display_centroids(liste_etoiles_image, image_original, size):
    new_layer = Image.new("RGBA", size, (255, 255, 255, 0)) #crée une image transparente pour afficher le texte sur fond transparent avant de le combiner à image_original
    drawing_instance = ImageDraw.Draw(new_layer)
    image_copy = image_original.copy().convert('RGBA')

    for etoile in liste_etoiles_image: #affichage
        drawing_instance.rectangle((round(etoile.x)-1, round(etoile.y)-1, round(etoile.x)+1, round(etoile.y)+1), (0,0,255,255), (0,0,255,255), width=1)

    return Image.alpha_composite(image_copy, new_layer)


#RESULTATS

def draw_name_pillow(etoile, green_draw, red_draw):
    x, y = etoile.pos
    r = config.DISPLAY_CIRCLE_RADIUS
    
    if y-r-4-config.FONT_SIZE>0:
        textxy = (x, y-r-4)
        mode = "ms" #middle baseline
    else:
        textxy = (x, y+r+4)
        mode = "mt" #middle top
    circle_x0y0x1y1 = (x-r, y-r, x+r, y+r)

    if etoile.starmatch != None:
        green_draw.text(textxy, etoile.starmatch.displayname(), font = config.FONT, fill = 1, anchor = mode)
        green_draw.ellipse(circle_x0y0x1y1, outline = 1, width = 1)
    else:
        red_draw.ellipse(circle_x0y0x1y1, outline = (255,0,0,255), width = 1)

def affiche_resultat_pillow(liste_etoiles_image, image_original, size, constellations): #affiche le nom de l'etoile identifiee pour chaque etoile de l'image
    green_layer = Image.new("1", size, 0) #cree une image binaire pour que le texte soit pixellise et pas flou donc lisible
    green_draw = ImageDraw.Draw(green_layer)
    red_layer = Image.new("RGBA", size, (0,0,0,0))
    red_draw = ImageDraw.Draw(red_layer)
    image_copy = image_original.copy().convert('RGBA')
    
    for et1 in liste_etoiles_image:
        draw_name_pillow(et1, green_draw, red_draw)
        #affichage constellations:
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
                            green_draw.line((et1.x, et1.y, et2.x, et2.y), 1, width = 1)

    green_layer2 = Image.new("RGBA", size, (0,0,0,0))
    for x in range(size[0]):
        for y in range(size[1]):
            if green_layer.getpixel((x,y))==1: green_layer2.putpixel((x,y),(0,255,0,255))
    return Image.alpha_composite(image_copy, Image.alpha_composite(red_layer, green_layer2))
