
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

def affiche_etoiles(L):
    layer = Image.new("RGBA", SIZE, (255, 255, 255, 0)) #cree une image transparente pour afficher le texte sur fond transparent avant de le combiner à image_original
    drawing_instance = ImageDraw.Draw(layer)
    image_copy = IMAGE_ORIGINAL.copy().convert('RGBA')
    r = 2
    (a,b)=L[0]
    print(a,b)
    drawing_instance.ellipse((a-4,b-4,a+4,b+4), outline = (255, 0, 255, 255), width = 3)
    for (x,y) in L[1:]:
        drawing_instance.ellipse((x-r,y-r,x+r,y+r), outline = (0, 255, 255, 255), width = 1)
    return Image.alpha_composite(image_copy, layer)


def changement_normalise_vers_image(M0, M1, map):
    #Calcul des coordonnées sur l'image d'une liste de coordonnées normalsiées
    E1 = vectpp(M0.xy, M1.xy)
    E2 = (-E1[1], E1[0])

    L=[]
    for starandvect in map:
        x = E1[0]*starandvect[1][0] + E2[0]*starandvect[1][1] + M0.x
        y = E1[1]*starandvect[1][0] + E2[1]*starandvect[1][1] + M0.y
        L.append((x,y))
    return L

