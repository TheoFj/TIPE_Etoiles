from PIL import Image
import config
import types_perso
import geometry

IMAGE_ORIGINAL = Image.open(config.IMAGE_PATH)
IMAGE_GRAYSCALE = IMAGE_ORIGINAL.convert('L')

SIZE = WIDTH, HEIGHT = IMAGE_ORIGINAL.size

def block(var, min, max):
    if var<min: return 0
    elif var>=max: return max-1
    return var

def convolve(image, width, height, mat):
    output = Image.new('L', (width, height))
    for x in range(width):
        for y in range(height):
            acc=0
            for i in range(3):
                for j in range(3):
                    xbis = block(x+i-1, 0, width)
                    ybis = block(y+j-1, 0, height)
                    acc += image.getpixel((xbis,ybis))*mat[i][j]
            acc = block(acc, 0, 256)
            output.putpixel((x, y), acc)
    return output

if config.CONVOLUTION_OR_NOT == True:
    MAT = [[ 0,-1, 0],
           [-1, 4,-1],
           [ 0,-1, 0]]
    IMAGE_TREATED = convolve(IMAGE_GRAYSCALE, WIDTH, HEIGHT, MAT)
    IMAGE_TREATED.show()
else:
    IMAGE_TREATED = IMAGE_GRAYSCALE

def pixel_to_NB(pixel): #permet de choisir le contraste de la conversion en noir et blanc
    return 1 if pixel > config.BLACK_WHITE_THRESHOLD else 0

IMAGE_NB = IMAGE_TREATED.point(pixel_to_NB, mode='1')

IMAGE_NB.show()

def average_pix(pix_list): #renvoie le centroide d'une liste de pixels pris par etoile
    return (sum([pix[0] for pix in pix_list])/len(pix_list), sum([pix[1] for pix in pix_list])/len(pix_list))

def pix_list_of_star(image, pos, pixel_list): #détecte tous les pixels blancs "collés" à pos, les met en noir, et renvoie la liste de ces pixels
    pixel_list.append(pos)
    image.putpixel(pos, 0)
    x,y = pos
    width,height = image.size
    for i in [+1,0,-1]:
        for j in [+1,0,-1]:
            if 0<=x+i<width and 0<=y+j<height and image.getpixel((x+i,y+j)) == 1: #le pixel central a déjà été colorié en noir
                pix_list_of_star(image, (x+i,y+j), pixel_list)
    return pixel_list

def new_star(image, pos, star_list):
    star_list.append(types_perso.Etoile_image(average_pix(pix_list_of_star(image, pos, [])))) #ajoute à star_list le centre de la liste des pixels d'une étoile détectée en pos

image_temp = IMAGE_NB.copy() #utile car la fonction "new_star" a besoin de modifier l'image
LISTE_ETOILES_IMAGE = []
for y in range(HEIGHT): #remplit la liste des coordonnéees des étoiles sur l'image
    for x in range(WIDTH):
        if image_temp.getpixel((x,y)) == 1:
            new_star(image_temp, (x,y), LISTE_ETOILES_IMAGE)

