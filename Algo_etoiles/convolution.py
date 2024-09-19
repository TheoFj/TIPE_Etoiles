
from math import *
from PIL import Image, ImageDraw, ImageFont
import numpy as np

IMAGE_PATH = "Algo_etoiles/images/test.png"
SAVE_PATH = "resultats/abcde.png"

IMAGE_ORIGINAL = Image.open(IMAGE_PATH)
IMAGE_GRAYSCALE = IMAGE_ORIGINAL.convert('L')

SIZE = WIDTH, HEIGHT = IMAGE_ORIGINAL.size

MAT = [[ 0,-1, 0],
       [-1, 4,-1],
       [ 0,-1, 0]]

MAT2 = [[ 0, 0,-1, 0, 0],
        [ 0,-2,-3,-2, 0],
        [-1,-3,24,-3,-1],
        [ 0,-2,-3,-2, 0],
        [ 0, 0,-1, 0, 0]]

def scalmult(k, mat):
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            mat[i][j] *= k

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

RESULT1 = convolve(IMAGE_GRAYSCALE, WIDTH, HEIGHT, MAT)

BLACK_WHITE_THRESHOLD = 200
def pixel_to_NB(pixel): #permet de choisir le contraste de la conversion en noir et blanc
    return 1 if pixel > BLACK_WHITE_THRESHOLD else 0
RESULT2 = RESULT1.point(pixel_to_NB, mode='1')
RESULT3 = IMAGE_GRAYSCALE.point(pixel_to_NB, mode='1')


IMAGE_GRAYSCALE.show()
RESULT1.show()
RESULT2.show()
RESULT3.show()

