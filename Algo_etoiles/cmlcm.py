from PIL import Image
import numpy as np
import config
import types_perso

img_original = Image.open(config.IMAGE_PATH)
img_grayscale = img_original.convert('L')

img_size = img_width, img_height = img_original.size

W4 = (1/70)*np.array([[ 2, 2, 2, 2, 2],
                      [-1,-1,-1,-1,-1],
                      [-2,-2,-2,-2,-2],
                      [-1,-1,-1,-1,-1],
                      [ 2, 2, 2, 2, 2]])

W5 =(1/100)*np.array([[ 4, 2, 0,-2,-4],
                      [ 2, 1, 0,-1,-2],
                      [ 0, 0, 0, 0, 0],
                      [-2,-1, 0, 1, 2],
                      [-4,-2, 0, 2, 4],])

W6 = (1/70)*np.array([[ 2,-1,-2,-1, 2],
                      [ 2,-1,-2,-1, 2],
                      [ 2,-1,-2,-1, 2],
                      [ 2,-1,-2,-1, 2],
                      [ 2,-1,-2,-1, 2]])

mat0 = 2*W6
mat45 = W4+W5+W6
mat90 = 2*W4
mat135 = W4-W5+W6

mean = 1/9*np.array([[1, 1, 1],
                     [1, 1, 1],
                     [1, 1, 1]])

sigmat = 1/8*np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]])

def block(var, min, max):
    if var<min: return 0
    elif var>=max: return max-1
    return var

def convolve(image, width, height, mat, matsize):
    output = Image.new('L', (width, height))
    mid = matsize//2
    for x in range(width):
        for y in range(height):
            sum=0
            for i in range(matsize):
                for j in range(matsize):
                    xbis = block(x+i-mid, 0, width)
                    ybis = block(y+j-mid, 0, height)
                    sum += image.getpixel((xbis,ybis))*mat[i][j]
            color = int(block(sum, 0, 256))
            output.putpixel((x, y), color)
    return output

def convolve_imgtonp(image, width, height, mat, matsize):
    output = np.array([[0. for j in range(height)]for i in range(width)])
    mid = matsize//2
    for x in range(width):
        for y in range(height):
            sum=0
            for i in range(matsize):
                for j in range(matsize):
                    xbis = block(x+i-mid, 0, width)
                    ybis = block(y+j-mid, 0, height)
                    sum += image.getpixel((xbis,ybis))*mat[i][j]
            output[x][y] = sum
    return output

def convolve_nptonp(image, width, height, mat, matsize):
    output = np.array([[0. for j in range(height)] for i in range(width)])
    mid = matsize//2
    for x in range(width):
        for y in range(height):
            sum=0
            for i in range(matsize):
                for j in range(matsize):
                    xbis = block(x+i-mid, 0, width)
                    ybis = block(y+j-mid, 0, height)
                    sum += image[xbis][ybis]*mat[i][j]
            output[x][y] = sum
    return output

def nptoimg(input, width, height):
    output = Image.new('L', (width, height))
    for x in range(width):
        for y in range(height):
            output.putpixel((x, y), int(input[x][y]))
    return output
'''
result = convolve_imgtonp(img_grayscale, img_width, img_height, mat0, 5)
#result2 = convolve_nptonp(result, img_width, img_height, mean, 3)
img = nptoimg(result, img_width, img_height)
img.show()
'''
def calculate_SIG(image, width, height):
    f_m = convolve_imgtonp(img_grayscale, img_width, img_height, sigmat, 3)
    SIG = np.array([[0. for j in range(height)] for i in range(width)])
    for x in range(width):
        for y in range(height):
            SIG[x][y] = (f_m[x][y]/image.getpixel((x,y))) if (f_m[x][y]<image.getpixel((x,y))) else (1/10)
    return SIG

def calculate_cmlcm(image, width, height):
    P_0 = convolve_nptonp(convolve_imgtonp(img_grayscale, img_width, img_height, mat0, 5), img_width, img_height, mean, 3)
    P_45 = convolve_nptonp(convolve_imgtonp(img_grayscale, img_width, img_height, mat45, 5), img_width, img_height, mean, 3)
    P_90 = convolve_nptonp(convolve_imgtonp(img_grayscale, img_width, img_height, mat90, 5), img_width, img_height, mean, 3)
    P_135 = convolve_nptonp(convolve_imgtonp(img_grayscale, img_width, img_height, mat135, 5), img_width, img_height, mean, 3)
    SIG = calculate_SIG(image, width, height)
    En0 = np.array([[0. for j in range(height)] for i in range(width)])
    En45 = np.array([[0. for j in range(height)] for i in range(width)])
    En90 = np.array([[0. for j in range(height)] for i in range(width)])
    En135 = np.array([[0. for j in range(height)] for i in range(width)])
    for x in range(width):
        for y in range(height):
            P0 = P_0[x][y]
            P1 = P_0[x][y-5] if (y-5>=0) else 0
            P2 = P_0[x][y-3] if (y-3>=0) else 0
            P3 = P_0[x][y+3] if (y+3<height) else 0
            P4 = P_0[x][y+5] if (y+5<height) else 0
            if (P2>0 and P0<0 and P3>0):
                En0[x][y] = (P2-P0+P3)*SIG[x][y]
            elif (P1>0 and P0<0 and P4>0):
                En0[x][y] = (P1-P0+P4)*SIG[x][y]
            else:
                En0[x][y] = (P2-P0+P3)/100 if (P2-P0+P3>0) else 0

            P0 = P_45[x][y]
            P1 = P_45[x-5][y-5] if (y-5>=0 and x-5>=0) else 0
            P2 = P_45[x-3][y-3] if (y-3>=0 and x-3>=0) else 0
            P3 = P_45[x+3][y+3] if (y+3<height and x+3<width) else 0
            P4 = P_45[x+5][y+5] if (y+5<height and x+5<width) else 0
            if (P2>0 and P0<0 and P3>0):
                En45[x][y] = (P2-P0+P3)*SIG[x][y]
            elif (P1>0 and P0<0 and P4>0):
                En45[x][y] = (P1-P0+P4)*SIG[x][y]
            else:
                En45[x][y] = (P2-P0+P3)/100 if (P2-P0+P3>0) else 0
            
            P0 = P_90[x][y]
            P1 = P_90[x-5][y] if (x-5>=0) else 0
            P2 = P_90[x-3][y] if (x-3>=0) else 0
            P3 = P_90[x+3][y] if (x+3<width) else 0
            P4 = P_90[x+5][y] if (x+5<width) else 0
            if (P2>0 and P0<0 and P3>0):
                En90[x][y] = (P2-P0+P3)*SIG[x][y]
            elif (P1>0 and P0<0 and P4>0):
                En90[x][y] = (P1-P0+P4)*SIG[x][y]
            else:
                En90[x][y] = (P2-P0+P3)/100 if (P2-P0+P3>0) else 0
            
            P0 = P_135[x][y]
            P1 = P_135[x-5][y+5] if (y+5<height and x-5>=0) else 0
            P2 = P_135[x-3][y+3] if (y+3<height and x-3>=0) else 0
            P3 = P_135[x+3][y-3] if (y-3>=0 and x+3<width) else 0
            P4 = P_135[x+5][y-5] if (y-5>=0 and x+5<width) else 0
            if (P2>0 and P0<0 and P3>0):
                En135[x][y] = (P2-P0+P3)*SIG[x][y]
            elif (P1>0 and P0<0 and P4>0):
                En135[x][y] = (P1-P0+P4)*SIG[x][y]
            else:
                En135[x][y] = (P2-P0+P3)/100 if (P2-P0+P3>0) else 0
            
    En0 = En0/sum(En0)
    En45 = En45/sum(En45)
    En90 = En90/sum(En90)
    En135 = En135/sum(En135)
    f_t = En0*En45*En90*En135
    return f_t

def threshold_segmentation(cmlcm, width, height, regionsize, s, l):
    I = width//regionsize
    J = height//regionsize
    result = np.array([[0 for y in range(height)] for x in range(width)])
    for i in range(I):
        for j in range(J):
            mean = 0
            npixels = (min(width, (i+1)*regionsize)-(i*regionsize)) * (min(height, (j+1)*regionsize)-(j*regionsize)) 
            for x in range(i*regionsize, min(width, (i+1)*regionsize)):
                for y in range(j*regionsize, min(height, (j+1)*regionsize)):
                    mean += cmlcm[x][y]
            mean /= npixels
            var = 0
            for x in range(i*regionsize, min(width, (i+1)*regionsize)):
                for y in range(j*regionsize, min(height, (j+1)*regionsize)):
                    var += (cmlcm[x][y] - mean)**2
            var /= npixels
            std = np.sqrt(var)
            threshold = s*mean + l*std
            for x in range(i*regionsize, min(width, (i+1)*regionsize)):
                for y in range(j*regionsize, min(height, (j+1)*regionsize)):
                    result[x][y] = 1 if cmlcm[x][y] >= threshold else 0
    return result



cmlcm = calculate_cmlcm(img_grayscale, img_width, img_height)
result = threshold_segmentation(cmlcm, img_width, img_height, 256, 1, 1)
img = nptoimg(255*result, img_width, img_height)
img.show()


def pixel_to_NB(pixel): #permet de choisir le contraste de la conversion en noir et blanc
    return 1 if pixel > config.BLACK_WHITE_THRESHOLD else 0

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
