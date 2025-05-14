from PIL import Image
import numpy as np
import types_perso

def average_pix(pix_list): #renvoie le centroide d'une liste de pixels pris par etoile
    return (sum([pix[0] for pix in pix_list])/len(pix_list), sum([pix[1] for pix in pix_list])/len(pix_list))

def pix_list_of_star(image, pos, pixel_list): #détecte tous les pixels blancs "collés" à pos, les met en noir, et renvoie la liste de ces pixels
    pixel_list.append(pos)
    image.putpixel(pos, 0)
    x,y = pos
    width,height = image.size
    for i in [+1,0,-1]:
        for j in [+1,0,-1]:
            if 0<=x+i<width and 0<=y+j<height and image.getpixel((x+i,y+j)) == 1:
                pix_list_of_star(image, (x+i,y+j), pixel_list)
    return pixel_list

def new_star(image, pos, star_list):
    star_list.append(types_perso.Etoile_image(average_pix(pix_list_of_star(image, pos, [])))) #ajoute à star_list le centre de la liste des pixels d'une étoile détectée en pos

def star_list(imgnb, width, height):
    imgtemp = imgnb.copy() #utile car la fonction "new_star" a besoin de modifier l'image
    liste_etoiles_image = []
    for y in range(height): #remplit la liste des coordonnéees des étoiles sur l'image
        for x in range(width):
            if imgtemp.getpixel((x,y)) == 1:
                new_star(imgtemp, (x,y), liste_etoiles_image)
    return liste_etoiles_image

# https://www.mdpi.com/2304-6732/9/1/13 CMLCM

W4 = (1/70)*np.array([[ 2,-1,-2,-1, 2], #transpose en apparence par rapport à l'article car on note les pixels des images (x,y) et x<->j et y<->i
                      [ 2,-1,-2,-1, 2],
                      [ 2,-1,-2,-1, 2],
                      [ 2,-1,-2,-1, 2],
                      [ 2,-1,-2,-1, 2]])

W5 =(1/100)*np.array([[ 4, 2, 0,-2,-4],
                      [ 2, 1, 0,-1,-2],
                      [ 0, 0, 0, 0, 0],
                      [-2,-1, 0, 1, 2],
                      [-4,-2, 0, 2, 4],])

W6 = (1/70)*np.array([[ 2, 2, 2, 2, 2],
                      [-1,-1,-1,-1,-1],
                      [-2,-2,-2,-2,-2],
                      [-1,-1,-1,-1,-1],
                      [ 2, 2, 2, 2, 2]])

MEAN = 1/9*np.array([[1, 1, 1],
                     [1, 1, 1],
                     [1, 1, 1]])

SIGMAT = 1/8*np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]])

def block(var, min, max):
    if var<min: return 0
    elif var>=max: return max-1
    return var

def is_in(x,y,width,height):
    return (0<=x<width and 0<=y<height)

def convolve_np(image, width, height, mat, matsize):
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

def imgtonp(input, width, height):
    return np.array([[input.getpixel((x,y)) for y in range(height)] for x in range(width)])

def nptoimg(input, width, height):
    output = Image.new('L', (width, height))
    for x in range(width):
        for y in range(height):
            output.putpixel((x, y), int(input[x][y]))
    return output

def nptoimgnb(input, width, height):
    output = Image.new('1', (width, height))
    for x in range(width):
        for y in range(height):
            output.putpixel((x, y), int(input[x][y]))
    return output

def calculate_SIG(image, width, height):
    f_m = convolve_np(image, width, height, SIGMAT, 3)
    sig = np.array([[0. for j in range(height)] for i in range(width)])
    for x in range(width):
        for y in range(height):
            sig[x][y] = (f_m[x][y]/image[x][y]) if (f_m[x][y]<image[x][y]) else (1/10)
    return sig

def calculate_En(P, sig, width, height, dx, dy):
    En = np.array([[0. for y in range(height)] for x in range(width)])
    for x in range(width):
        for y in range(height):
            P0 = P[x][y]
            P1 = P[x-dx*5][y-dy*5] if is_in(x-dx*5,y-dy*5,width,height) else 0
            P2 = P[x-dx*3][y-dy*3] if is_in(x-dx*3,y-dy*3,width,height) else 0
            P3 = P[x+dx*3][y+dy*3] if is_in(x+dx*3,y+dy*3,width,height) else 0
            P4 = P[x+dx*5][y+dy*5] if is_in(x+dx*5,y+dy*5,width,height) else 0
            if (P2>0 and P0<0 and P3>0):
                En[x][y] = (P2-P0+P3)*sig[x][y]
            elif (P1>0 and P0<0 and P4>0):
                En[x][y] = (P1-P0+P4)*sig[x][y]
            else:
                En[x][y] = (P2-P0+P3)/100 if (P2-P0+P3>0) else 0
    return En

def calculate_cmlcm(image, width, height): #L'image importee doit être en grayscale numpy
    Q4 = convolve_np(convolve_np(image, width, height, W4, 5), width, height, MEAN, 3) #On met le flou 3x3 dès maintenant car il y en a besoin après
    Q5 = convolve_np(convolve_np(image, width, height, W5, 5), width, height, MEAN, 3)
    Q6 = convolve_np(convolve_np(image, width, height, W6, 5), width, height, MEAN, 3)
    sig = calculate_SIG(image, width, height)
    P_0 = 2*Q6              #formule (7)
    P_45 = Q4+Q5+Q6
    P_90 = 2*Q4
    P_135 = Q4-Q5+Q6
    En0 = calculate_En(P_0, sig, width, height, 1, 0)
    En45 = calculate_En(P_45, sig, width, height, 1, 1)
    En90 = calculate_En(P_90, sig, width, height, 0, 1)
    En135 = calculate_En(P_135, sig, width, height, -1, 1)
    En0 = En0/sum(En0)
    En45 = En45/sum(En45)
    En90 = En90/sum(En90)
    En135 = En135/sum(En135)
    f_t = En0*En45*En90*En135
    return f_t

def threshold_segmentation(f_t, width, height, blocksize, s, l):
    result = np.array([[0 for y in range(height)] for x in range(width)])
    for X in range(0, width, blocksize):
        for Y in range(0, height, blocksize):
            block = np.array([col[Y:Y+blocksize] for col in f_t[X:X+blocksize]])
            threshold = s*np.mean(block) + l*np.std(block)
            for x in range(X, min(width, X+blocksize)):
                for y in range(Y, min(height, Y+blocksize)):
                    result[x][y] = 1 if f_t[x][y] >= threshold else 0
    return result

def cmlcmtotal(imggrayscale, width, height, blocksize, s, l):
    f_t = calculate_cmlcm(imgtonp(imggrayscale, width, height), width, height)
    result = threshold_segmentation(f_t, width, height, blocksize, s, l)
    return nptoimgnb(result, width, height)