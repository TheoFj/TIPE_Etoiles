
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
