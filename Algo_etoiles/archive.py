'''
PDF
'''

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

'''
KDTREES
'''

def tripartition(lst, pivot, dir):
    lows = [x for x in lst if x.pos[dir] < pivot.pos[dir]]
    pivots = [x for x in lst if x.pos[dir] == pivot.pos[dir]]
    highs = [x for x in lst if x.pos[dir] > pivot.pos[dir]]
    return lows, pivots, highs


def bipartition(lst, pivot, dir):
    lows = [x for x in lst if (x.pos[dir] <= pivot.pos[dir] and x.pos != pivot.pos)]
    highs = [x for x in lst if (x.pos[dir] <= pivot.pos[dir] and x.pos != pivot.pos)]
    return lows, highs


def quickselect(lst, k, dir): #Implementation de quickselect utilisant la methode de la mediane des medianes pour le choix du pivot
    if len(lst) <= 5:
        return sorted(lst, key = lambda star: star.pos[dir])[k]
    
    groups = [lst[i:i + 5] for i in range(0, len(lst), 5)]
    medians = [sorted(group, key = lambda star: star.pos[dir])[len(group) // 2] for group in groups]

    median_of_medians = quickselect(medians, len(medians) // 2, dir)

    lows, pivots, highs = tripartition(lst, median_of_medians, dir)

    if k < len(lows):
        return quickselect(lows, k, dir)
    elif k < len(lows) + len(pivots):
        return median_of_medians
    else:
        return quickselect(highs, k - len(lows) - len(pivots), dir)


def build_tree1(lst, dim, dir): #O(nlog(n))
    if lst==[]:
        return None
    
    medianstar = quickselect(lst, len(lst)//2, dir)
    lows, highs = bipartition(lst, medianstar, dir)

    return KDtree(medianstar, build_tree1(lows, dim, (dir+1)%dim), build_tree1(highs, dim, (dir+1)%dim))

'''
GEOMETRY
'''

def closest_nstars(starandvect_list, n, with_zero, dim):
    '''
    Renvoie les n vecteurs les plus proches de 0, en dimension 2 ou 3, en incluant le 0 ou non.
    O(len(starandvect_list)*n) ATTENTION: compare par rapport à (0,0)/(0,0,0), pour comparer à un autre point, mettre en entree les vecteurs "P0->P1"
    '''
    if dim==2: L_dist = [(norm2_sqr(starandvect[1]), starandvect) for starandvect in starandvect_list]
    elif dim==3: L_dist = [(norm3_sqr(starandvect[1]), starandvect) for starandvect in starandvect_list]
    else: print("????????????")
    best_n = []

    for (dist, starandvect) in L_dist:
        i=0
        while i<(len(best_n)):
            if dist <= best_n[i][0]:
                break
            i+=1
        best_n.insert(i, (dist, starandvect))
        best_n = best_n[:(n+1)]

    if with_zero:
        return [best_n[i][1] for i in range(0,n)] #on ne renvoie pas le dernier, indicé n, car on ne renvoie que les n meilleures
    else:
        return [best_n[i][1] for i in range(1,n+1)] #on ne renvoie pas le numero 0 car il s'agit de l'etoile "pos"
def calcul_map_dtf_tfl_2darchive(M0, image_star_list, tree):
    '''
    Calculs à faire sur chaque étoile (ou une partie des étoiles) de l'image
    '''

    starandvect_list = [(M, vectpp(M0.pos, M.pos)) for M in image_star_list]
    best_4 = closest_nstars(starandvect_list, 4, False, 2)
    M0.closest_star = best_4[0][0]
    E1 = best_4[0][1]
    E2 = (-E1[1], E1[0])
    k = norm2_sqr(E1)

    #PROJECTIONS DANS LA NOUVELLE BASE
    L = [(starandvect[0], (dot_prod2(starandvect[1], E1)/k, dot_prod2(starandvect[1], E2)/k)) for starandvect in starandvect_list]

    M0xy = (0., 0.)       
    M1xy = (1., 0.)
    M2xy = (dot_prod2(best_4[1][1], E1)/k, dot_prod2(best_4[1][1], E2)/k)
    M3xy = (dot_prod2(best_4[2][1], E1)/k, dot_prod2(best_4[2][1], E2)/k)
    M4xy = (dot_prod2(best_4[3][1], E1)/k, dot_prod2(best_4[3][1], E2)/k)

    M0.F1 = F1 = calcul_double_triangle_feature(M0xy, M1xy, M2xy, M3xy)
    M0.F2 = F2 = calcul_double_triangle_feature(M0xy, M2xy, M3xy, M4xy)
    M0.F1_length = calcul_total_feature_length(F1)
    M0.F2_length = calcul_total_feature_length(F2)
    M0.normalized_map = L
    return
