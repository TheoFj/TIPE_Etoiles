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
