import config

def dtf_diff(dtf1, dtf2): #cf. eq(6) de https://www.mdpi.com/1424-8220/20/11/3027
    return sum([abs(dtf1[i] - dtf2[i]) for i in range(12)])

def match_maps(mapimg, mapbdd, threshold):
    '''
        Une fois que des etoiles semblent correspondre via la methode des doubles triangles,
        on compare leurs projections normalisees afin d'établir des correspondances, le r_score est le nombre d'étoiles identifiées.
    '''
    r_score = 0
    matchlist = []
    for (etoile_image, (x2,y2)) in mapimg:
        for (catalog_star_id, (x1,y1)) in mapbdd:
            if abs(x1-x2) <= threshold and abs(y1-y2) <= threshold:
                r_score += 1
                matchlist.append((catalog_star_id, etoile_image))
                break #a ameliorer
    return r_score, matchlist

def closest_dtf(dtf, data_base): #renvoie l'objet Star de data_base avec la feature la plus proche de dtf
    best = data_base[0]
    mindiff = dtf_diff(best.double_triangle_feature, dtf)
    for star in data_base:
        diff = dtf_diff(star.double_triangle_feature, dtf)
        if diff < mindiff:
            best, mindiff = star, diff
    return best
