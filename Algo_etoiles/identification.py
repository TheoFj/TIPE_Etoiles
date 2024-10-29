import config

def dtf_diff(dtf1, dtf2): #cf. eq(6) de https://www.mdpi.com/1424-8220/20/11/3027
    return sum([abs(dtf1[i] - dtf2[i]) for i in range(12)])

def identify(D0, M0):
    '''
        Une fois que des etoiles D0 et M0 semblent correspondre via la methode des doubles triangles,
        on compare leurs projections normalisees afin d'établir des correspondances,
        le résuultat final est l'identification avec le meilleur "r_score": le nombre d'étoiles identifiées.
    '''
    #ATTENTION L'ORDRE EST PLUS OU MOINS ALEATOIRE: LA PREMIRE ETOILE DE D0.gnomic_projection_map N'EST EN GENERAL PAS D0!!
    D0.load_gnomic() #besoin d'importer ces valeurs seulement pour les etoiles choisies par la "premiere selction"
    r_score = 0
    matchlist = []
    for (etoile_image, (x2,y2)) in M0.normalized_map:
        for (catalog_star_id, (x1,y1)) in D0.gnomic_projection_map:
            if abs(x1-x2)<=config.ID_THRESHOLD and abs(y1-y2)<=config.ID_THRESHOLD:
                r_score += 1
                matchlist.append((catalog_star_id, etoile_image))
                break #a ameliorer
    return r_score, matchlist, M0

def closest_dtf(dtf, data_base): #renvoie l'objet Star de data_base avec la feature la plus proche de dtf
    best = data_base[0]
    mindiff = dtf_diff(best.double_triangle_feature, dtf)
    for star in data_base:
        diff = dtf_diff(star.double_triangle_feature, dtf)
        if diff < mindiff:
            best, mindiff = star, diff
    return best
