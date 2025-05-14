from math import *
import csv
import geometry
import types_perso
import kdtree
import heapq

#LIRE https://www.mdpi.com/1424-8220/20/11/3027

"""
CONSTANTES
"""
DATA_BASE_NAME = "athyg_modified_vmagmax6.csv"
DATA_BASE_PATH = "Algo_etoiles/databasecsv/"+DATA_BASE_NAME
DATA_BASE_TREATED_PATH = "Algo_etoiles/databasecsv/treated_"+DATA_BASE_NAME

def intbis(s):
    return None if (s=='' or s=='#N/A') else int(s)
def strbis(s):
    return None if (s=='' or s=='#N/A') else str(s)
def fltbis(s):
    return None if (s=='' or s=='#N/A') else float(s)

csv_file = open(DATA_BASE_PATH, "r", encoding="utf-8")
next(csv_file)

reader = csv.reader(csv_file, delimiter=',')
DATA_BASE = [types_perso.Star(                      #construction de DATA_BASE: liste d'objets de type Star, conversions un peu inutiles mais copiees
                    id=int(row[0]),
                    hip=strbis(row[1]),
                    bayer=strbis(row[2]),
                    flam=strbis(row[3]),
                    con=str(row[4]),
                    proper=strbis(row[5]),
                    ra=float(row[6]),
                    dec=float(row[7]),
                    mag=float(row[8]),
                    full=strbis(row[9]),
                    gen=strbis(row[10]),
                    greek_bay=strbis(row[11]),
                    x=None, y=None, z=None, 
                    tfl=None, 
                    dtf=None, 
                    ) for row in reader]
csv_file.close()


for star in DATA_BASE:
    star.calcxyz()
TREE = kdtree.build_tree(DATA_BASE.copy(), dim=3, dir=0)

csv_file2 = open(DATA_BASE_PATH, "r", encoding="utf-8")
csv_save = open(DATA_BASE_TREATED_PATH, 'w', newline='', encoding="utf-8")
reader2 = csv.reader(csv_file2, delimiter=',')
writer = csv.writer(csv_save, delimiter=',')

row1 = next(reader2)
writer.writerow(row1+["x","y","z","tfl","th01","th021","th12","th023","th03","th23","l01","l02","l12","l02","l03","l23"])

for star in DATA_BASE:
    row = next(reader2)
    assert(star.id == int(row[0]))
    geometry.calcul_gnomic_dtf_tfl(star, DATA_BASE, TREE)
    writer.writerow(row+[star.x,star.y,star.z,star.total_feature_length]+star.double_triangle_feature)
    star.save_gnomic()

csv_save.close()
csv_file2.close()