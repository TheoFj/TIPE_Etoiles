import csv
import config
import types_perso
from random import randint

def intbis(s):
    return None if (s=='' or s=='#N/A') else int(s)
def strbis(s):
    return None if (s=='' or s=='#N/A') else str(s)
def fltbis(s):
    return None if (s=='' or s=='#N/A') else float(s)

def parse_constellation_line(line):
    links = []
    content = line.split(" ")
    if len(content) < 3:
        print("warning: malformated line when reading constellation")
        return None
    constellation_name = content[0]

    link_count = int(content[1])
    if len(content) < 3 + link_count*2:
        print("warning: malformated line when reading constellation")
        return None
    for i in range(3, link_count*2+3, 2):
        links.append((int(content[i]), int(content[i+1])))
    return constellation_name, link_count, links

def parse_constellation_file(filename):
    data = {}
    with open(filename, "r") as f:
        for line in f.readlines():
            x = parse_constellation_line(line)
            if x == None: continue
            name, n, links = x
            data[name] = (n, links)
    return data

def parse_database_file(filename):
    csv_data_base = open(filename, "r", encoding="utf-8")
    next(csv_data_base) #skip la premiere ligne
    reader = csv.reader(csv_data_base, delimiter=',')

    data_base = [types_perso.Star(                      #construction de DATA_BASE: liste d'objets de type Star
                    id=int(row[0]),
                    hip=intbis(row[1]),
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
                    x=float(row[12]),
                    y=float(row[13]),
                    z=float(row[14]),
                    tfl=str(row[15]),
                    dtf=[float(row[i]) for i in range(16,28)],
                    ) for row in reader]

    csv_data_base.close()
    return data_base

def get_by_attribute(data_base, attribute, val):
    if attribute == "bayer":
        return next((star for star in data_base if star.bayer == val), [None])
    elif attribute == "id":
        return next((star for star in data_base if star.id == val), [None])
    elif attribute == "hip":
        return next((star for star in data_base if star.hip == val), [None])


def choose_random(liste_etoiles, lambd): #a coder: si trop d'etoiles sur l'image, on en choisit lambda pour reduire le temps
    pass
