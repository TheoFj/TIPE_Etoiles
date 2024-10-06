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

print(parse_constellation_file("Algo_etoiles/constellations_graph.txt"))