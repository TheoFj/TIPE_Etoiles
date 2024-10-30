import csv

# Ouvrir le fichier CSV en mode lecture avec encodage UTF-8
with open("Algo_etoiles/greek_data.csv", "r", encoding="utf-8") as file:
    reader = csv.reader(file)
    
    # Lire et afficher les données
    for row in reader:
        print(row)

def unidecode(s):
    return None if (s=='' or s=="#N/A") else bytes(s).decode("utf-8")