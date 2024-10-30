from PIL import ImageFont

LAMBD = 50 #nombre maximum d'etoiles de "reference"
DISPLAY_MODE = "mixt" #bayerflam / hip / mixt
FONT_SIZE = 11
DISPLAY_CIRCLE_RADIUS = 3 #rayon des cercles lors de l'affichage des etoiles trouvees
GREEK_DISPLAY = False

DATA_BASE_TREATED_PATH = "Algo_etoiles/databasecsv/treated_athyg_modified_vmagmax6.csv"
GNOMIC_PATH = "Algo_etoiles/treated_data/"
FONT_PATH = "Algo_etoiles/arial.ttf"
CON_DIC_PATH = "Algo_etoiles/constellations_graph.txt"
FONT = ImageFont.truetype(font = FONT_PATH, size = FONT_SIZE)
ID_THRESHOLD = 0.15

IMAGE_PATH = "Algo_etoiles/images/Orion.jpg"
CONVOLUTION_OR_NOT = False
BLACK_WHITE_THRESHOLD = 125

SAVE_CENTROIDS = False
CENTROIDS_SAVE_PATH = "resultats/aaaaa.png"
SAVE_IMAGE = False
RESULTS_SAVE_PATH = "resultats/aaaaa.png"
SAVE_PDF = False
PDF_SAVE_PATH = "resultatspdf/aaaaa.pdf"