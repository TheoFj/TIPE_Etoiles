from PIL import ImageFont

LAMBD = 50 #nombre maximum d'etoiles de "reference"
DISPLAY_MODE = "mixt" #bayerflam / hip / mixt
FONT_SIZE = 11
DISPLAY_CIRCLE_RADIUS = 3 #rayon des cercles lors de l'affichage des etoiles trouvees
GREEK_DISPLAY = False

DATA_BASE_TREATED_PATH = "Algo_etoiles/databasecsv/treated_athyg_modified_vmagmax6.csv"
GNOMIC_PATH = "Algo_etoiles/star_maps/"
FONT_PATH = "Algo_etoiles/arial.ttf"
CON_DIC_PATH = "Algo_etoiles/constellations_graph.txt"
FONT = ImageFont.truetype(font = FONT_PATH, size = FONT_SIZE)

BLACK_WHITE_THRESHOLD = 110
CMLCM_OR_NOT = True
BLOCKSIZE = 256
S = 1
L = 2

ID_THRESHOLD = 0.15
ID_THRESHOLD2 = ID_THRESHOLD/10
N_ITE = 1000

IMAGE_PATH = "Algo_etoiles/allthesky/ursamajor28v-b.jpg"

SAVE_IMAGE = True
RESULTS_SAVE_PATH = "resultats/ursamajorprime3.png"
SAVE_CENTROIDS = False
CENTROIDS_SAVE_PATH = "resultats/aaaaa.png"
SAVE_PDF = False
PDF_SAVE_PATH = "resultatspdf/aaaaa.pdf"