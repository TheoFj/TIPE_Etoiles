from PIL import ImageFont

LAMBD = 50 #nombre maximum d'etoiles de "reference"
DISPLAY_MODE = "mixt" #bayerflam / hip / mixt
FONT_SIZE = 10
DISPLAY_CIRCLE_RADIUS = 3 #rayon des cercles lors de l'affichage des etoiles trouvees

IMAGE_PATH = "Algo_etoiles/images/orion2.jpg"
DATA_BASE_PATH = "Algo_etoiles/databasecsv/treated_athyg_modified_vmagmax6.csv"
CENTROIDS_SAVE_PATH = "resultats/centroids.png"
RESULTS_SAVE_PATH = "resultats/results7.png"
RESULTSPDF_SAVE_PATH = "resultatspdf/results7.pdf"
GNOMIC_PATH = "Algo_etoiles/treated_data/"
FONT_PATH = "Algo_etoiles/arial.ttf"
CON_DIC_PATH = "Algo_etoiles/constellations_graph.txt"
FONT = ImageFont.truetype(font = FONT_PATH, size = FONT_SIZE)
ID_THRESHOLD = 0.15


CONVOLUTION_OR_NOT = False
BLACK_WHITE_THRESHOLD = 130
