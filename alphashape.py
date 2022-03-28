import sys
import logging
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from geometry import alpha_shape, coords_to_polys, simplify_coords, merc_to_lonlat, polys_to_coords
from utils import save_pkl

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
if len(sys.argv) > 1:
    hotspots_path = sys.argv[1]
else:
    hotspots_path = 'arcgis_hotspots.pkl'
now = datetime.utcnow().replace(tzinfo=timezone.utc)
alpha = 375
r = pd.read_pickle(hotspots_path)
#r = r[r['confidence'] != 'low']
x = r['SHAPE'].map(lambda row: row['x'])
y = r['SHAPE'].map(lambda row: row['y'])
points = np.vstack([x,y]).T
if len(points) < 4:
    logging.warning('{} hotspots found, less than 4 so using previous hotspots...'.format(len(points)))
    sys.exit()
t = r['esritimeutc'].sort_values(ascending=False).iloc[0].to_pydatetime(warn=False).replace(tzinfo=timezone.utc) 
# Computing the alpha shape
perims = alpha_shape(points, alpha=alpha, only_outer=True)
if len(perims):
    # set parameters to reduce complexity of multipolygons
    params = {
        'min_inner_coords': 20,     # minimum number of coordinates for an inner polygon
        'max_inner_coords': 1000,   # maximum number of coordinates for an inner polygon
        'min_outer_coords': 4,      # minimum number of coordinates for an outer polygon
        'max_outer_coords': 10000   # maximum number of coordinates for an outer polygon
    }
    perims = simplify_coords(perims,**params)
    polys = coords_to_polys(perims)
    transf_polys = merc_to_lonlat(polys)
    transf_perims = polys_to_coords(transf_polys)
    # Save perimeters
    if len(transf_perims):
        save_pkl((t,transf_perims),'perim2.pkl')
else:
    transf_perims = perims
save_pkl((t,transf_perims),'hotspots_perims_{:02d}{:02d}_{:02d}z.pkl'.format(t.month,t.day,t.hour))
