from datetime import datetime,timezone
import numpy as np
import logging
import pickle
import sys
from geometry import alpha_shape, coords_to_polys, simplify_coords, transform_polys, merc_to_lonlat

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    hotspots_path = 'arcgis_hotspots.pkl'
    now = datetime.utcnow().replace(tzinfo=timezone.utc)
    alpha = 375
    with open(hotspots_path,'rb') as f:
        r = pickle.load(f)
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
    # set parameters to reduce complexity of multipolygons
    params = {
        'min_inner_coords': 20,     # minimum number of coordinates for an inner polygon
        'max_inner_coords': 1000,   # maximum number of coordinates for an inner polygon
        'min_outer_coords': 4,      # minimum number of coordinates for an outer polygon
        'max_outer_coords': 10000   # maximum number of coordinates for an outer polygon
    }
    polys = coords_to_polys(perims)
    # Transform perimeters into WGS84
    transf_polys = transform_polys(polys, merc_to_lonlat, **params)
    # Save perimeters
    if len(transf_polys):
        with open('perim2.pkl','wb') as f:
            pickle.dump((t,transf_polys),f,protocol=3)
    with open('hotspots_perims_{:02d}{:02d}_{:02d}z.pkl'.format(now.month,now.day,now.hour),'wb') as f:
        pickle.dump((t,transf_polys),f,protocol=3)
