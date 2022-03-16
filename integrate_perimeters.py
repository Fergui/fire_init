import sys
import pickle
import logging
import os.path as osp
from process_perimeter_masks import integrate_perims,add_smoke

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    if len(sys.argv) != 1:
        print('ERROR: {} pkl_path'.format(sys.argv[0]))
    elif not osp.exists(sys.argv[1]):
        print('ERROR: pkl_path {} not existent'.format(sys.argv[1]))
    pkl_path = sys.argv[1]
    wrfinput_path = 'wrfinput_d03'
    with open(pkl_path,'rb') as f:
        r = pickle.load(f)
    TIGN_G = r['TIGN_G']
    FUEL_MASK = r['FUEL_MASK']
    params = r['params']
    # integrate perimeter interpolation
    integrate_perims(wrfinput_path, TIGN_G, FUEL_MASK, **params)
    # add smoke from previous forecast
    add_smoke()
