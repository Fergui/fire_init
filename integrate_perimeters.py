import sys
import logging
import os.path as osp
from utils import load_pkl,integrate_perims,add_smoke

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
if len(sys.argv) != 2:
    print('ERROR: {} pkl_path'.format(sys.argv[0]))
elif not osp.exists(sys.argv[1]):
    print('ERROR: pkl_path {} not existent'.format(sys.argv[1]))
pkl_path = sys.argv[1]
wrfinput_path = 'wrfinput_d03'
data = load_pkl(pkl_path)
TIGN_G = data['TIGN_G']
FUEL_MASK = data['FUEL_MASK']
params = data['params']
# integrate perimeter interpolation
integrate_perims(wrfinput_path, TIGN_G, FUEL_MASK, **params)
# add smoke from previous forecast
add_smoke()