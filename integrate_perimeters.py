try:
    from .tools import load_pkl,integrate_init,add_smoke
except:
    from tools import load_pkl,integrate_init,add_smoke
import os.path as osp
import logging
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
if len(sys.argv) != 2:
    print('ERROR: {} pkl_path'.format(sys.argv[0]))
elif not osp.exists(sys.argv[1]):
    print('ERROR: pkl_path {} not existent'.format(sys.argv[1]))
pkl_path = sys.argv[1]
data = load_pkl(pkl_path)
# integrate perimeter interpolation and fuel mask
integrate_init('wrfinput_d03', data['TIGN_G'], data['FUEL_MASK'], data['params']['outside_time'])
# add smoke from previous forecast
add_smoke(
    ['wrfout_d01','wrfout_d02','wrfout_d03'],
    ['wrfinput_d01','wrfinput_d02','wrfinput_d03']
)
