import sys
import os.path as osp
from utils import load_pkl
from plot_tools import plot_perims,plot_coords,plot_tign,plot_fmask

if len(sys.argv) != 2:
    print('ERROR: {} pkl_path'.format(sys.argv[0]))
elif not osp.exists(sys.argv[1]):
    print('ERROR: pkl_path {} not existent'.format(sys.argv[1]))
pkl_path = sys.argv[1]
# plotting fire arrival time result
plot_tign(pkl_path)
try:
    plot_perims(pkl_path,show=True)
except:
    p = load_pkl(pkl_path)
    plot_coords(p['PERIMS']['perim1'],color='g')
    plot_coords(p['PERIMS']['perim2'],color='k',show=True)
# plotting fuel mask result
plot_fmask(pkl_path)
try:
    plot_perims(pkl_path,show=True)
except:
    p = load_pkl(pkl_path)
    plot_coords(p['PERIMS']['perim1'],color='g')
    plot_coords(p['PERIMS']['perim2'],color='k',show=True)