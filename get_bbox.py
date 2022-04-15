try:
    from .tools import read_bbox
except:
    from tools import read_bbox
import os.path as osp
import sys

if len(sys.argv) != 2:
    print('usage: python {} wrf_path'.format(sys.argv[0]))
    sys.exit(1)
elif not osp.exists(sys.argv[1]):
    print('error: wrf_path {} not existent'.format(sys.argv[1]))
wrfout_path = sys.argv[1]
bbox = read_bbox(wrfout_path)
print('{},{},{},{}'.format(*bbox))