from process_perimeter_masks import read_wrfout
import os.path as osp
import sys

if __name__ == '__main__':
    if len(sys.argv) != 1:
        print('ERROR: {} wrf_path'.format(sys.argv[0]))
    elif not osp.exists(sys.argv[1]):
        print('ERROR: wrf_path {} not existent'.format(sys.argv[1]))
    wrfout_path = sys.argv[1]
    ifm,ifn,ffm,ffn,fxlon,fxlat = read_wrfout(wrfout_path)
    bbox = (fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max())
    print('{},{},{},{}'.format(*bbox))
