import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
from utils import read_bbox, load_pkl

_colors = get_cmap('tab10')

def plot_bbox(bbox, show=False):
    """
    Plot bounding box.
    :param bbox: bounding box coordinates [min_lon,max_lon,min_lat,max_lat]
    :param show: optional, show the plot
    """
    xbox = [bbox[0],bbox[0],bbox[1],bbox[1],bbox[0]]
    ybox = [bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]]
    plt.plot(xbox,ybox,'k-')
    if show:
        plt.show()
    
def set_bbox(bbox, show=False):
    """
    Set bounding box in an existing plot.
    :param bbox: bounding box coordinates [min_lon,max_lon,min_lat,max_lat]
    :param show: optional, show the plot
    """
    plt.xlim([bbox[0],bbox[1]])
    plt.ylim([bbox[2],bbox[3]])
    if show:
        plt.show()

def plot_coords(coords, show=False, **kargs):
    """
    Plot polygons in list of coordinates format.
    :param coords: list of lists with outer and inner coordinates
    :param show: optional, show the plot
    :param kargs: other paramaters passed to matplotlib.pyplot.plot
    """
    gkargs = kargs.copy()
    for p in coords:
        if len(p):
            pp = plt.plot(np.array(p[0])[:,0],np.array(p[0])[:,1],**gkargs)
            kargs['color'] = pp[-1].get_color()
            if len(p) > 1:
                for ring in p[1:]:
                    plt.plot(np.array(ring)[:,0],np.array(ring)[:,1],**kargs)
    if show:
        plt.show()

def plot_perim(perim, show=False, **kargs):
    """
    Plot polygons in Perimeter format.
    :param perim: polygon in Perimeter format
    :param show: optional, show the plot
    :param kargs: other paramaters passed to matplotlib.pyplot.plot
    """
    gkargs = kargs.copy()
    for p in perim.coords:
        if len(p):
            pp = plt.plot(np.array(p[0])[:,0],np.array(p[0])[:,1],**gkargs)
            kargs['color'] = pp[-1].get_color()
            if len(p) > 1:
                for ring in p[1:]:
                    plt.plot(np.array(ring)[:,0],np.array(ring)[:,1],**kargs)
    if show:
        plt.show()

def plot_detections(path, fmt='r.', show=False, **kargs):
    """
    Plot fire detections stored in a pickle file.
    :param path: path to a pickle file with fire detections
    :param fmt: optional, format of markers
    :param show: optional, show the plot
    :param kargs: other paramaters passed to matplotlib.pyplot.plot
    """
    if 'ms' not in kargs:
        kargs['ms'] = 2
    data = load_pkl(path)
    plt.plot(data['longitude'], data['latitude'], fmt, **kargs)
    if show:
        plt.show()

def plot_perims(path, show=False):
    """
    Plot perimeters from a pickle file.
    :param path: path to a pickle file with perimeters
    :param show: optional, show the plot
    """
    data = load_pkl(path)
    fxlon = data['fxlon']
    fxlat = data['fxlat']
    perims = data['PERIMS']
    perims['perim1'].plot(color='g')
    perims['perim2'].plot(color='k')
    bbox = [fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max()]
    set_bbox(bbox,show=True)
    if show:
        plt.show()
    
def plot_tign(path, show=False):
    """
    Plot fire arrival time from a pickle file.
    :param path: path to a pickle file with fire arrival time
    :param show: optional, show the plot
    """
    data = load_pkl(path)
    params = data['params']
    fxlon = data['fxlon']
    fxlat = data['fxlat']
    TIGN_G = data['TIGN_G']
    TIGN_G[TIGN_G==params['outside_time']] = np.nan
    TIGN_G[TIGN_G==params['time_step']] = np.nan
    plt.pcolormesh(fxlon,fxlat,TIGN_G)
    bbox = [fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max()]
    set_bbox(bbox)
    if show:
        plt.show()

def plot_fmask(path, show=False):
    """
    Plot fuel mask from a pickle file.
    :param path: path to a pickle file with fuel mask
    :param show: optional, show the plot
    """
    data = load_pkl(path)
    fxlon = data['fxlon']
    fxlat = data['fxlat']
    FUEL_MASK = data['FUEL_MASK'].astype(float)
    FUEL_MASK[FUEL_MASK == 0] = np.nan
    plt.pcolormesh(fxlon,fxlat,FUEL_MASK,cmap='Greys',vmin=0,vmax=2)
    bbox = [fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max()]
    set_bbox(bbox)
    if show:
        plt.show()
    
def plot_pmask(fxlon, fxlat, mask, show=False):
    """
    Plot mask in domain.
    :param fxlon: longitude coordinates grid
    :param fxlat: latitude coordinates grid
    :param mask: mask to plot
    :param show: optional, show the plot
    """
    plt.pcolormesh(fxlon,fxlat,mask,cmap='Greys',vmin=0,vmax=2)
    bbox = [fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max()]
    set_bbox(bbox)
    if show:
        plt.show()

if __name__ == '__main__':
    from perimeters import Perimeters
    import sys
    if len(sys.argv) < 2:
        print('ERROR: {} perim_path [perim_path ...]'.format(sys.argv[0]))
    perims_path = sys.argv[1:]
    plot_detections('arcgis_hotspots.pkl')
    perims = Perimeters(perims_path)
    perims.plot()
    bbox = read_bbox('wrf_init')
    set_bbox(bbox,show=True)