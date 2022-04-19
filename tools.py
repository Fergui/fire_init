import logging
import pickle
import netCDF4 as nc
import numpy as np
from lxml import etree
import os.path as osp

def save_pkl(data, path):
    """
    Save data to pickle file.
    :param data: python object to save
    :param path: pickle file path to save the data on
    """
    with open(path,'wb') as f:
        pickle.dump(data,f,protocol=3)

def load_pkl(path):
    """
    Load data from pickle file.
    :param path: pickle file path to load the data from
    :return: python object in the pickle file
    """
    with open(path,'rb') as f:
        r = pickle.load(f)
    return r

def read_wrfinfo(filePath):
    """
    Read WRF netCDF file information.
    :param filePath: path to WRF netCDF file
    :return: indices without strip, and longitude and latitude coordinates grids 
    """
    logging.info('reading wrfout fire grid')
    # open netCDF file
    with nc.Dataset(filePath) as d:
        # atmospheric dimensions
        m,n = d.variables['XLONG'].shape[1:]
        # fire dimensions
        fm,fn = d.variables['FXLONG'].shape[1:]
        # dimensions corrected for extra strips
        ifm=int(fm/(m+1)/2)
        ifn=int(fn/(n+1)/2)
        ffm=int(fm-1.5*fm/(m+1))
        ffn=int(fn-1.5*fm/(m+1))
        # read coordinates masking extra strip
        fxlon = d.variables['FXLONG'][0,ifm:ffm,ifn:ffn].data
        fxlat = d.variables['FXLAT'][0,ifm:ffm,ifn:ffn].data
    return ifm,ifn,ffm,ffn,fxlon,fxlat

def read_bbox(path):
    """
    Read bounding box from WRF netCDF file.
    :param path: path to WRF netCDF file
    :return: bounding box coordinates [min_lon,max_lon,min_lat,max_lat]
    """
    logging.info('reading wrfout bbox')
    _,_,_,_,fxlon,fxlat = read_wrfinfo(path)
    bbox = (fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max())
    return bbox

def read_kml(path):
    """
    Read KML file with polygon coordinates.
    :param path: path to KML file with polygon coordinates
    :return: polygons in format list of lists with outer and inner coordinates
    """
    # parse XLM file
    tree = etree.parse(path)
    # get its root to start finding elements
    root = tree.getroot()
    # get namespace map that each tag is going to contain
    nsmap = root.nsmap.get('kml',root.nsmap.get(None))
    # create xpath lambda function to generate paths to elements with the namespace map
    xpath = lambda tag: './/{{{}}}{}'.format(nsmap,tag) if nsmap else './/{}'.format(tag)
    polys = []
    for pm in root.iterfind(xpath('Placemark')): 
        for pp in pm.iterfind(xpath('Polygon')):
            # find outerBoundaryIs coordinates element
            out_elem = pp.find(xpath('outerBoundaryIs')).find(xpath('coordinates'))
            # parse outer coordinates from the text of the outer coordinates element
            out_coords = np.array([np.array(coord.split(',')[:2]).astype(float) for coord in out_elem.text.split()])
            poly = [out_coords]
            # for each innerBoundaryIs
            for inn_elems in pp.iterfind(xpath('innerBoundaryIs')):
                # for each coordinates element
                for inn_elem in inn_elems.iterfind(xpath('coordinates')):
                    # parse inner coordinates from the text of the inner coordinates element
                    inn_coords = np.array([np.array(coord.split(',')[:2]).astype(float) for coord in inn_elem.text.split()])
                    poly.append(inn_coords)
            polys.append(poly)
    return polys

def integrate_init(wrfinput_path, TIGN_G, FUEL_MASK, outside_time=360000.):
    """
    Integrate fire arrival time and fuel mask in WRF netCDF file.
    :param wrfinput_path: path to WRF netCDF file
    :param TIGN_G: fire arrival time matrix
    :param FUEL_MASK: fuel mask matrix
    :param outside_time: optional, time outside of fire arrival time
    """
    logging.info('integrate initialization')
    ifm,ifn,ffm,ffn,_,_ = read_wrfinfo(wrfinput_path)
    # open wrfinput file in append mode
    with nc.Dataset(wrfinput_path,'a') as d:
        # read existing fire arrival time and fuel
        TIGN_G_ = d.variables['TIGN_G'][0].data
        NFUEL_CAT_ = d.variables['NFUEL_CAT'][0].data
        # masking unburnable places in constructed fire arrival time
        TIGN_G[NFUEL_CAT_[ifm:ffm,ifn:ffn] == 14.] = outside_time 
        # fill fire arrival time in wrfinput file
        TIGN_G_ = np.ones(TIGN_G_.shape)*outside_time
        TIGN_G_[ifm:ffm,ifn:ffn] = TIGN_G
        d.variables['TIGN_G'][0] = TIGN_G_
        # fill fuel mask to category 14 in wrfinput file
        NFUEL_CAT = NFUEL_CAT_[ifm:ffm,ifn:ffn]
        NFUEL_CAT[FUEL_MASK] = 14.
        NFUEL_CAT_[ifm:ffm,ifn:ffn] = NFUEL_CAT
        d.variables['NFUEL_CAT'][0] = NFUEL_CAT_
    
def add_smoke(wrfout_paths, wrfinput_paths):
    """
    Add smoke from wrfout paths to wrfinput paths.
    :param wrfout_paths: list of paths to WRF netCDF output files
    :param wrfinput_paths: list of paths to WRF netCDF input files
    """
    logging.info('add smoke')
    if any([not osp.exists(wp) for wp in wrfout_paths]):
        logging.warning('missing some wrfout file, so skipping')
        return
    assert len(wrfout_paths) == len(wrfinput_paths), 'missing some wrfout file, so skipping'
    # smoke from previous forecast
    for i in range(len(wrfout_paths)):
        logging.info('smoke in domain {}'.format(i+1))
        with nc.Dataset(wrfout_paths[i]) as d:
            tr17_1 = d['tr17_1'][...]
        if d['tr17_1'][:].shape != tr17_1.shape:
            logging.warning('size of domains do not correspond, add_smoke skipped')
            return
        with nc.Dataset(wrfinput_paths[i],'a') as d:
            d['tr17_1'][:] = tr17_1
