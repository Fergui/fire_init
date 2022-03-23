import logging
import pickle
import os.path as osp
import netCDF4 as nc
import numpy as np
from lxml import etree

def save_pkl(data, path):
    with open(path,'wb') as f:
        pickle.dump(data,f,protocol=3)

def load_pkl(path):
    with open(path,'rb') as f:
        r = pickle.load(f)
    return r

def read_wrfinfo(filePath):
    logging.info('reading wrfout fire grid')
    # open netCDF file
    d = nc.Dataset(filePath)
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
    # close netCDF file
    d.close()
    return ifm,ifn,ffm,ffn,fxlon,fxlat

def read_bbox(path):
    logging.info('reading wrfout bbox')
    _,_,_,_,fxlon,fxlat = read_wrfinfo(path)
    bbox = (fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max())
    return bbox

def read_kml(path):
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
    logging.info('integrate initialization')
    ifm,ifn,ffm,ffn,_,_ = read_wrfinfo(wrfinput_path)
    # open wrfinput file in append mode
    d = nc.Dataset(wrfinput_path,'a')
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
    # close netCDF file
    d.close()
    
def add_smoke():
    logging.info('add smoke')
    if not (osp.exists('wrfout_d01') and osp.exists('wrfout_d02') and osp.exists('wrfout_d03')):
        logging.warning('missing some wrfout file, so skipping')
        return
    # smoke from previous forecast
    logging.info('smoke in domain 1')
    d = nc.Dataset('wrfout_d01')
    tr17_1 = d['tr17_1'][...]
    d.close()
    d = nc.Dataset('wrfinput_d01','a')
    d['tr17_1'][:] = tr17_1
    d.close()
    logging.info('smoke in domain 2')
    d = nc.Dataset('wrfout_d02')
    tr17_1 = d['tr17_1'][...]
    d.close()
    d = nc.Dataset('wrfinput_d02','a')
    d['tr17_1'][:] = tr17_1
    d.close()
    logging.info('smoke in domain 3')
    d = nc.Dataset('wrfout_d03')
    tr17_1 = d['tr17_1'][...]
    d.close()
    d = nc.Dataset('wrfinput_d03','a')
    d['tr17_1'][:] = tr17_1
    d.close()