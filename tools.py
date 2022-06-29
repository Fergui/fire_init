import logging
import pickle
import netCDF4 as nc
import numpy as np
from lxml import etree
import os.path as osp
from datetime import datetime, timezone

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

def parse_placemark(root, xpath):
    """
    Parse polygon coordinates from a Placemark KML object.
    :param root: Placemark KML object parsed with lxml
    :param xpath: Lambda with namespace to add to object name
    :return: list of times and polys to concatenate
    """
    times = []
    polys = []
    # look for time in file
    tspan = root.find(xpath('TimeSpan')) 
    tstamp = root.find(xpath('TimeStamp')) 
    if tspan is not None:
        begin = tspan.find(xpath('begin'))
        end = tspan.find(xpath('end'))
        if begin is not None:
            try:
                time = datetime.strptime(begin.text,'%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=timezone.utc)
            except:
                time = begin.text
        elif end is not None:
            try:
                time = datetime.strptime(end.text,'%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=timezone.utc)
            except:
                time = end.text
        else:
            time = None
    elif tstamp is not None:
        when = tstamp.find(xpath('when'))
        if when is not None:
            try:
                time = datetime.strptime(when.text,'%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=timezone.utc)
            except:
                time = when.text
    else:
        time = None
    for pp in root.iterfind(xpath('Polygon')):
        # find outerBoundaryIs coordinates element
        out_elem = pp.find(xpath('outerBoundaryIs')).find(xpath('coordinates'))
        # parse outer coordinates from the text of the outer coordinates element
        out_coords = np.array([np.char.strip(coord.split(',')[:2]).astype(float) for coord in out_elem.text.strip().replace(' ',',').split('0,') if coord.strip() != ''])
        if len(out_coords) >= 50:
            poly = [out_coords]
            # for each innerBoundaryIs
            for inn_elems in pp.iterfind(xpath('innerBoundaryIs')):
                # for each coordinates element
                for inn_elem in inn_elems.iterfind(xpath('coordinates')):
                    # parse inner coordinates from the text of the inner coordinates element
                    inn_coords = np.array([np.array(coord.split(',')[:2]).astype(float) for coord in inn_elem.text.split()])
                    if len(inn_coords) >= 10:
                        poly.append(inn_coords)
            times.append(time)
            polys.append(poly)
    return times,polys

def read_kml(path):
    """
    Read KML file with polygon coordinates.
    :param path: path to KML file with polygon coordinates
    :return: polygons in format list of lists with outer and inner coordinates
    """
    # parse XML file
    parser = etree.XMLParser(recover=True, remove_blank_text=True)
    tree = etree.parse(path, parser)
    # get its root to start finding elements
    root = tree.getroot()
    # get namespace map that each tag is going to contain
    nsmap = root.nsmap.get('kml',root.nsmap.get(None))
    # create xpath lambda function to generate paths to elements with the namespace map
    xpath = lambda tag: './/{{{}}}{}'.format(nsmap,tag) if nsmap else './/{}'.format(tag)
    polys = []
    times = []
    proc_perims = False
    # start by trying to process placemarks that has some key decription parameters
    for pm in root.iterfind(xpath('Placemark')):
        desc = pm.find(xpath('description'))
        desc = '' if desc is None else desc.text.lower()
        if 'gis_acres' in desc or 'polygon' in desc and 'fc ir polygon type' not in desc or 'ir heat perimeter' in desc and 'area covered' not in desc:
            proc_perims = True
            ptimes,ppolys = parse_placemark(pm, xpath)
            times += ptimes
            polys += ppolys
    # if any perimeter was processed in the previous stage, try to get all placemarks with document name having the word "perimeter"
    if not proc_perims:
        for doc in root.iterfind(xpath('Document')):
            name = doc.find(xpath('name'))
            name = '' if name is None else name.text.lower()
            if 'perimeter' in name:
                proc_perim = True
                for pm in doc.iterfind(xpath('Placemark')):
                    ptimes,ppolys = parse_placemark(pm, xpath)
                    times += ptimes
                    polys += ppolys
    # if any perimeter was processed in the previous two stages, try to get all placemarks having name the word "perimeter"
    if not proc_perims:
        for pm in root.iterfind(xpath('Placemark')):
            name = pm.find(xpath('name'))
            name = '' if name is None else name.text.lower()
            if 'perimeter' in name:
                proc_perim = True
                ptimes,ppolys = parse_placemark(pm, xpath)
                times += ptimes
                polys += ppolys
    # clean unique times and re-organize coordinates accordingly
    times = np.array(times)
    valid_times = times != None
    invalid_polys = []
    for t in np.where(~valid_times)[0]:
        invalid_polys.append(polys[t])
    result_times, indxs = np.unique(times[np.where(valid_times)[0]], return_inverse=True)
    result_times = list(result_times)
    result_polys = []
    for t in np.unique(indxs):
        c = []
        for i in np.where(indxs == t)[0]:
            c.append(polys[i])
        result_polys.append(c)
    result_ids = [id(p) for p in result_polys]
    result_polys = sorted(result_polys, key=lambda x: result_times[result_ids.index(id(x))])
    result_times = sorted(result_times)
    if len(invalid_polys):
        result_times.append(None)
        result_polys.append(invalid_polys)
    return result_times, result_polys

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
    assert len(wrfout_paths) == len(wrfinput_paths), 'missing some wrfout file, so skipping'
    if any([not osp.exists(wp) for wp in wrfout_paths]):
        logging.warning('missing some wrfout file, so skipping')
        return
    with nc.Dataset(wrfout_paths[-1]) as d:
        xlon = d.variables['XLONG'][...]
        xlat = d.variables['XLAT'][...]
        bbox_out = (xlon.min(), xlon.max(), xlat.min(), xlat.max())
    with nc.Dataset(wrfinput_paths[-1]) as d:
        xlon = d.variables['XLONG'][...]
        xlat = d.variables['XLAT'][...]
        bbox_in = (xlon.min(), xlon.max(), xlat.min(), xlat.max())
    if any(abs(np.array(bbox_out)-np.array(bbox_in)) > 1e-5):
        logging.warning('bounding box for wrfout_d03 \n {} \n which is different than wrfinput_d03 \n {} \n add_smoke skipped'.format(bbox_out,bbox_in))
        return
    # smoke from previous forecast
    for i in range(len(wrfout_paths)):
        logging.info('smoke in domain {}'.format(i+1))
        with nc.Dataset(wrfout_paths[i]) as d:
            tr17_1 = d['tr17_1'][...]
        with nc.Dataset(wrfinput_paths[i],'a') as d:
            if d['tr17_1'][...].data.shape != tr17_1.data.shape:
                logging.warning('size of domains do not correspond, add_smoke skipped')
                return
            d['tr17_1'][:] = tr17_1
