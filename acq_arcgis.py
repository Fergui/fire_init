try:
    from .tools import read_wrfinfo, save_pkl
    from .geometry import alpha_shape, coords_to_polys, polys_to_coords, featureset_to_coords, simplify_coords, merc_to_lonlat, lonlat_to_merc
except:
    from tools import read_wrfinfo, save_pkl
    from geometry import alpha_shape, coords_to_polys, polys_to_coords, featureset_to_coords, simplify_coords, merc_to_lonlat, lonlat_to_merc
from datetime import datetime, timedelta, timezone
import pandas as pd
import numpy as np
import arcgis
import pickle
import logging
import time

def retry_query(layer, where=None, geometry_filter=None, max_retries=5):
    for retry in range(max_retries):
        try:
            if where is None:
                feature_set = layer.query(geometry_filter=geometry_filter)
            else:
                feature_set = layer.query(where=where,geometry_filter=geometry_filter)
        except Exception as e:
            if retry == max_retries-1:
                print('Max retries to the server with exception {}'.format(e))
                return
            time.sleep(5)
            pass
    return feature_set

def acq_arcgis_past_perims(bbox, now=datetime.utcnow().replace(tzinfo=timezone.utc), max_retries=5):
    # number of years
    #nyears = 3
    # transform bbox to mercator
    xmin,ymin = lonlat_to_merc(bbox[0],bbox[2])
    xmax,ymax = lonlat_to_merc(bbox[1],bbox[3])
    # define AOI
    aoi_latlon = {
            'xmax': float(bbox[1]), 'xmin': float(bbox[0]), 
            'ymax': float(bbox[3]), 'ymin': float(bbox[2])
    }
    aoi_merc = {
            'xmax': xmax, 'xmin': xmin, 
            'ymax': ymax, 'ymin': ymin
    }
    # set parameters to reduce complexity of multipolygons
    params = {
        'min_inner_coords': 50,     # minimum number of coordinates for an inner polygon
        'max_inner_coords': 1000,   # maximum number of coordinates for an inner polygon
        'min_outer_coords': 500,    # minimum number of coordinates for an outer polygon
        'max_outer_coords': 10000   # maximum number of coordinates for an outer polygon 
    }
    # define temporal interval
    #time_from = now-timedelta(days=nyears*365+90)
    #time_to = now-timedelta(days=90)
    # create temporal filter
    #time_fil = '(poly_DateCurrent > DATE \'{:04d}-{:02d}-{:02d}\') AND (poly_DateCurrent < DATE \'{:04d}-{:02d}-{:02d}\')'.format(time_from.year,time_from.month,time_from.day,time_to.year,time_to.month,time_to.day)
    # query IR fire perimeters in bounding box
    # GIS object
    gis = arcgis.gis.GIS()
    # get Historic Geomac Perimeters Combined 2000-2018 dataset
    dataset = gis.content.get('72a928a667634be1b795aa76a61e95f8')
    layer = dataset.layers[19]
    feature_set = retry_query(layer, geometry_filter=arcgis.geometry.filters.intersects(aoi_merc), max_retries=max_retries)
    feature_set_1 = feature_set.sdf.reset_index(drop=True)
    # get Historic Geomac Perimeters 2019 dataset
    dataset = gis.content.get('a829aefbe4e5471490d8f3d47ca5410d')
    layer = dataset.layers[0]
    feature_set = retry_query(layer, geometry_filter=arcgis.geometry.filters.intersects(aoi_merc), max_retries=max_retries)
    feature_set_2 = feature_set.sdf.reset_index(drop=True)
    feature_set = pd.concat((feature_set_1, feature_set_2))
    coords = featureset_to_coords(feature_set)
    coords = simplify_coords(coords, **params)
    polys = coords_to_polys(coords)
    transf_polys = merc_to_lonlat(polys)
    transf_coords = polys_to_coords(transf_polys)
    # get Fire History Perimeters dataset
    dataset = gis.content.get('585b8ff97f5c45fe924d3a1221b446c6')
    layer = dataset.layers[0]
    feature_set = retry_query(layer, geometry_filter=arcgis.geometry.filters.intersects(aoi_latlon), max_retries=max_retries)
    #feature_set = retry_query(layer, where=time_fil, geometry_filter=arcgis.geometry.filters.intersects(aoi_latlon), max_retries=max_retries)
    feature_set = feature_set.sdf.reset_index(drop=True)
    coords = featureset_to_coords(feature_set)
    coords = simplify_coords(coords, **params)
    coords += transf_coords
    if len(coords):
        # save pickle file 
        with open('arcgis_past_perims.pkl','wb') as f:
            pickle.dump((now,coords), f, protocol=3)
    else: 
        logging.warning('any past perimeter found in the search, past perimeters will be skipped')
    with open('arcgis_past_perims_{:02d}{:02d}_{:02d}z.pkl'.format(now.month,now.day,now.hour),'wb') as f:
        pickle.dump((now,coords), f, protocol=3)
    return now, coords

def acq_arcgis_perims(bbox, now=datetime.utcnow().replace(tzinfo=timezone.utc), max_retries=5):
    time = now-timedelta(days=2)
    # set parameters to reduce complexity of multipolygons
    params = {
        'min_inner_coords': 20,     # minimum number of coordinates for an inner polygon
        'max_inner_coords': 1000,   # maximum number of coordinates for an inner polygon
        'min_outer_coords': 100,    # minimum number of coordinates for an outer polygon
        'max_outer_coords': 10000   # maximum number of coordinates for an outer polygon 
    }
    # GIS object
    gis = arcgis.gis.GIS()
    # define AOI
    aoi = {'xmax': float(bbox[1]), 'xmin': float(bbox[0]), 'ymax': float(bbox[3]), 'ymin': float(bbox[2])}
    # get dataset
    dataset = gis.content.get('2191f997056547bd9dc530ab9866ab61')
    # get layer
    layer = dataset.layers[0]  
    # create temporal filter
    time_fil = '(poly_DateCurrent >= DATE \'{0:04d}-{1:02d}-{2:02d}\')'.format(time.year,time.month,time.day)  
    # query IR fire perimeters in bounding box
    feature_set = retry_query(layer, where=time_fil, geometry_filter=arcgis.geometry.filters.intersects(aoi), max_retries=max_retries)
    # reset index
    feature_set = feature_set.sdf.reset_index(drop=True)
    if 'poly_DateCurrent' in feature_set.keys():
        # modified datetime
        perim_time = feature_set['poly_DateCurrent'].sort_values(ascending=False).iloc[0].to_pydatetime(warn=False).replace(tzinfo=timezone.utc)
    else:
        perim_time = None
    # create SHAPE without geometry
    coords = featureset_to_coords(feature_set)
    coords = simplify_coords(coords, **params)
    if len(coords):
        # save pickle file 
        with open('perim1.pkl','wb') as f:
            pickle.dump((perim_time,coords), f, protocol=3)
    else: 
        logging.warning('any current perimeter found in the search, perim1.pkl will be the previous perimeter')
    with open('arcgis_perims_{:02d}{:02d}_{:02d}z.pkl'.format(now.month,now.day,now.hour),'wb') as f:
        pickle.dump((perim_time,coords), f, protocol=3)
    return perim_time, coords

def acq_arcgis_viirs(bbox, now=datetime.utcnow().replace(tzinfo=timezone.utc), max_retries=5):
    time = now-timedelta(days=2)
    # set parameters to reduce complexity of multipolygons
    params = {
        'min_inner_coords': 20,     # minimum number of coordinates for an inner polygon
        'max_inner_coords': 1000,   # maximum number of coordinates for an inner polygon
        'min_outer_coords': 4,      # minimum number of coordinates for an outer polygon
        'max_outer_coords': 10000   # maximum number of coordinates for an outer polygon
    }
    # GIS object
    gis = arcgis.gis.GIS()
    # transform bbox to mercator
    xmin,ymin = lonlat_to_merc(bbox[0],bbox[2])
    xmax,ymax = lonlat_to_merc(bbox[1],bbox[3])
    # define AOI
    aoi = {'xmax': xmax, 'xmin': xmin, 'ymax': ymax, 'ymin': ymin}
    # get dataset
    dataset = gis.content.get('dece90af1a0242dcbf0ca36d30276aa3')
    # get layer
    layer = dataset.layers[0]
    # create temporal filter
    time_fil = '(acq_time >= DATE \'{0:04d}-{1:02d}-{2:02d}\')'.format(time.year,time.month,time.day)
    # query fire detections in bounding box
    feature_set = retry_query(layer, where=time_fil, geometry_filter=arcgis.geometry.filters.intersects(aoi), max_retries=max_retries)
    # reset index
    feature_set = feature_set.sdf.reset_index(drop=True)
    # extract coordinates
    coords = [[shape['x'],shape['y']] for shape in feature_set.get('SHAPE',feature_set.get('shape', []))]
    # create shape of type object
    shape = [{'x': s[0], 'y': s[1], 'spatialReference': {'wkid': 4326}} for s in coords]
    # new column SHAPE with the final transformed coordinates
    feature_set['SHAPE'] = shape
    feature_set['SHAPE'] = feature_set['SHAPE'].astype(object)
    logging.info('     found {} hotspots'.format(len(feature_set)))
    # save pickle file 
    out_file = 'arcgis_hotspots_{:02d}{:02d}_{:02d}z.pkl'.format(now.month,now.day,now.hour)
    feature_set.to_pickle(out_file, protocol=4)
    feature_set.to_pickle('arcgis_hotspots.pkl', protocol=4)
    alpha = 375
    x = feature_set['SHAPE'].map(lambda row: row['x'])
    y = feature_set['SHAPE'].map(lambda row: row['y'])
    points = np.vstack([x,y]).T
    if len(points) < 4:
        logging.warning('{} hotspots found, less than 4 so skipping alphashape'.format(len(points)))
        return out_file
    # last detection time
    t = feature_set['esritimeutc'].sort_values(ascending=False).iloc[0].to_pydatetime(warn=False).replace(tzinfo=timezone.utc)
    # alpha shape
    perims = alpha_shape(points, alpha=alpha, only_outer=True)
    # clean the perimeters
    if len(perims):
        perims = simplify_coords(perims,**params)
        polys = coords_to_polys(perims)
        transf_polys = merc_to_lonlat(polys)
        transf_perims = polys_to_coords(transf_polys)
        # save perimeters
        if len(transf_perims):
            save_pkl((t,transf_perims),'perim2.pkl')
    else:
        transf_perims = []
    save_pkl((t,transf_perims),'hotspots_perims_{:02d}{:02d}_{:02d}z.pkl'.format(t.month,t.day,t.hour))
    return out_file
 
if __name__ == '__main__':
    import os.path as osp
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    wrf_path = 'wrf_init'
    ifm,ifn,ffm,ffn,fxlon,fxlat = read_wrfinfo(wrf_path) 
    bbox = (fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max())
    logging.info('bbox {}'.format(bbox))
    now = datetime.utcnow().replace(tzinfo=timezone.utc)
    if not osp.exists('arcgis_past_perims.pkl'):
        # acquisiton of last fire perimeter
        logging.info('acquisition of past fire perimeters')
        acq_arcgis_past_perims(bbox, now)
    # acquisiton of last fire perimeter
    logging.info('acquisition of current fire perimeters')
    acq_arcgis_perims(bbox, now)
    # acquisition of VIIRS fire detections
    logging.info('acquisition of fire detections')
    acq_arcgis_viirs(bbox, now)
