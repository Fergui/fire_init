try:
    from .geometry import mask_perim, fire_interp, lonlat_to_merc, merc_to_lonlat, polys_to_coords, coords_to_polys
    from .tools import read_wrfinfo, load_pkl, save_pkl, integrate_init, add_smoke
except:
    from geometry import mask_perim, fire_interp, lonlat_to_merc, merc_to_lonlat, polys_to_coords, coords_to_polys
    from tools import read_wrfinfo, load_pkl, save_pkl, integrate_init, add_smoke
from datetime import datetime, timedelta, timezone
import os.path as osp
import numpy as np
import logging

def perims_interp(perim1, perim2, fxlon, fxlat, **params):
    """
    Perimeter interpolation and fuel masking.
    :param perim1: first perimeter coordinates (list of lists with outer and inner coordinates)
    :param perim2: second perimeter coordinates (list of lists with outer and inner coordinates)
    :param fxlon: longitude coordinates grid
    :param fxlat: latitude coordinates grid
    :param params: dictonary of parameters to interpolate fire arrival time and fuel mask
    :return: fire arrival time, fuel mask, and perimeter masks
    """
    # resolving parameters
    simplify_tol = params.get('simplify_tol', 1e-4)
    past_perims = params.get('past_perims', None)
    prev_perims = params.get('prev_perims', None)
    fuel_method = params.get('fuel_method', 1)
    boundary_margin = params.get('boundary_margin', 1)
    ir_buffer_dist = params.get('ir_buffer_dist', 100.)
    sat_buffer_dist = params.get('sat_buffer_dist', 400.)
    scars_mask_path = params.get('scars_mask_path', 'scars_mask.pkl')
    # defining points in the grid
    points = np.c_[fxlon.ravel(),fxlat.ravel()]
    # process past perimeters
    processed_past_perims = False
    if osp.exists(scars_mask_path): 
        logging.info(f'found scars mask in {scars_mask_path}')
        insidepastperims = load_pkl(scars_mask_path)
        processed_past_perims = True
    elif past_perims is not None:
        # find points inside and outside past perimeters (defining masks)
        logging.info('finding mask for past perimeters')
        past_perims = polys_to_coords(coords_to_polys(past_perims.coords).simplify(simplify_tol))
        insidepastperims = mask_perim(past_perims,points)
        save_pkl(insidepastperims,scars_mask_path)
        processed_past_perims = True
    # process previous perimeters and masks
    processed_prev_perims = False
    if prev_perims is not None:
        prev_PERIMS = {}
        prev_MASKS = {}
        prev_BUFF_MASKS = {}
        try:
            data = load_pkl(prev_perims)
            prev_PERIMS = data.get('PERIMS',{})
            prev_MASKS = data.get('PERIM_MASKS',{})
            prev_BUFF_MASKS = data.get('BUFF_MASKS',{})
            processed_prev_perims = True
        except Exception as e:
            logging.warning('something wrong when getting previous perimeters with exception {}'.format(e))
    # find points inside and outside perimeters (defining masks)
    logging.info('finding mask for perimeter 1')
    perim1 = polys_to_coords(coords_to_polys(perim1).simplify(simplify_tol))
    insideperim1 = mask_perim(perim1,points)
    insideperim1 = np.reshape(insideperim1,fxlon.shape)
    if prev_perims is not None:
        logging.info('found mask for previous perimeter 1')
        insideprevperim1 = prev_MASKS['mask1']
        insideperim1 = np.logical_or(insideperim1,insideprevperim1)
    logging.info('finding mask for perimeter 2')
    perim2 = polys_to_coords(coords_to_polys(perim2).simplify(simplify_tol))
    insideperim2 = mask_perim(perim2,points)
    insideperim2 = np.reshape(insideperim2,fxlon.shape)
    insideperim2 = np.logical_or(insideperim1,insideperim2)
    PERIM_MASKS = {'mask1': insideperim1, 'mask2': insideperim2}
    # calculate tign_g
    logging.info('calculate tign_g')
    TIGN_G = fire_interp(insideperim1, insideperim2, perim1, perim2, fxlon, fxlat, **params)
    # calculate fuel mask
    BUFF_MASKS = {}
    logging.info('calculate fuel mask')
    if fuel_method == 0:
        FUEL_MASK = insideperim1.copy()
        FUEL_MASK = FUEL_MASK.astype(bool)
        inds = np.argwhere(FUEL_MASK)
        for i in range(len(inds)):
            FUEL_MASK[inds[i,0]-boundary_margin:inds[i,0]+boundary_margin,inds[i,1]-boundary_margin:inds[i,1]+boundary_margin] = True
        if processed_past_perims:
            FUEL_MASK[insidepastperims.reshape(FUEL_MASK.shape)] = True
    elif fuel_method == 1:
        logging.info('finding mask for perimeter 1 buffer')
        perim1buff = polys_to_coords(merc_to_lonlat(lonlat_to_merc(coords_to_polys(perim1)).buffer(ir_buffer_dist)).simplify(simplify_tol))
        insideperim1buff = mask_perim(perim1buff,points)
        insideperim1buff = np.reshape(insideperim1buff,fxlon.shape)
        if processed_prev_perims:
            insideprevperim1buff = prev_BUFF_MASKS['mask1']
            insideperim1buff = np.logical_or(insideperim1buff,insideprevperim1buff)
            logging.info('finding mask for previous perimeter 2 buffer')
            prev_perim2 = prev_PERIMS['perim2'].coords
            prev_perim2_buff = polys_to_coords(merc_to_lonlat(lonlat_to_merc(coords_to_polys(prev_perim2)).buffer(sat_buffer_dist)).simplify(simplify_tol))
            insideprevperim2buff = mask_perim(prev_perim2_buff,points)
            insideprevperim2buff = np.reshape(insideprevperim2buff,fxlon.shape)
        else:
            insideprevperim2buff = np.zeros(fxlon.shape)
        BUFF_MASKS.update({'mask1': insideperim1buff, 'mask2': insideprevperim2buff})
        FUEL_MASK = np.logical_or(insideperim1buff.astype(bool),insideprevperim2buff.astype(bool))
        FUEL_MASK = FUEL_MASK.astype(bool)
        if processed_past_perims:
            FUEL_MASK[insidepastperims.reshape(FUEL_MASK.shape)] = True
    else:
        logging.error()
    return TIGN_G,FUEL_MASK,PERIM_MASKS,BUFF_MASKS

if __name__ == '__main__':
    from perimeters import Perimeter
    import sys
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    if len(sys.argv) < 4:
        logging.error('{} perim1_path perim2_path YYYYMMDDHH [prev_results_path]'.format(sys.argv[0]))
        sys.exit()
    perim1_path = sys.argv[1]
    perim2_path = sys.argv[2]
    start_time = sys.argv[3]
    wrf_path = 'wrf_init'
    wrfinput_path = 'wrfinput_d03'
    
    stime = datetime.strptime(start_time,'%Y%m%d%H').replace(tzinfo=timezone.utc)
    dtime = stime+timedelta(hours=3)
    result_path = 'wszystko_{:02d}{:02d}_{:02d}z'.format(dtime.month,dtime.day,dtime.hour)
    logging.info(f'interpolating between {perim1_path} and {perim2_path} into a result file {result_path}.pkl')
    # parameters (normally static)
    params = {
        'wrf_path': wrf_path,  # wrf path for getting grid of points
        'simplify_tol': 1e-4,  # simplify tolerance for optimization of masking 
        'time_step': 60.,  # minimal time of the fire arrival time for WRF-SFIRE
        'perim1_time': -36000.,  # time of the first perimeter 10h before the sim start (07.22.2021 00UTC); 10h = -36000s
        'perim2_time': 7200.,  # time of the second perimeter from the start of the run (07.22.2021 12UTC); 2h = 7200s
        'outside_time': 360000.,  # a large number outside of the last perimeter should be greater than fire_perimeter_time in namelist.input
        'fuel_method': 1, # fuel removal method: 0 - grid boundary margin, 1 - buffer perimeter margin
        'boundary_margin': 5,  # margin in grid points defining the region outside of the first perimeter where the fuel will be removed (only used for fuel_method 0)
        'ir_buffer_dist': 200., # distance in meters to create the buffer for IR perimeters to define masking of fuel (only used for fuel_method 1)
        # maybe change it to min_ros*(perim2.time-perim1.time).total_seconds()? With 2h=7200s and ros_min=0.01m/s => ir_buffer_dist=72m
        'sat_buffer_dist': 500., # distance in meters to create the buffer for IR perimeters to define masking of fuel (only used for fuel_method 1)
        # maybe change it to maximal fire detection resolution? So, VIIRS are 350m but adding angle can get to 500m easily
        'past_perims_path': 'arcgis_past_perims.pkl',  # path to past perims to create scars mask
        'scars_mask_path': 'scars_mask.pkl',  # path to past scars mask
        'integrate_now': False  # wrfinput already provided and ready to integrate 
    }
    # adding previous results if given
    if len(sys.argv) > 4:
        params.update({'prev_perims': sys.argv[4]}) # previous results to use perimeters and masks
    # read fire grid and bounding box from it
    ifm,ifn,ffm,ffn,fxlon,fxlat = read_wrfinfo(wrf_path)
    params.update({'ifm': ifm, 'ifn': ifn, 'ffm': ffm, 'ffn': ffn})
    bbox = [fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max()]
    # getting IR perimeters
    logging.info('getting IR fire perimeters')
    perim1 = Perimeter(perim1_path)
    perim2 = Perimeter(perim2_path)
    if 'scars_mask_path' in params and not osp.exists(params['scars_mask_path']) and 'past_perims_path' in params:
        past_perim = Perimeter(params['past_perims_path'])
        params.update({'past_perims': past_perim}) 
    # if the time is included in the perimeters
    c1 = perim1.time is not None
    c2 = perim2.time is not None
    if c1 and c2:
        # difference in time from both perimeters, always larger or equal than 0 and set a maximum of -10h=-36000s
        params['perim1_time'] = -min(max(0,(perim2.time-perim1.time).total_seconds()),36000.)
    PERIMS = {'perim1': perim1, 'perim2': perim2}
    # calculate tign_g and fuel_mask
    TIGN_G,FUEL_MASK,PERIM_MASKS,BUFF_MASKS = perims_interp(perim1.coords, perim2.coords, fxlon, fxlat, **params)
    save_pkl({'PERIMS': PERIMS, 'PERIM_MASKS': PERIM_MASKS, 
              'TIGN_G': TIGN_G, 'FUEL_MASK': FUEL_MASK, 
              'BUFF_MASKS': BUFF_MASKS, 'params': params, 
              'fxlon': fxlon, 'fxlat': fxlat}, result_path+'.pkl')
    # integrate perimeters and add smoke
    if params['integrate_now']:
        integrate_init(wrfinput_path, TIGN_G, FUEL_MASK, **params)
        add_smoke(
            ['wrfout_d01','wrfout_d02','wrfout_d03'],
            ['wrfinput_d01','wrfinput_d02','wrfinput_d03']
        )
    logging.info('DONE')
