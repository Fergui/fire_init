import logging
import numpy as np
import os.path as osp
import pickle
from datetime import datetime, timedelta, timezone
from utils import read_wrfinfo
from geometry import mask_perim
from utils import integrate_init, add_smoke

def perims_interp(perim1, perim2, fxlon, fxlat, **params):
    # resolving parameters
    time_step = params.get('time_step', 60.)
    prev_PERIMS = params.get('prev_perims',{})
    perim1_time = params.get('perim1_time', -79200.)
    perim2_time = params.get('perim2_time', 7200.)
    outside_time = params.get('outside_time', 360000.)
    boundary_margin = params.get('boundary_margin', 1)
    past_perims = params.get('past_perims', None)
    scars_mask_path = params.get('scars_mask_path', 'scars_mask.pkl')
    # defining points in the grid
    points = np.c_[fxlon.ravel(),fxlat.ravel()]
    # process past perimeters
    processed_past_perims = False
    if osp.exists(scars_mask_path): 
        logging.info(f'found scars mask in {scars_mask_path}')
        with open(scars_mask_path, 'rb') as f:
            insidepastperims = pickle.load(f)
        processed_past_perims = True
    elif past_perims is not None:
        # find points inside and outside past perimeters (defining masks)
        logging.info('finding mask for past perimeters')
        insidepastperims = mask_perim(past_perims,points)
        with open(scars_mask_path, 'wb') as f:
            pickle.dump(insidepastperims,f)
        processed_past_perims = True
    # find points inside and outside perimeters (defining masks)
    logging.info('finding mask for perimeter 1')
    insideperim1 = mask_perim(perim1,points)
    insideperim1 = np.reshape(insideperim1,fxlon.shape)
    if len(prev_PERIMS):
        logging.info('finding mask for previous perimeter 1')
        insideprevperim1 = mask_perim(prev_PERIMS['perim1'],points)
        insideprevperim1 = np.reshape(insideprevperim1,fxlon.shape)
        insideperim1 = np.logical_or(insideperim1,insideprevperim1)
    logging.info('finding mask for perimeter 2')
    insideperim2 = mask_perim(perim2,points)
    insideperim2 = np.reshape(insideperim2,fxlon.shape)
    insideperim2 = np.logical_or(insideperim1,insideperim2)
    if len(prev_PERIMS):
        logging.info('finding mask for previous perimeter 2')
        insideprevperim2 = mask_perim(prev_PERIMS['perim2'],points)
        insideprevperim2 = np.reshape(insideprevperim2,fxlon.shape)
        insideperim2 = np.logical_and(insideperim2,~insideprevperim2)
    outsideperim1 = ~insideperim1
    outsideperim2 = ~insideperim2
    # define active and inactive regions and coordinates
    logging.info('define active and inactive regions')
    active = np.logical_and(outsideperim1,insideperim2)
    num_active = (active==1).sum()
    fxlon_active = fxlon[active]
    fxlat_active = fxlat[active]
    # initialize fire arrival time
    TIGN_G = insideperim1*perim1_time+outsideperim2*perim2_time
    # concatenate all coordinates from all outer boundary polygons
    perim1_lons = np.array([c[0] for p in perim1 for c in p[0]])
    perim1_lats = np.array([c[1] for p in perim1 for c in p[0]])
    perim2_lons = np.array([c[0] for p in perim2 for c in p[0]])
    perim2_lats = np.array([c[1] for p in perim2 for c in p[0]])
    # compute distance between each active point compared to the closest first perimeter point
    pp = np.tile(perim1_lons,(num_active,1))
    d1 = (pp-fxlon_active[:,np.newaxis])**2
    pp = np.tile(perim1_lats,(num_active,1))
    d1 = (d1 + (pp-fxlat_active[:,np.newaxis])**2).min(axis=1)
    # compute distance between each active point compared to the closest second perimeter point
    pp = np.tile(perim2_lons,(num_active,1))
    d2 = (pp-fxlon_active[:,np.newaxis])**2
    pp = np.tile(perim2_lats,(num_active,1))
    d2 = (d2 + (pp-fxlat_active[:,np.newaxis])**2).min(axis=1)
    # compute fire arrival time at the active points as a linear interpolation between perimeter times weighted by distances
    eps = 1e-8
    TIGN_G_active = d2/(d1+d2+eps)*perim1_time+d1/(d1+d2+eps)*perim2_time
    # substitute the points in the active region, using the interpolated data
    TIGN_G_final = TIGN_G.copy()
    TIGN_G_final[active] = TIGN_G_active                    # use interpolated data in the active region
    TIGN_G_final[TIGN_G_final<perim1_time] = perim1_time    # set points less then the first perimeter time to the first perimeter time
    TIGN_G_final[outsideperim2] = outside_time              # set large value outside of the second perimeter so this region won't get ignited
    TIGN_G_final[TIGN_G_final>outside_time] = outside_time  # if interpolation gave crazy high numbers set them to outside time.
    TIGN_G_final[TIGN_G_final<time_step] = time_step        # if interpolation gave crazy small numbers set them to time step.
    # calculate fuel mask
    logging.info('calculate fuel mask')
    FUEL_MASK = insideperim1.copy()
    FUEL_MASK = FUEL_MASK.astype(bool)
    inds = np.argwhere(FUEL_MASK)
    for i in range(len(inds)):
        FUEL_MASK[inds[i,0]-boundary_margin:inds[i,0]+boundary_margin,inds[i,1]-boundary_margin:inds[i,1]+boundary_margin] = True
    if processed_past_perims:
        FUEL_MASK[insidepastperims.reshape(FUEL_MASK.shape)] = True
    return TIGN_G_final,FUEL_MASK

if __name__ == '__main__':
    from perimeters import Perimeter
    import sys
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    perim1_path = sys.argv[1]
    perim2_path = sys.argv[2]
    start_time = sys.argv[3]
    wrf_path = 'wrf_init'
    wrfinput_path = 'wrfinput_d03'
    prev_PERIMS = {}
    if len(sys.argv) > 4:
        prev_results = sys.argv[4]
        try:
            with open(prev_results,'rb') as f:
                r = pickle.load(f)
            prev_PERIMS = r.get('PERIMS',{})
        except Exception as e:
            logging.warning('something wrong when getting previous perimeters with exception {}'.format(e))
    stime = datetime.strptime(start_time,'%Y%m%d%H').replace(tzinfo=timezone.utc)
    dtime = stime+timedelta(hours=3)
    result_path = 'wszystko_{:02d}{:02d}_{:02d}z'.format(dtime.month,dtime.day,dtime.hour)
    logging.info(f'interpolating between {perim1_path} and {perim2_path} into a result file {result_path}.pkl')
    # parameters (normally static)
    params = {
        'wrf_path': wrf_path,  # wrf path for getting grid of points
        'prev_perims': prev_PERIMS,  # previous perimeters from previous estimation
        'time_step': 60.,  # minimal time of the fire arrival time for WRF-SFIRE
        'perim1_time': -36000.,  # time of the first perimeter 10h before the sim start (07.22.2021 00UTC); 10h = -36000s
        'perim2_time': 7200.,  # time of the second perimeter from the start of the run (07.22.2021 12UTC); 2h = 7200s
        'outside_time': 360000.,  # a large number outside of the last perimeter should be greater than fire_perimeter_time in namelist.input
        'boundary_margin': 5,  # margin in grid points defining the region outside of the first perimeter where the fuel will be removed
        'include_boundary': False,  # include boundary in processing fire grid
        'past_perims_path': 'arcgis_past_perims.pkl',  # path to past perims to create scars mask
        'scars_mask_path': 'scars_mask.pkl',  # path to past scars mask
        'integrate_now': False  # wrfinput already provided and ready to integrate 
    }
    # read fire grid and bounding box from it
    ifm,ifn,ffm,ffn,fxlon,fxlat = read_wrfinfo(wrf_path,boundary=params.get('include_boundary', False))
    params.update({'ifm': ifm, 'ifn': ifn, 'ffm': ffm, 'ffn': ffn})
    bbox = [fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max()]
    # getting IR perimeters
    logging.info('getting IR fire perimeters')
    perim1 = Perimeter(perim1_path)
    perim2 = Perimeter(perim2_path)
    if 'scars_mask_path' in params and not osp.exists(params['scars_mask_path']) and 'past_perims_path' in params:
        past_perim = Perimeter(params['past_perims_path'])
        params.update({'past_perim': past_perim.coords}) 
    # if the time is included in the perimeters
    c1 = perim1.time is not None
    c2 = perim2.time is not None
    if c1 and c2:
        # difference in time from both perimeters, always larger or equal than 0 and set a maximum of -10h=-36000s
        params['perim1_time'] = -min(max(0,(perim2[0]-perim1[0]).total_seconds()),36000.)
    PERIMS = {'perim1': perim1, 'perim2': perim2}
    # calculate tign_g and fuel_mask
    TIGN_G,FUEL_MASK = perims_interp(perim1.coords, perim2.coords, fxlon, fxlat, **params)
    with open(result_path+'.pkl','wb') as f: 
        pickle.dump({'PERIMS': PERIMS, 'TIGN_G': TIGN_G, 'FUEL_MASK': FUEL_MASK, 'params': params, 'fxlon': fxlon, 'fxlat': fxlat},f)
    # integrate perimeters and add smoke
    if params['integrate_now']:
        integrate_init(wrfinput_path, TIGN_G, FUEL_MASK, **params)
        add_smoke()
    logging.info('DONE')
