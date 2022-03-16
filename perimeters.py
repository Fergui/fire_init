import logging
import pickle
import simplekml
import numpy as np
import pandas as pd
from collections.abc import Iterable
from geometry import coords_to_polys, poly_area
from plot_tools import plot_perim
from utils import read_kml

class PerimeterError(Exception):
    pass

class Perimeter(object):
    """
    Define and process a perimeter
    """
    def __init__(self, info):
        self.path = None
        self.time = None
        if isinstance(info, str):  
            self._from_file(info)
        elif isinstance(info, dict):
            logging.debug('Perimeter.__init__ - processing perimeter')
            if 'poly' in info.keys():
                self.time = info.get('time', None)
                self.coords = info['poly']
                self.poly = self.polygonize()
            elif 'array' in info.keys():
                self.time = info.get('time', None)
                self.coords = self._from_array(info['array'])
                self.poly = self.polygonize()
            elif 'path' in info.keys():
                self._from_file(info['path'])
                self.time = info.get('time', self.time)   
        else:
            raise PerimeterError('Perimeter - info provided {} not recognized'.format(info))
        self.area = poly_area(self.poly)
        
    def plot(self):
        plot_perim(self, show=True)

    def _polygonize(self):
        return coords_to_polys(self.coords)

    def _from_file(self, path):
        logging.debug('Perimeter._from_file - processing perimeter {}'.format(path))
        self.path = path
        sinfo = self.path.split('.')
        if len(sinfo) > 1:
            ext = sinfo[-1]
        else:
            raise PerimeterError('Perimeter - string provided {} without extension'.format(path))
        if ext == 'pkl':
            self._from_pkl(path)
        elif ext == 'kml':
            self._from_kml(path)
        else:
            raise PerimeterError('Perimeter - string provided {} not recognized'.format(path))
    
    def _from_kml(self, path):
        self.time = None
        self.coords = read_kml(path)
        self.poly = self._polygonize()

    def _from_pkl(self, path):
        with open(path,'rb') as f:
            data = pickle.load(f)
        if isinstance(data,tuple) and len(data) == 2:
            self.time, self.coords = data
        else:
            self.coords = data
        self.poly = self._polygonize()
    
    def _from_array(self, array):
        logging.error('Not implemented.')
        pass

class Perimeters(object):
    """
    Define and process perimeters
    """
    def __init__(self, info):
        self._perims = []
        self._index = -1
        if not isinstance(info, Iterable) or isinstance(info,str):
            self.num_perims = 1
            logging.info('Perimeters.__info__ - processing {} perimeters'.format(self.num_perims))
            self._perims.append(Perimeter(info))
        else:
            self.num_perims = len(info)
            logging.info('Perimeters.__info__ - processing {} perimeters'.format(self.num_perims))
            for it in info:
                self._perims.append(Perimeter(it))   
        self.path = [p.path for p in self._perims]
        self.time = [p.time for p in self._perims]
        self.area = [p.area for p in self._perims]

    def __iter__(self):
        return self
    
    def __next__(self):
        self._index += 1
        if self._index >= len(self._perims):
            self._index = -1
            raise StopIteration
        else:
            return self._perims[self._index]

    def __getitem__(self, i):
        return self._perims[i]

    def __reversed__(self):
        return self._perims[::-1]
    
    def __len__(self):
        return self.num_perims
    
    def plot(self):
        for p in self._perims[:-1]:
            plot_perim(p)
        plot_perim(self._perims[-1], show=True)
    
    def to_csv(self, path):
        logging.info('Perimeters.to_csv - creating CSV file {}'.format(path))
        df = pd.DataFrame({'file': self.path, 'time': self.time, 'area': self.area})
        if not all([_ == None for _ in self.time]):
            df = df.groupby('time').agg({'file': 'first', 'area': 'mean'}).reset_index()
        df.to_csv(path, index=False)

    def to_kml(self, path):
        logging.info('Perimeters.to_kml - creating KML file {}'.format(path))
        kml = simplekml.Kml()
        kml.document.name = "Perimeters"
        for n,perim in enumerate(self._perims):
            if perim.time is not None:
                time = perim.time.strftime('%Y-%m-%dT%H:%M:%SZ')
                multipoly = kml.newmultigeometry(name='PERIM_'+time)
            else:
                multipoly = kml.newmultigeometry(name='PERIM_'+n)
            if isinstance(perim.poly,Iterable):
                rings = perim.poly.geoms
            else:
                rings = [perim.poly]
            for p in rings:
                if p.exterior is None:
                    logging.warning('Perimeters.to_kml - perimeter {} contains a polygon without an outer ring'.format(perim.path))
                    continue
                try:
                    p = p.simplify(1e-6)
                except:
                    logging.debug('Perimeters.to_kml - not possible to apply buffer to simplify multi-polygon')
                outer = np.c_[p.exterior.xy]
                if len(p.interiors) > 0:
                    inner = [np.c_[i.xy] for i in p.interiors]
                    multipoly.newpolygon(
                          outerboundaryis=outer, 
                          innerboundaryis=inner)
                else:
                    multipoly.newpolygon(
                          outerboundaryis=outer)             
            multipoly.timestamp.when = time
            multipoly.style.polystyle.color = '990000ff'
            multipoly.style.linestyle.color = simplekml.Color.red
        kml.save(path)
            
 
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')    
    import os.path as osp
    import glob
    dir_path = '/home/afarguell/scratch/forecasts/caldor/fire_init_data'
    files = sorted(glob.glob(osp.join(dir_path,'arcgis_perims_*')))
    perims = Perimeters(files)
    perims.to_csv('tests/caldor_perims.csv')
    perims.to_kml('tests/caldor_perims.kml')
    