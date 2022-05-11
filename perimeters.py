try:
    from .geometry import coords_to_polys, polys_to_coords, poly_area
    from .plot_tools import plot_perim,_colors
    from .tools import load_pkl,read_kml
except:
    from geometry import coords_to_polys, polys_to_coords, poly_area
    from plot_tools import plot_perim,_colors
    from tools import load_pkl,read_kml 
from collections.abc import Iterable
import pandas as pd
import numpy as np
import simplekml
import logging

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
            logging.debug('Perimeter.__init__ - processing perimeter from {}'.format(info))
            self._from_file(info)
        elif isinstance(info, dict):
            logging.debug('Perimeter.__init__ - processing perimeter from dictionary')
            if 'poly' in info.keys():
                self.time = info.get('time', None)
                self.coords = info['poly']
                self.poly = self._polygonize()
            elif 'array' in info.keys():
                self.time = info.get('time', None)
                self.coords = self._from_array(info['array'])
                self.poly = self._polygonize()
            elif 'path' in info.keys():
                self._from_file(info['path'])
                self.time = info.get('time', self.time)   
        else:
            raise PerimeterError('Perimeter - info provided {} not recognized'.format(info))
        self._update_params()
        
    def __len__(self):
        return len(self.coords)
    
    def plot(self, show=False, **args):
        return plot_perim(self, show=show, **args)

    def simplify(self, simplify_tol=1e-4):
        self.poly = self.poly.simplify(simplify_tol)
        self.coords = polys_to_coords(self.poly)
        self._update_params()

    def _update_params(self):
        self.area = poly_area(self.poly)
        self.bounds = self._bounds()

    def _bounds(self):
        x = [_[0] for coord in self.coords for c in coord for _ in c]
        y = [_[1] for coord in self.coords for c in coord for _ in c]
        if len(x) and len(y):
            return min(x),max(x),min(y),max(y)
        else:
            return None

    def _polygonize(self):
        return coords_to_polys(self.coords)

    def _from_file(self, path):
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
            raise PerimeterError('Perimeter - string provided {} extension not recognized'.format(path))
    
    def _from_kml(self, path):
        times,coords = read_kml(path)
        if len(times) > 1:
            self.time = None
            self.coords = [_ for c in coords for _ in c]
        elif len(times) == 1:
            self.time = times[0]
            self.coords = coords[0]
        else:
            self.time = None
            self.coords = []
        self.poly = self._polygonize()

    def _from_pkl(self, path):
        data = load_pkl(path)
        if isinstance(data,tuple) and len(data) == 2:
            self.time, self.coords = data
        else:
            self.coords = data
        self.poly = self._polygonize()
    
    def _from_array(self, array):
        raise NotImplementedError('Perimeter._from_array not implemented')

class Perimeters(object):
    """
    Define and process perimeters
    """
    def __init__(self, info):
        self._perims = []
        self._index = -1
        if not isinstance(info, Iterable) or isinstance(info, dict):
            self.num_perims = 1
            logging.info('Perimeters.__init__ - processing {} perimeters'.format(self.num_perims))
            self._perims.append(Perimeter(info))
        elif isinstance(info, str):
            self._from_file(info)
        else:
            self.num_perims = len(info)
            logging.info('Perimeters.__init__ - processing {} perimeters'.format(self.num_perims))
            for n,it in enumerate(info):
                logging.debug('Perimeters.__init__ - processing perimeter {}/{}'.format(n+1, self.num_perims))
                if isinstance(it, Perimeter):
                    self._perims.append(it)   
                else:
                    self._perims.append(Perimeter(it))   
        self.sort()

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
    
    def plot(self, show=False, **args):
        i = -1
        for i,p in enumerate(self._perims[:-1]):
            if 'color' not in args:
                color = _colors(i % 10) 
                args['color'] = color
            plot_perim(p, show=False, **args)
        if 'color' not in args:
            color = _colors((i+1) % 10)
            args['color'] = color
        plot_perim(self._perims[-1], show=show, **args)
    
    def sort(self):
        valid_perims = [p for p in self._perims if p.time is not None]
        invalid_perims = [p for p in self._perims if p.time is None]
        self._perims = sorted(valid_perims, key=lambda x: x.time) + invalid_perims
        self._update_params()
    
    def _update_params(self):
        self.path = [p.path for p in self._perims]
        self.time = [p.time for p in self._perims]
        self.area = [p.area for p in self._perims]
        self.bounds = self._bounds()
    
    def _bounds(self):
        x_min_opts = [perim.bounds[0] for perim in self._perims if perim.bounds is not None]
        x_max_opts = [perim.bounds[1] for perim in self._perims if perim.bounds is not None]
        y_min_opts = [perim.bounds[2] for perim in self._perims if perim.bounds is not None]
        y_max_opts = [perim.bounds[3] for perim in self._perims if perim.bounds is not None]
        if len(x_min_opts)+len(x_max_opts)+len(y_min_opts)+len(y_max_opts) == 0:
            return None
        x_min = min(x_min_opts)
        x_max = max(x_max_opts)
        y_min = min(y_min_opts)
        y_max = max(y_max_opts)
        return (x_min,x_max,y_min,y_max)
    
    def _from_file(self, path):
        logging.debug('Perimeters._from_file - processing perimeters {}'.format(path))
        self.path = path
        sinfo = self.path.split('.')
        if len(sinfo) > 1:
            ext = sinfo[-1]
        else:
            raise PerimeterError('Perimeters - string provided {} without extension'.format(path))
        if ext == 'pkl':
            self._from_pkl(path)
        elif ext == 'kml':
            self._from_kml(path)
        else:
            raise PerimeterError('Perimeters - string provided {} extension not recognized'.format(path))
    
    def _from_kml(self, path):
        times,coords = read_kml(path)
        self.num_perims = len(times)
        logging.info('Perimeters._from_kml - processing {} perimeters'.format(self.num_perims))
        for n,(t,c) in enumerate(zip(times,coords)):
            logging.debug('Perimeters._from_kml - processing perimeter {}/{}'.format(n+1, self.num_perims))
            self._perims.append(Perimeter({'time': t, 'poly': c}))  

    def _from_pkl(self, path):
        logging.error('Not implemented.')
        pass
    
    def _from_array(self, array):
        logging.error('Not implemented.')
        pass
    
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
    