import logging
import numpy as np
from packaging.version import parse as vparse
from scipy.spatial import Delaunay
from shapely.ops import transform
from shapely.geometry import Polygon, MultiPolygon
from matplotlib.path import Path
import pyproj
if vparse(pyproj.__version__) < vparse('2.2'):
    from functools import partial
    proj_merc = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
    proj_latlon = "+proj=latlong +ellps=WGS84 +datum=WGS84 +units=degree +no_defs"
else:
    proj_merc = "epsg:3857"
    proj_latlon = "epsg:4326"

def is_clockwise(ring):
    """
    Check if ring is clock-wise.
    :param ring: list of coordinates
    :return: if the order is clock-wise or not
    """
    xc = [c[0] for c in ring]
    yc = [c[1] for c in ring]
    xx = np.array(xc[1:])-np.array(xc[:-1])
    yy = np.array(yc[1:])+np.array(yc[:-1])
    return sum(xx*yy) > 0.0

def validate(poly):
    """
    Correct self-intersection in a polygon.
    :param poly: shapely polygon object
    :return: same polygon without self intersections
    """
    if not poly.is_valid:
        poly = poly.buffer(0)
    return poly

def coords_to_polys(coords):
    """
    Polygonize list coordinates format.
    :param coords: list of lists with outer and inner coordinates
    :return: shapely MultiPolygon object
    """
    mp = Polygon()
    for coord in coords:
        outer = Polygon(coord[0]).buffer(1e-6)
        inners = validate(MultiPolygon([Polygon(p).buffer(1e-6) for p in coord[1:]]))
        n_poly = validate(outer.difference(inners))
        mp = mp.union(n_poly)
    return mp

def polys_to_coords(polys):
    """
    Convert shapely polygons object to coords list format.
    :param polys: shapely MultiPolygon object
    :return: list of lists with outer and inner coordinates 
    """
    if isinstance(polys, Polygon):
        polys = [polys]
    elif isinstance(polys, MultiPolygon):
        pass
    else:
        logging.error(f'polys_to_coords - format of {polys} not recognized')
        return []
    coords = []
    for poly in polys:
        outer = [[x,y] for x,y in np.c_[poly.exterior.xy]]
        coord = [outer]
        for interior in poly.interiors:
            inner = [[x,y] for x,y in np.c_[interior.xy]]
            coord.append(inner)
        coords.append(coord)
    return coords

def lonlat_to_merc(*args):
    """
    Tranform coordinates or shapely object from WGS84 degrees to mercator meters.
    When coordinates, pass two arguments: lonlat_to_merc(lons,lats).
    :param args: two numbers, two lists, or a shapely object in WGS84 degrees
    :return: two numbers, two lists, or a shapely object transformed in mercator meters
    """
    if vparse(pyproj.__version__) < vparse('2.2'):
        transformer = partial(pyproj.transform, pyproj.Proj(proj_latlon), pyproj.Proj(proj_merc))
        if len(args) != 2:
            return transform(transformer, *args)  
        else:
            return transformer(args[0], args[1])
    else:
        transformer = pyproj.Transformer.from_crs(proj_latlon, proj_merc).transform
        if len(args) != 2:
            return transform(transformer, *args)
        else:
            return transformer(args[1], args[0])

def merc_to_lonlat(*args):
    """
    Tranform coordinates from mercator meters to WGS84 degrees.
    :param args: two numbers, two lists, or a shapely object in mercator meters
    :return: two numbers, two lists, or a shapely object transformed in WGS84 degrees 
    """
    if vparse(pyproj.__version__) < vparse('2.2'):
        transformer = partial(pyproj.transform, pyproj.Proj(proj_merc), pyproj.Proj(proj_latlon))
        return transformer(*args)
    else:
        transformer = pyproj.Transformer.from_crs(proj_merc, proj_latlon).transform
        if len(args) != 2:
            return transformer(*args)
        else:
            return transformer(*args)[::-1]

def featureset_to_coords(feature_set):
    """
    Feature set to perimeters.
    Feature sets need to have a SHAPE or shape field (geocoded) with a MultiPolygon object.
    Returning perimeter is in format:
    [[outer, inners], [outer, inners], ...]
    :param feature_set: feature set with a SHAPE or shape MultiPolygon.
    :return: return list of lists with outer and inner coordinates 
    """
    coord = []
    coords = []
    for rings in np.array(feature_set.get('SHAPE',feature_set.get('shape', []))):
        for poly in rings['rings']:
            poly = np.array(poly)
            if is_clockwise(poly):
                if len(coord):
                    coords.append(coord)
                coord = [poly]
            else:
                coord.append(poly)
        if len(coord):
            coords.append(coord)
    return coords

def simplify_coords(coords, min_outer_coords=3, max_outer_coords=1e10, min_inner_coords=3, max_inner_coords=1e10):
    """
    Simplify perimeter by setting min and max coordinates depending on outer or inner ring.
    If ring is smaller than minimal coordinates, the ring is removed.
    If ring is larger than maximal coordinates, the ring is thinned.
    :param coords: input coords to simplify
    :param min_outer_coords: min outer coordinates
    :param max_outer_coords: max outer coordinates
    :param min_inner_coords: min inner coordinates
    :param max_inner_coords: max inner coordinates
    :return: return coords simplified
    """
    result = []
    for coord in coords:
        outer = coord[0]
        inners = coord[1:]
        if len(outer) > min_outer_coords:
            if len(outer) > max_outer_coords:
                thinning = int(np.ceil(len(outer)/max_outer_coords))
                outer = outer[::thinning]
            pp = [outer]
            for inner in inners:
                if len(inner) > min_inner_coords:
                    if len(inner) > max_inner_coords:
                        thinning = int(np.ceil(len(inner)/max_inner_coords))
                        inner = inner[::thinning]
                    pp.append(inner)
            result.append(pp)
    return result
        
def create_coords(points, rings):
    """
    Create list of lists format with outer and inner coordinates 
    from points and rings.
    :param points: np.array of shape (n,2) points
    :param rings: closed rings created from the edges
    :return: list of lists format with outer and inner coordinates
    """
    pp = [Polygon(points[r]).buffer(0) for r in rings]
    polys = []
    for p in pp:
        if isinstance(p,Polygon):
            polys.append(p)
        elif isinstance(p,MultiPolygon):
            for _ in p:
                polys.append(_)
    indxs = [[] for _ in range(len(polys))]
    for i,p1 in enumerate(polys):
        for j,p2 in enumerate(polys):
            if i != j and p1.contains(p2):
                indxs[i].append(j)
    inner = [i for indx in indxs for i in indx]
    mp = []
    for i,indx in enumerate(indxs):
        if len(indx):
            mp.append([np.c_[polys[i].exterior.xy]] + [np.c_[polys[_].exterior.xy] for _ in indx])
        elif not i in inner:
            exterior = polys[i].exterior.coords
            if len(exterior):
                mp.append([np.c_[exterior.xy]])
    return mp
    
def edges_to_rings(edges):
    """
    Rings creation from pairs of edges
    :param edges: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array
    :return: closed rings created from the edges
    """
    edges_list = list(edges)
    rings = []
    while len(edges_list):
        edge = edges_list.pop(0)
        ring = list(edge)
        next_edge = [i for i, e in enumerate(edges_list) if e[0] == edge[1]]
        while len(next_edge):
            edge = edges_list.pop(next_edge[0])
            ring.append(edge[1])
            next_edge = [i for i, e in enumerate(edges_list) if e[0] == edge[1]]
        rings.append(ring)
    return rings

def alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    The set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    Then, the perimeters using list of lists format with outer and inner 
    coordinates is returned.
    :param points: np.array of shape (n,2) points
    :param alpha: alpha value
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges
    :return: list of lists format with outer and inner coordinates
    """
    assert points.shape[0] > 3, "Need at least four points"
    def add_edge(edges, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it is not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))
    logging.info(f'points -> edges, using alpha shape with alpha={alpha}')
    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    logging.info('edges -> linear rings')
    rings = edges_to_rings(edges)
    logging.info('rings -> coords')
    coords = create_coords(points, rings)
    return coords

def transform_polys(polys, transformer, min_outer_coords=3, max_outer_coords=1e10, min_inner_coords=3, max_inner_coords=1e10):
    """
    Transform polys
    :param polys:
    :param tgt_proj:
    :param src_proj:
    :param min_outer_coords: 
    :param max_outer_coords: 
    :param min_inner_coords: 
    :param max_inner_coords:
    :return: 
    """
    transf_polys = []
    for poly in polys:
        transf_poly = []
        p = poly[0]
        if len(p) > min_outer_coords: 
            if len(p) > max_outer_coords:
                thinning = int(np.ceil(len(p)/max_outer_coords))
                p = p[::thinning]
            mp = Polygon(p)
            transf_mp = transformer(mp)
            transf_poly.append(np.c_[transf_mp.exterior.coords.xy])
            for p in poly[1:]:
                if len(p) > min_inner_coords: 
                    if len(p) > max_inner_coords:
                        thinning = int(np.ceil(len(p)/max_outer_coords))
                        p = p[::thinning]
                    mp = Polygon(p)
                    transf_mp = transformer(mp)
                    transf_poly.append(np.c_[transf_mp.exterior.coords.xy])
        if len(transf_poly):
            transf_polys.append(transf_poly)
    return transf_polys

def mask_perim(perims, points):
    """
    Mask perimeter
    :param perims:
    :param points:
    :return: 
    """
    # create False mask with dimension of the points
    mask = np.zeros(len(points),dtype=bool)
    # for each polygon
    for perim in perims:
        # create Path element with the outerboundary coordinates
        path = Path(perim[0],closed=True)
        # mask outer polygon
        mask_out = path.contains_points(points)
        # set mask to True inside outer polygon
        mask[mask_out] = True
        # for each innerboundary polygon
        for p in perim[1:]:
            # create Path element with the innerboundary coordinates
            path = Path(p,closed=True)
            # mask inner polygon only at points that are True (the rest are already False)
            mask_inn = path.contains_points(points)
            # set mask to False inside the inner polygon
            mask[mask_inn] = False
    return mask

def poly_area(poly):
    """
    Polygon area
    :param poly:
    :return: 
    """
    return lonlat_to_merc(poly).area/4047.

''' # set parameters to reduce complexity of multipolygons (from transform_poly)
    min_inner_coords = 10     # minimum number of coordinates for an inner polygon
    max_inner_coords = 1000   # maximum number of coordinates for an inner polygon
    min_outer_coords = 50     # minimum number of coordinates for an outer polygon
    max_outer_coords = 10000  # maximum number of coordinates for an outer polygon'''