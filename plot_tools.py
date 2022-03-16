import matplotlib.pyplot as plt
from process_perimeter_masks import read_wrfinfo

def plot_bbox(path):
    _,_,_,_,fxlon,fxlat = read_wrfinfo(path)
    bbox = (fxlon.min(),fxlon.max(),fxlat.min(),fxlat.max())
    xbox = [bbox[0],bbox[0],bbox[1],bbox[1],bbox[0]]
    ybox = [bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]]
    plt.plot(xbox,ybox,'k-')

def plot_perim(perim, show=False, **kargs):
    gkargs = kargs.copy()
    for p in perim.coords:
        if len(p):
            pp = plt.plot(p[0][:,0],p[0][:,1],**gkargs)
            kargs['color'] = pp[-1].get_color()
            if len(p) > 1:
                for ring in p[1:]:
                    plt.plot(ring[:,0],ring[:,1],**kargs)
    if show:
        plt.show()

if __name__ == '__main__':
    from perimeters import Perimeter
    import sys
    wrf_path = './wrf_init'
    ir_path = sys.argv[1]
    perims_path = sys.argv[2]
    ir_perim = Perimeter(sys.argv[1])
    past_perim = Perimeter(sys.argv[1])
    plot_perim(ir_perim,'g')
    plot_perim(ir_perim,'b')
    plt.show()
