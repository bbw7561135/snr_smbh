def plot_single(in_file, zfunc, zname, out_file):

    import pylab
    import numpy
    import h5py
    from matplotlib.collections import PolyCollection
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    print(out_file)

    with h5py.File(in_file,'r+') as f:
        vert_idx_list = numpy.concatenate(([0],
                                           numpy.cumsum(f['Number of vertices in cell'])))
        verts = []
        for i in range(len(f['density'])):
            lowbound = int(vert_idx_list[i])
            upbound = int(vert_idx_list[i+1])
            verts.append([[x,y] for x,y
                          in zip(f['x position of vertices'][lowbound:upbound],
                                 f['y position of vertices'][lowbound:upbound])])
        coll = PolyCollection(verts, 
                              array=zfunc(f),
                              cmap = mpl.cm.jet,
                              edgecolors = 'none')
        fig, ax = plt.subplots()
        fig.suptitle(zname+' @ t = '+str(numpy.array(f['time'])[0]))
        ax.add_collection(coll)
        ax.autoscale_view()
        ax.set_aspect('equal')
        fig.colorbar(coll,ax=ax)
        print(out_file)
        if out_file==None:
            plt.show()
        else:
            plt.savefig(out_file)

def plot_all(zfunc, zname):

    import glob
    import numpy
    import joblib

    flist = glob.glob('snapshot_*.h5')

    joblib.Parallel(n_jobs=8)(joblib.delayed(plot_single)
                              (fname,
                               zfunc,
                               zname,
                               fname.replace('snapshot',zname).replace('.h5','.png')) for fname in flist)
    #[plot_single(fname,zfunc,zname,
    #             fname.replace('snapshot',zname).replace('.h5','.png'))
    # for fname in flist]

def log10_number_density_cgs(f):

    import numpy

    parsec2cm = 3.08e18
    solar_mass2g = 1.9891e33
    density_scale = solar_mass2g/parsec2cm**3
    proton_mass_g = 1.67e-24
    return numpy.log10(numpy.array(f['density'])*density_scale/proton_mass_g)

def log10_temperature(f):

    import numpy

    return numpy.log10(f['temperature'])

def x_velocity(f):

    import numpy

    return numpy.array(f['x_velocity'])

def y_velocity(f):

    import numpy

    return numpy.array(f['y_velocity'])

def main():

    import numpy

    #plot_single('snapshot_405.h5',
    #            log10_temperature,
    #            'log10_temperature',
    #            'log10_temperature_405.png')

    plot_all(log10_number_density_cgs, 'log10_number_density')
    plot_all(log10_temperature, 'log10_temperature')
    plot_all(x_velocity, 'x_velocity')
    plot_all(y_velocity, 'y_velocity')

if __name__ == '__main__':

    main()
