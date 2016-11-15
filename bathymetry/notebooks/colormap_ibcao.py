import  numpy as np
import  matplotlib.cm as cm
import scipy as sc, scipy.io
import matplotlib.pyplot as plt

def Colormap():
    COLORMAP = """\
    # downloaded from IBCAO homepage
    #Discrete color table for Ocean and continous for land in RGB for the Arctic bathymetry and topography
    -6000	18	10	59	-5000	18	10	59
    -5000	22	44	103	-4000	22	44	103
    -4000	22	88	135	-3000	22	88	135
    -3000	22	138	170	-2000	22	138	170
    -2000	22	154	184	-1500	22	154	184
    -1500	23	170	198	-1000	23	170	198
    -1000	23	186	212	-500	23	186	212
    -500	24	196	223	-250	24	196	223
    -250	25	206	234	-100	25	206	234
    -100	27	216	245	-75	27	216	245
    -75	38	223	241	-50	38	223	241
    -50	49	230	236	-25	49	230	236
    -25	105	242	233	-10	105	242	233
    -10	161	255	230	0	161	255	230
    0	40	158	38	25	44	176	42
    25	44	176	42	50	49	195	46
    50	49	195	46	75	145	208	80
    75	145	208	80	100	242	202	90
    100	242	202	90	200	227	170	48
    200	227	170	48	300	190	140	40
    300	190	140	40	400	151	109	31
    400	151	109	31	500	114	80	23
    500	114	80	23	600	95	63	12
    600	95	63	12	700	81	57	16
    700	81	57	16	800	114	97	71
    800	114	97	71	1000	105	105	105
    1000	105	105	105	1500	170	170	170
    1500	170	170	170	5000	200	200	200
    """
    cmap = np.empty ((0,4))
    c = 0
    for l in COLORMAP.split("\n"):
      l = l.strip()
      if len(l) == 0 or l[0] == '#':
        continue
      ls = np.array([float (v) for v in l.split ()])
      if ls.shape[0] < 8:
        continue
      c += 1
      cmap.resize (c, 4)
      cmap[c-1,:] = ls[:4]
    c += 1
    cmap.resize (c, 4)
    cmap[c-1,:] = ls[4:]
    cmap[:,[1, 2, 3]] = cmap[:,[1, 2, 3]] / 255.
    cmap_out = cm.colors.ListedColormap (cmap[:,1:4], 'ibcao', c)
    norm     = cm.colors.BoundaryNorm (cmap[:,0], c)
    return (cmap_out, norm)

def def_region():
    ibcao_file = scipy.io.netcdf_file('/ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/grid/IBCAO_V3_500m_RR.grd')
    x_ibcao = ibcao_file.variables['x'][:]
    y_ibcao = ibcao_file.variables['y'][:]
    z_ibcao = ibcao_file.variables['z'][:]
    xl=-1750e3; xr=-1000e3; yb=1300e3; yt=2050e3
    xmin = np.where(x_ibcao==xl)[0][0]
    xmax = np.where(x_ibcao==xr)[0][0]
    ymin = np.where(y_ibcao==yb)[0][0]
    ymax = np.where(y_ibcao==yt)[0][0]
    x_region = x_ibcao[xmin:xmax]
    y_region = y_ibcao[ymin:ymax]
    z_region = z_ibcao[ymin:ymax, xmin:xmax]
    return x_region, y_region, z_region
    
def plot_region(fig, ax):
    ax.contour(x_region, y_region, z_region, 25, colors='k', linestyles='solid', alpha=0.6)
    ax.contour(x_region, y_region, z_region, levels = [-80, -40.1], colors='k', linestyles='solid', alpha=0.6)
    return fig, ax
