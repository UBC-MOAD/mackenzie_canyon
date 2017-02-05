import  numpy as np
import netCDF4 as nc
import  matplotlib.cm as cm
import  os
import glob
import scipy as sc, scipy.io
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

#--------------------------------------------------------------------------------------

def get_variables(projection):
    ''' Loads the file for either the Stereographic
    or Geographic projection of the Arctic Ocean
    (IBCAO) and returns the file's data for x, y,
    and z. This data can be plotted using the
    functions below.

    Stereographic projection IBCAO_V3_500m_RR.
    Geographic projection IBCAO_V3_30arcsec_RR.
    '''
    if projection == 'S':
        ibcao_grid_name = 'IBCAO_V3_500m_RR.grd'
    elif projection == 'G':
        ibcao_grid_name = 'IBCAO_V3_30arcsec_RR.grd'
    ibcao_grid_dir = '/ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/grid'
    ibcao_grid = os.path.join(ibcao_grid_dir, ibcao_grid_name)
    ibcao_nc = scipy.io.netcdf_file (ibcao_grid)

    x = ibcao_nc.variables['x'][:]
    y = ibcao_nc.variables['y'][:]
    z = ibcao_nc.variables['z'][:]
    return x, y, z

#--------------------------------------------------------------------------------------

def def_regionG(xl=-145, xr=-133, yb=68.6, yt=72.5):
    ''' Returns an extract of the complete IBCAO bathymetric
    grid using xl, xr, yb, and yt, which correpond to the left,
    right, bottom, and top boundaries when plotted.
    '''

    ibcao_file = scipy.io.netcdf_file('/ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/grid/IBCAO_V3_30arcsec_RR.grd')
    x = ibcao_file.variables['x'][:]
    y = ibcao_file.variables['y'][:]
    z = ibcao_file.variables['z'][:]
    xmin = np.where(np.round(x,2)==xl)[0][0]
    xmax = np.where(np.round(x,2)==xr)[0][0]
    ymin = np.where(np.round(y,2)==yb)[0][0]
    ymax = np.where(np.round(y,2)==yt)[0][0]
    x_region = x[xmin:xmax]
    y_region = y[ymin:ymax]
    z_region = z[ymin:ymax, xmin:xmax]
    return x_region, y_region, z_region

#--------------------------------------------------------------------------------------

def def_regionS(xl=-1750e3, xr=-1000e3, yb=1300e3, yt=2050e3):
    ''' Returns an extract of the complete IBCAO bathymetric
    grid using xl, xr, yb, and yt, which correpond to the left,
    right, bottom, and top boundaries when plotted.
    '''

    ibcao_file = scipy.io.netcdf_file('/ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/grid/IBCAO_V3_500m_RR.grd')
    x = ibcao_file.variables['x'][:]
    y = ibcao_file.variables['y'][:]
    z = ibcao_file.variables['z'][:]
    xmin = np.where(x==xl)[0][0]
    xmax = np.where(x==xr)[0][0]
    ymin = np.where(y==yb)[0][0]
    ymax = np.where(y==yt)[0][0]
    x_region = x[xmin:xmax]
    y_region = y[ymin:ymax]
    z_region = z[ymin:ymax, xmin:xmax]
    return x_region, y_region, z_region

#--------------------------------------------------------------------------------------

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

#--------------------------------------------------------------------------------------

def plot_region(fig, ax, x_region, y_region, z_region):
    ax.contour(x_region, y_region, z_region, 25, colors='k', linestyles='solid', alpha=0.6)
    ax.contour(x_region, y_region, z_region, levels = [-80, -40.1], colors='k', linestyles='solid', alpha=0.6)
    return fig, ax

#--------------------------------------------------------------------------------------

#def mask_output(var):
    #var_m = np.ma.masked_values(var, 0.0)
    #var_m = np.ma.masked_values(var, np.isnan(var))
    #return var_m

#--------------------------------------------------------------------------------------

def load_model_output(path, cfg):
    ''' This function is used to load the important
    variables stored in the different output files
    produced by a NEMO run. The path ends at the CONFIG
    directory while the cfg path starts at the specific
    configuration directory and ends in the final directory
    without a slash.
    '''
    gridT = nc.Dataset(glob.glob(path + cfg + '/GYRE_*_grid_T.nc')[0])
    gridU = nc.Dataset(glob.glob(path + cfg + '/GYRE_*_grid_U.nc')[0])
    gridV = nc.Dataset(glob.glob(path + cfg + '/GYRE_*_grid_V.nc')[0])
    gridW = nc.Dataset(glob.glob(path + cfg + '/GYRE_*_grid_W.nc')[0])
    mesh_mask = nc.Dataset(glob.glob(path + cfg + '/mesh_mask.nc')[0])

    lon = gridT.variables['nav_lon']
    lat = gridT.variables['nav_lat']
    tem = gridT.variables['votemper'][:]
    sal = gridT.variables['vosaline'][:]
    ssh = gridT.variables['sossheig'][:]
    U = gridU.variables['vozocrtx'][:]
    V = gridV.variables['vomecrty'][:]
    W = gridW.variables['vovecrtz'][:]
    
    tmask0 = 1 - mesh_mask['tmask'][:]
    time_len = tem.shape[0]
    tmask = np.tile(tmask0, (time_len, 1, 1, 1))
    
    tem_masked = np.ma.array(tem, mask=tmask)
    sal_masked = np.ma.array(sal, mask=tmask)
    ssh_masked = np.ma.array(ssh, mask=tmask)
    U_masked = np.ma.array(U, mask=tmask)
    V_masked = np.ma.array(V, mask=tmask)
    W_masked = np.ma.array(W, mask=tmask)

    return gridT, lon, lat, tem_masked, sal_masked, ssh_masked, U_masked, V_masked, W_masked
