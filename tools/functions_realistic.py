# This module extracts and prepares the 
# realistic Mackenzie Canyon bathymetry.

# Uses coordinate grid from functions_grid.py

# Imports
import numpy as np
import time
from netCDF4 import Dataset
from scipy.interpolate import griddata
from salishsea_tools import bathy_tools
from scipy import interpolate
from scipy.interpolate import Rbf

def extract_canyon(lon_s_grid, lat_s_grid, x_region, y_region, z_region):
    
    ''' Extracts a bathymetry domain defined by a grid
    of coordinates from a larger bathymetry dataset.
    
    All values are in IBCAO's Stereographic projection.
    Use cubic interpolation method for the extraction.
    
    :arg lon_s_grid: Longitudes of subdomain grid
    :arg lat_s_grid: Latitudes of subdomain grid
    :arg x_region: X values of larger dataset
    :arg y_region: Y values of larger dataset
    :arg z_region: Bathymetry of larger dataset
    :arg interp_method: 'nearest', linear', 'cubic'
    :returns: Subdomain bathymetry
    '''
    
    lon_s_grid_rotate = np.fliplr(np.rot90(lon_s_grid, 3))
    lat_s_grid_rotate = np.fliplr(np.rot90(lat_s_grid, 3))

    X_region, Y_region = np.meshgrid(x_region, y_region) 
    points = np.c_[np.ravel(X_region), np.ravel(Y_region)]
    values = np.ravel(z_region)

    z_canyon = griddata(points, values, (lon_s_grid_rotate, lat_s_grid_rotate), method = 'cubic')
    
    return z_canyon

# -----------------------------------------------------------------------------------------

def cyclic_canyon_extend(extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region):
    
    ''' Extracts the bathymetric subdomain, 
    flips and rotates it by 270 deg, and extrapolates 
    so that the end column looks like the start. 
    This is to facilitate periodic boundaries.
    
    :arg extension: Number of grid cells for extrapolation.
    :arg lon_s_grid: Longitudes of subdomain grid
    :arg lat_s_grid: Latitudes of subdomain grid
    :arg x_region: X values of larger dataset
    :arg y_region: Y values of larger dataset
    :arg z_region: Bathymetry of larger dataset
    :arg interp_method: 'nearest', linear', 'cubic'
    :returns: Subdomain bathymetry
    '''
    
    z_rotated = extract_canyon(lon_s_grid_rotate, lat_s_grid_rotate, x_region, y_region, z_region)
    
    x_extended = np.arange(z_rotated.shape[1] + extension)
    x_original = np.append(np.arange(z_rotated.shape[1]), x_extended[-1])
    y_original = np.arange(z_rotated.shape[0])
    z_original = np.concatenate([z_rotated, z_rotated[:,:1]], axis=1)
    
    X_original, Y_original = np.meshgrid(x_original, y_original)
    points = np.c_[np.ravel(X_original), np.ravel(Y_original)]
    values = np.ravel(z_original)

    X_extended, Y_extended = np.meshgrid(x_extended, y_original)

    z_cyclic = griddata(points, values, (X_extended, Y_extended), method = 'linear')
    
    return z_cyclic

# -----------------------------------------------------------------------------------------

def cyclic_canyon_truncate(index, lon_s_grid, lat_s_grid, x_region, y_region, z_region):
    
    ''' Extracts the bathymetric subdomain, 
    flips and rotates it by 270 deg, and extrapolates 
    so that the end column looks like the start. 
    This is to facilitate periodic boundaries.
    
    :arg extension: Number of grid cells for extrapolation.
    :arg lon_s_grid: Longitudes of subdomain grid
    :arg lat_s_grid: Latitudes of subdomain grid
    :arg x_region: X values of larger dataset
    :arg y_region: Y values of larger dataset
    :arg z_region: Bathymetry of larger dataset
    :arg interp_method: 'nearest', linear', 'cubic'
    :returns: Subdomain bathymetry
    '''
    
    z_rotated = extract_canyon(lon_s_grid, lat_s_grid, x_region, y_region, z_region)

    x_original = np.zeros(z_rotated.shape[1]-index+1)
    x_original[1:] = np.arange(index, z_rotated.shape[1])
    x_extended = np.arange(z_rotated.shape[1])
    y_original = np.arange(z_rotated.shape[0])
    z_original = np.concatenate([z_rotated[:,-1:], z_rotated[:,index:]], axis=1)

    X_original, Y_original = np.meshgrid(x_original, y_original)
    points = np.c_[np.ravel(X_original), np.ravel(Y_original)]
    values = np.ravel(z_original)

    X_extended, Y_extended = np.meshgrid(x_extended, y_original)

    z_cyclic = griddata(points, values, (X_extended, Y_extended), method = 'linear')
    
    return z_cyclic

# -----------------------------------------------------------------------------------------

def match_shape(bathymetry):
    '''The shape of the canyon bathymetry must be equal
    to that of the coordinates set. The coordinates set
    has 3 less rows and columns than the original. 
    In the x direction, the extension was reduced by 3.
    In the y direction, the last three rows were simply
    removed from the extended bathymetry. This will not
    remove any bathymetric features since that area is
    already flattened at 1300 m depth.
    '''
    
    bathy = bathymetry[:-3,:-3]
    
    return bathy
    
# -----------------------------------------------------------------------------------------    

def make_positive(bathymetry):
    
    ''' Makes the canyon bathymetry positive.
    This is required for the NEMO model.
    '''
    
    bathy0 = bathymetry * -1
    bathy = np.array(bathy0, copy=True)
    bathy[bathy == -0.] = 0.
    
    return bathy

# -----------------------------------------------------------------------------------------

def canyon_for_model(fluid_depth, index, lon_s_grid, lat_s_grid, x_region, y_region, z_region):
    
    ''' Limits bathymetry values from 0 m to the
    fluid depth, fixes the bathymetry array shape, 
    and converts values to positive. This is applied
    to the rotated, extended, cyclic bathymetry.
    
    :arg fluid_depth: Maximum subdomain depth (positive)
    :arg extension: Number of grid cells for extrapolation.
    :arg lon_s_grid: Longitudes of subdomain grid
    :arg lat_s_grid: Latitudes of subdomain grid
    :arg x_region: X values of larger dataset
    :arg y_region: Y values of larger dataset
    :arg z_region: Bathymetry of larger dataset
    :arg interp_method: 'nearest', linear', 'cubic'
    :returns: Subdomain bathymetry
    '''
    
    #z_cyclic = cyclic_canyon(extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region)
    z_cyclic = cyclic_canyon_truncate(index, lon_s_grid, lat_s_grid, x_region, y_region, z_region)

    z_cyclic[z_cyclic < -fluid_depth] = -fluid_depth
    z_cyclic[z_cyclic > 0] = 0
    
    z_shape = match_shape(z_cyclic)
    z_positive = make_positive(z_shape)
    
    return z_positive
    
# -----------------------------------------------------------------------------------------

def smooth_canyon(max_norm_depth_diff, smooth_factor, z_positive):
    
    #z_positive = canyon_for_model(fluid_depth, index, lon_s_grid, lat_s_grid, x_region, y_region, z_region)
    
    z_to_smooth = np.array(z_positive, copy=True)
    z_masked0 = np.ma.array(z_to_smooth)

    z_smoothed = bathy_tools.smooth(z_masked0, max_norm_depth_diff, smooth_factor)

    z_original = np.ma.array(z_positive)

    return z_original, z_smoothed

# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------

def calculate_slopes(z_smoothed):
    
    slopes = np.zeros([z_smoothed.shape[-2], z_smoothed.shape[-1] - 1])

    for j in range(z_smoothed.shape[-2]):
        for i in range(z_smoothed.shape[-1] - 1):
            slopes[j, i] = z_smoothed[j, i+1] - z_smoothed[j, i]

    slopes_row_avg = np.mean(slopes, axis=1)
    slopes_row_max = np.max(slopes, axis=1)
    slopes_row_min = np.min(slopes, axis=1)

    slopes_maxmin = np.zeros_like(slopes)

    for j in range(z_smoothed.shape[-2]):
        for i in range(z_smoothed.shape[-1] - 1):
            if slopes[j, i] == slopes_row_max[j]:
                slopes_maxmin[j, i] = 1
            if slopes[j, i] == slopes_row_min[j]:
                slopes_maxmin[j, i] = -1
   
    return slopes, slopes_row_max, slopes_row_min, slopes_maxmin, slopes_row_avg
            
# ----------------------------------------------------------------------------------------

def identify_interp_inds(z_smoothed, slopes, slopes_row_max, slopes_row_min, buffer):
    
    bathy_centre = z_smoothed[:, int(0.5 * z_smoothed.shape[-1])]
    frac_full = np.where(bathy_centre == 1300)[0][0] / z_smoothed.shape[-2]
    frac_y = frac_full * buffer
    ind_y = int(frac_y * z_smoothed.shape[-2])

    ind_x_left0 = np.zeros(ind_y)
    ind_x_right0 = np.zeros(ind_y)
    ind_x_diff = np.zeros(ind_y)
    for j in range(ind_y):
        ind_x_left0[j] = np.where(slopes[j, :] == slopes_row_max[j])[0][0]
        ind_x_right0[j] = np.where(slopes[j, :] == slopes_row_min[j])[0][0]
        ind_x_diff[j] = ind_x_right0[j] - ind_x_left0[j]
    
    ind_x_diff_max = np.where(ind_x_diff == ind_x_diff.max())[0][0]
    random_add = int(ind_x_diff[ind_x_diff_max]/7)
    ind_x_left = ind_x_left0[ind_x_diff_max] + random_add
    ind_x_right = ind_x_right0[ind_x_diff_max]

    print(ind_y, ind_x_left, ind_x_right)

    ind_x_left = 125
    ind_x_right = 190
    ind_y = 73

    print(ind_y, ind_x_left, ind_x_right)

    return ind_y, ind_x_left, ind_x_right

# ----------------------------------------------------------------------------------------
    
def gather_interp_args(z_smoothed, ind_y, ind_x_left, ind_x_right):
    
    x_orig = np.array([])
    y_orig = np.array([])
    
    for j in range(ind_y):     
        x_orig_left = np.arange(ind_x_left + 1)
        x_orig_right = np.arange(ind_x_right, z_smoothed.shape[-1])
        x_orig0 = np.concatenate((x_orig_left, x_orig_right))
        y_orig0 = np.ones_like(x_orig0) * j

        x_orig = np.append(x_orig, x_orig0)
        y_orig = np.append(y_orig, y_orig0)

    for j in range(ind_y, z_smoothed.shape[-2]):
        x_orig0 = np.arange(z_smoothed.shape[-1])
        y_orig0 = np.ones_like(x_orig0) * j

        x_orig = np.append(x_orig, x_orig0)
        y_orig = np.append(y_orig, y_orig0)

    z_orig = np.zeros_like(x_orig)
    for k in range(len(x_orig)):
        z_orig[k] = z_smoothed[y_orig[k], x_orig[k]]
        
    return x_orig, y_orig, z_orig
    
# ----------------------------------------------------------------------------------------

def perform_interp(z_smoothed, x_orig, y_orig, z_orig):

    rbf = Rbf(x_orig, y_orig, z_orig, function='linear')

    xi = np.arange(z_smoothed.shape[-1])
    yi = np.arange(z_smoothed.shape[-2])
    XI, YI = np.meshgrid(xi, yi)
    z_nocanyon_real = rbf(XI, YI)
    
    return z_nocanyon_real

# ----------------------------------------------------------------------------------------

def make_no_canyon_real(z_smoothed, buffer=0.78):
    slopes, slopes_row_max, slopes_row_min, slopes_maxmin, slopes_row_avg = calculate_slopes(z_smoothed)
    
    ind_y, ind_x_left, ind_x_right = identify_interp_inds(z_smoothed, slopes, slopes_row_max, slopes_row_min, buffer)
    
    x_orig, y_orig, z_orig = gather_interp_args(z_smoothed, ind_y, ind_x_left, ind_x_right)
    
    z_nocanyon_real = perform_interp(z_smoothed, x_orig, y_orig, z_orig)
    
    return z_nocanyon_real
    


# -----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------

def create_bathy_file(bathymetry, filename, title, description, ipynbname):
    
    ''' This function creates a netCDF4 file for
    the canyon bathymetry given the filename and 
    the x and y grid cell number.
    
    :arg X: Alongshore indices (from set_domain_grid)
    :arg Y: Cross-shore indices (from set_domain_grid)
    :arg bathymetry: Canyon bathymetry (from make_topo_smooth)
    :arg filename: Directory and name of netcdf file
    :arg title: Title of bathymetry version
    :arg description: Details about bathymetry version
    :arg ipynbname: Name of source ipython notebook
    '''
    
    directory = '/ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/NEMO_files/realistic/'
    dataset = Dataset(directory + filename, 'w')
    file_x = dataset.createDimension('x', bathymetry.shape[1])
    file_y = dataset.createDimension('y', bathymetry.shape[0])

    Bathymetry = dataset.createVariable('Bathymetry', 'f8', ('y','x'))

    dataset.title = title
    dataset.author = 'Idalia A. Machuca'
    dataset.institution = 'Dept of Earth, Ocean & Atmospheric Sciences, University of British Columbia'
    dataset.source = 'bitbucket.org/CanyonsUBC/mackenzie_canyon/bathymetry/notebooks/' + ipynbname
    dataset.description = description
    dataset.timeStamp = time.ctime(time.time())
    Bathymetry.standard_name = 'Bathymetry'
    Bathymetry.units = 'm'
    Bathymetry.positive = 'upward'

    Bathymetry[:] = bathymetry[:]

    dataset.close()
    
# -----------------------------------------------------------------------------------------
