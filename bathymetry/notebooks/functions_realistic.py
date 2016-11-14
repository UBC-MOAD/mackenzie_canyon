# This module extracts and prepares the 
# realistic Mackenzie Canyon bathymetry.

# Uses coordinate grid from functions_grid.py

# Imports
import numpy as np
from scipy.interpolate import griddata
from salishsea_tools import bathy_tools

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
    
    X_region, Y_region = np.meshgrid(x_region, y_region) 
    points = np.c_[np.ravel(X_region), np.ravel(Y_region)]
    values = np.ravel(z_region)

    z_canyon = griddata(points, values, (lon_s_grid, lat_s_grid), method = 'cubic')
    
    return z_canyon

# -----------------------------------------------------------------------------------------

def cyclic_canyon(extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region):
    
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
    
    lon_s_grid_rotate = np.fliplr(np.rot90(lon_s_grid, 3))
    lat_s_grid_rotate = np.fliplr(np.rot90(lat_s_grid, 3))
    
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
    
    bathy = bathymetry[:-3,:]
    
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

def canyon_for_model(fluid_depth, extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region):
    
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
    
    z_cyclic = cyclic_canyon(extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region)

    z_cyclic[z_cyclic < -fluid_depth] = -fluid_depth
    z_cyclic[z_cyclic > 0] = 0
    
    z_shape = match_shape(z_cyclic)
    z_positive = make_positive(z_shape)
    
    return z_positive
    
# -----------------------------------------------------------------------------------------

def smooth_canyon(max_norm_depth_diff, smooth_factor, fluid_depth, extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region):
    
    z_positive = canyon_for_model(fluid_depth, extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region)
    
    z_to_smooth = np.array(z_positive, copy=True)
    z_masked0 = np.ma.array(z_to_smooth)
    z_smoothed = bathy_tools.smooth(z_masked0, max_norm_depth_diff, smooth_factor)

    z_original = np.ma.array(z_positive)

    return z_original, z_smoothed
    
# -----------------------------------------------------------------------------------------