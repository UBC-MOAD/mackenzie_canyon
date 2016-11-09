# This module extracts and prepares the 
# realistic Mackenzie Canyon bathymetry.

# Uses coordinate grid from functions_grid.py

# Imports
import numpy as np
from scipy.interpolate import griddata

def extract_canyon(lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method):
    
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

    z_canyon = griddata(points, values, (lon_s_grid, lat_s_grid), method = interp_method)
    
    return z_canyon

# -----------------------------------------------------------------------------------------

def rotate_extracted_canyon(lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method):
    
    ''' Extracts the bathymetric subdomain and 
    flips and rotates it by 270 deg.
    
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
    
    z_rotated = extract_canyon(lon_s_grid_rotate, lat_s_grid_rotate, x_region, y_region, z_region, interp_method)
    
    return z_rotated

# -----------------------------------------------------------------------------------------

def extend_rotated_canyon(extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method):
    
    ''' Extracts the bathymetric subdomain, 
    flips and rotates it, and extrapolates so that the
    end column looks like the start. This is to facilitate
    periodic boundaries.
    
    All values are in IBCAO's Stereographic projection.
    Use cubic interpolation method for the extraction.
    
    :arg extension: Number of grid cells for extrapolation.
    :arg lon_s_grid: Longitudes of subdomain grid
    :arg lat_s_grid: Latitudes of subdomain grid
    :arg x_region: X values of larger dataset
    :arg y_region: Y values of larger dataset
    :arg z_region: Bathymetry of larger dataset
    :arg interp_method: 'nearest', linear', 'cubic'
    :returns: Subdomain bathymetry
    '''
    
    z_rotated = rotate_extracted_canyon(lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method)
    x_extended = np.arange(z_rotated.shape[1] + extension)
    x_original = np.append(np.arange(z_rotated.shape[1]), x_extended[-1])
    y_original = np.arange(z_rotated.shape[0])
    z_original = np.concatenate([z_rotated, z_rotated[:,:1]], axis=1)
    
    X_original, Y_original = np.meshgrid(x_original, y_original)
    points = np.c_[np.ravel(X_original), np.ravel(Y_original)]
    values = np.ravel(z_original)

    X_extended, Y_extended = np.meshgrid(x_extended, y_original)

    z_extended = griddata(points, values, (X_extended, Y_extended), method = interp_method)
    
    return z_extended

# -----------------------------------------------------------------------------------------

def flatten_extended_canyon(fluid_depth, extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method):
    
    ''' Removes bathymetry over 0 m and under
    fluid depth from the rotated and extended
    subdomain.
    
    All values are in IBCAO's Stereographic projection.
    Use cubic interpolation method for the extraction.
    
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
    
    z_extended = extend_rotated_canyon(extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method)

    z_fluid_depth = np.array(z_extended, copy=True)
    z_fluid_depth[z_fluid_depth < -fluid_depth] = -fluid_depth

    z_flattened = np.array(z_fluid_depth, copy=True)
    z_flattened[z_flattened > 0] = 0
    
    return z_flattened
    
# -----------------------------------------------------------------------------------------

def match_shape_canyon(fluid_depth, extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method):
    
    ''' The shape of the canyon bathymetry must be equal
    to that of the coordinates set. The coordinates set
    has 3 less rows and columns than the original. 
    In the x direction, the extension was reduced by 3.
    In the y direction, the last three rows were simply
    removed from the extended bathymetry. This will not
    remove any bathymetric features since that area is
    already flattened at 1300 m depth.
    
    All values are in IBCAO's Stereographic projection.
    Use cubic interpolation method for the extraction.
    
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
    
    z_flattened = flatten_extended_canyon(fluid_depth, extension - 3, lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method)
    
    z_shape = z_flattened[:-3,:]
    
    return z_shape
    
# -----------------------------------------------------------------------------------------    

def positive_canyon(fluid_depth, extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method):
    
    ''' Converts the extracted, rotated, extended,
    flattened, and shaped subdomain into positive 
    values for NEMO.
    
    All values are in IBCAO's Stereographic projection.
    Use cubic interpolation method for the extraction.
    
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
    
    z_shape = match_shape_canyon(fluid_depth, extension, lon_s_grid, lat_s_grid, x_region, y_region, z_region, interp_method)
    
    z_positive0 = z_shape * -1
    z_positive = np.array(z_positive0, copy=True)
    z_positive[z_positive == -0.] = 0.
    
    return z_positive

# -----------------------------------------------------------------------------------------