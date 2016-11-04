import numpy as np
from netCDF4 import Dataset
import time

def tanktopo(y, y_base, y_break, y_coast,
             fluid_depth, z_bottom, z_break, z_wall):
    
    ''' This function generates the topographical profile of the continental
    slope and shelf without the canyon. The profile is created in parts using
    the equation of a line: topography = z2 = (m * y2) - (m* y1) + z1, where
    the values for y represent key distances along the cross-shore direction
    and the values for z2 are the calculated depths based on a known z1 depth.
    '''
    
    sls_ct = (z_wall - z_break) / (y_coast - y_break)
    sls_sb = (z_break - z_bottom) / (y_break - y_base)
    topo_sp = np.zeros(len(y))
    slope_profile = np.zeros(len(y))
   
    for jj in np.arange(len(y)):

        if y[jj] > y_base and y[jj] <= y_break:
            topo_sp[jj] = (sls_sb * y[jj]) - (sls_sb * y_base) + z_bottom
                    
        elif y[jj] > y_break and y[jj] < y_coast:
            topo_sp[jj] = (sls_ct * y[jj]) - (sls_ct * y_break) + z_break
                                  
        elif y[jj] >= y_coast:
            topo_sp[jj] = z_wall

        slope_profile[jj] = topo_sp[jj] - fluid_depth
        
    return slope_profile
#----------------------------------------------------------------------------------

def canyontopo(y, y_base, y_break, y_head, y_coast,
               fluid_depth, z_bottom, z_break, z_wall):
    
    ''' This function generates the topographical profile for the canyon along
    its axis (cross-shore direction). Similar to tanktopo, the profile is
    formed using a collection of lines.
    '''
    
    slc_L = (z_break - z_bottom) / (y_head - y_base)
    slc_ct = (z_wall - z_break) / (y_coast - y_head)
    topo_cp = np.zeros(len(y))
    canyon_profile = np.zeros(len(y))
    
    for ii in np.arange(len(y)):
        
        if y[ii] <= y_base:
            topo_cp[ii] = z_bottom
        
        elif y[ii] > y_base and y[ii] <= y_head:
            topo_cp[ii] = (slc_L * y[ii]) - (slc_L * y_base) + z_bottom
                    
        elif y[ii] > y_head and y[ii] <= y_coast :      
            topo_cp[ii] = (slc_ct * y[ii]) - (slc_ct * y_head) + z_break
          
        elif y[ii] > y_coast:
            topo_cp[ii] = z_wall
        
        canyon_profile[ii] = topo_cp[ii] - fluid_depth
 
    return canyon_profile
#----------------------------------------------------------------------------------
   
def widthprofile(y, y_base, y_break, y_head, y_coast, cR, L,
                 w_break, w_mid, w_head, p):
    
    ''' This function defines the width profile of the canyon (top-down view).
    The width of the canyon is defined for all distances in the cross-shore
    direction. 
    '''

    sigmaa = 1.0 / ((9e-7) * cR)
    half = -w_break / 2.0+ w_mid / 2.0
    e = (L / 2.0 - sigmaa * half**2) / half**p 
    sc = 1
    dG_dxh = p * e * (w_head - w_break / 2)**(p-1) + 2 * sigmaa * (w_head - w_break)
    dh = 0.5 / dG_dxh / sc
    Ah = (w_break - w_head) / (y_base - y_head)**2;
    wp = np.zeros(len(y))

    for l in np.arange(len(y)):

            if y[l] <= y_base:
                wp[l] = w_break

            elif y[l] > y_base and y[l] <= y_head:   
                wp[l] = Ah * (y[l] - y_head)**2 + dh * (y[l] - y_head) + w_head

            elif y[l] > y_head and y[l] <= y_coast:
                wp[l] = wp[l-1]
            
            elif y[l] > y_coast:
                wp[l] = w_break
                
        
    width_profile = wp 
    return width_profile
#----------------------------------------------------------------------------------

def make_topo_smooth(y, y_base, y_break, y_head, y_coast, cR, L,
                     x, x_wall, w_break, w_mid, w_head, p,
                     fluid_depth, z_bottom, z_break, z_wall):
    
    ''' This function returns the depth field of the continental slope and
    shelf with a sech-shaped canyon. It uses the functions tanktopo,
    canyontopo, and widthprofile.
    
    :arg y: Array of cross-shore distances
    :arg y_base: Distance to the base of the continental slope
    :arg y_break: Distance to the shelf break
    :arg y_head: Distance to the canyon head
    :arg y_coast: Distance beyond y_head where shelf flattens
    :arg cR: Radius of curvature at the shelf break depth
    :arg x: Array of alongshore distances
    :arg x_wall: Width of the Domain
    :arg w_break: Width of the canyon at the shelf break
    :arg w_mid: Width of the canyon half-way along its length
    :arg w_head: Width of the canyon head
    :arg p: Geometric parameter used to help shape of canyon
    :arg fluid_depth: Total height of the fluid in the domain.
    :arg z_bottom: Depth of the deep ocean (measured upward)
    :arg z_break: Depth of the shelf break (measured upward)
    :arg z_wall: Depth of shelf beyond y_coast (measured upward)
    '''
    
    # Topography without the canyon
    slope_profile = tanktopo(y, y_base, y_break, y_coast,
                             fluid_depth, z_bottom, z_break, z_wall)
    
    # Slope of the canyon
    canyon_profile = canyontopo(y, y_base, y_break, y_head, y_coast,
                                fluid_depth, z_bottom, z_break, z_wall)
  
    # Slope of the canyon as well as the shape
    width_profile = widthprofile(y, y_base, y_break, y_head, y_coast, cR, L,
                                 w_break, w_mid, w_head, p)
  
    # Depth of the canyon (negative values set to zero)
    canyondepth = slope_profile - canyon_profile
  
    canyondepth[canyondepth < 0] = 0
  
    # Sech shaped canyon
    topography = np.zeros((len(y),len(x)))
    for j in np.arange(len(x)):
        topography[:,j] = (slope_profile - canyondepth * 
                           (1.0 / (np.cosh(0.5 / width_profile * (x[j] - (0.5 * x_wall))))**50))
    #topography=np.transpose(topography)
    topo = -1* topography[0:-1, :]
   
    return topo
#----------------------------------------------------------------------------------

def set_domain_grid(xsize, ysize, x_wall, y_wall):
    
    ''' Sets up the domain dimensions and grid cells used
    to generate the canyon bathymetry.
    
    :arg xsize: X direction dimension (alongshore)
    :arg ysize: Y direction dimension (cross-shore)
    '''
    
    xgrd_all = np.arange(0, xsize, 1)
    xgrd_bounds = [0, xsize-1]
    xval_bounds = [0, x_wall]
    xval_all = np.interp(xgrd_all, xgrd_bounds, xval_bounds)
    x_edge = np.zeros(xsize)
    x_edge[:] = xval_all[:]
    x = ((x_edge[1:] + x_edge[0:-1])/2)

    ygrd_all = np.arange(0, ysize, 1)
    ygrd_bounds = [0, ysize-1]
    yval_bounds = [0, y_wall]
    yval_all = np.interp(ygrd_all, ygrd_bounds, yval_bounds)
    y_edge = np.zeros(ysize)
    y_edge[:] = yval_all[:]
    y = ((y_edge[1:] + y_edge[0:-1])/2)

    X, Y = np.meshgrid(x, y)
    
    return x, y, y_edge, X, Y
#----------------------------------------------------------------------------------

def create_bathy_file(X, Y, bathymetry, filename, title, description):
    
    """ This function creates a netCDF4 file for
    the canyon bathymetry given the filename and 
    the x and y grid cell number.
    
    :arg X: Alongshore indices (from set_domain_grid)
    :arg Y: Cross-shore indices (from set_domain_grid)
    :arg bathymetry: Canyon bathymetry (from make_topo_smooth)
    :arg filename: Directory and name of netcdf file
    :arg title: Title of bathymetry version
    :arg description: Details about bathymetry version
    """
    
    dataset = Dataset(filename, 'w')
    file_x = dataset.createDimension('x', X.shape[1])
    file_y = dataset.createDimension('y', X.shape[0])

    file_X = dataset.createVariable('X', 'f8', ('y','x'))
    file_Y = dataset.createVariable('Y', 'f8', ('y','x'))
    Bathymetry = dataset.createVariable('Bathymetry', 'f8', ('y','x'))

    dataset.title = title
    dataset.author = 'Idalia A. Machuca'
    dataset.institution = 'Dept of Earth, Ocean & Atmospheric Sciences, University of British Columbia'
    dataset.source = 'bitbucket.org/CanyonsUBC/mackenzie_canyon/bathymetry/notebooks/make_mackenzie.ipynb'
    dataset.description = description
    dataset.timeStamp = time.ctime(time.time())
    file_X.standard_name = 'Along-Shore Distance'
    file_X.units = 'm'
    file_Y.standard_name = 'Cross-Shore Distance'
    file_Y.units = 'm'
    Bathymetry.standard_name = 'Bathymetry'
    Bathymetry.units = 'm'
    Bathymetry.positive = 'upward'

    file_X[:] = X[:]
    file_Y[:] = Y[:]
    Bathymetry[:] = bathymetry[:]

    dataset.close()
    
#----------------------------------------------------------------------------------

def define_Mackenzie_measurements():
    ''' This function defines all measurements made
    for Mackenzie Canyon that are used to create 
    the idealized bathymetry profile.
    '''
    
    # Alongshore
    w_break = 62681.735776859277  
    w_mid = 46456.969337226466  
    w_head = 14142.13562373095 
    width_f = 62681.735776859277
    x_wall = width_f * 7

    # Adjustments
    mouth = 51865.209919559762
    length = 74607.305272339116
    y_wall_1300 = 174731.93755006552
    y_wall = (mouth + length) * 2.57
    adjust = y_wall - y_wall_1300

    # Cross-shore
    cR = 9246.0
    L = 93744.0     
    y_base = np.mean([16500.0,7000.0]) + adjust 
    y_break = np.mean([38000.0, 57500.0]) + adjust 
    y_coast = 148105.0 + adjust 
    y_head = y_break + L

    # Depths
    fluid_depth = 1300.0
    z_bottom = fluid_depth - fluid_depth
    z_break = fluid_depth - 80.0
    z_wall = fluid_depth - 40.0 
    p = 4.0
    
    return w_break, w_mid, w_head, width_f, x_wall,\
            mouth, length, y_wall_1300, y_wall, adjust,\
            cR, L, y_base, y_break, y_coast, y_head,\
            fluid_depth, z_bottom, z_break, z_wall, p

# --------------------------------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset
import time
from salishsea_tools import geo_tools

def step_forward(lat1, lon1, d, brng):
    ''' Given a start point, initial bearing, and distance, 
    this process will calculate the destination point
    and final bearing travelling along a great circle arc.
    
    φ2 = asin( sin φ1 ⋅ cos δ + cos φ1 ⋅ sin δ ⋅ cos θ )
    λ2 = λ1 + atan2( sin θ ⋅ sin δ ⋅ cos φ1, cos δ − sin φ1 ⋅ sin φ2 )
    
    φ is latitude
    λ is longitude 
    θ is the bearing (clockwise from north) 
    δ is the angular distance d/R
    d is the distance travelled
    R is the earth’s radius
    
    :arg lat1: Latitude of start point (deg)
    :arg lon1: Longitude of start point (deg)
    :arg d: Distance to destination point (km)
    :arg brng: Initial bearing (rad)
    :returns: Latitude and longitude of destination point (deg)
    '''
    
    R = 6378.1 

    lat1 = np.radians(lat1) 
    lon1 = np.radians(lon1)
    brng = np.radians(brng)

    lat2 = np.arcsin( np.sin(lat1)*np.cos(d/R) +
             np.cos(lat1)*np.sin(d/R)*np.cos(brng))

    lon2 = lon1 + np.arctan2(np.sin(brng)*np.sin(d/R)*np.cos(lat1),
             np.cos(d/R)-np.sin(lat1)*np.sin(lat2))

    lat2 = np.degrees(lat2)
    lon2 = np.degrees(lon2)

    return lat2, lon2

# ------------------------------------------------------------------------------------------------

def great_circle_points(lat1, lon1, lat2, lon2, npts):
    ''' The longitude and latitudes of intermediate points at 
    any fraction along the great circle path between two points 
    is calculated. The haversine in geo_tools calculates the 
    great circle distance between two points.
    
    A = sin((1−f)⋅δ) / sin δ
    B = sin(f⋅δ) / sin δ
    x = A ⋅ cos φ1 ⋅ cos λ1 + B ⋅ cos φ2 ⋅ cos λ2
    y = A ⋅ cos φ1 ⋅ sin λ1 + B ⋅ cos φ2 ⋅ sin λ2
    z = A ⋅ sin φ1 + B ⋅ sin φ2
    φi = atan2(z, √x² + y²)
    λi = atan2(y, x) 
    
    φ is latitude
    λ is longitude 
    f is fraction along great circle route
    (f=0 is point 1, f=1 is point 2)
    δ is the angular distance d/R between the two points.
    
    :arg lat1: Latitude of start point (deg)
    :arg lon1: Longitude of start point (deg)
    :arg lat2: Latitude of start point (deg)
    :arg lon2: Longitude of start point (deg)
    :arg npts: Number of intermediate points along the great circle
    :returns: Latitudes and longitudes of all intermediate points (deg)
    '''
    
    R = 6378.1 
    
    d = geo_tools.haversine(lon1, lat1, lon2, lat2) / R
    lats = np.zeros(npts)
    lons = np.zeros(npts)
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    for i in range(npts):
        f = float(i) / float(npts-1.)
        A = np.sin((1 - f) * d) / np.sin(d)
        B = np.sin(f * d) / np.sin(d)
        x = A * np.cos(lat1) * np.cos(lon1) +  B * np.cos(lat2) * np.cos(lon2)
        y = A * np.cos(lat1) * np.sin(lon1) +  B * np.cos(lat2) * np.sin(lon2)
        z = A * np.sin(lat1) +  B*np.sin(lat2)
        lats[i] = np.arctan2(z, np.sqrt(x*x + y*y))
        lons[i] = np.arctan2(y, x)

    lats = np.degrees(lats)
    lons = np.degrees(lons)

    return lats, lons

# ------------------------------------------------------------------------------------------------

def calculate_initial_compass_bearing(pointA, pointB):
    ''' The bearing between two points is calculated.
       
    θ = atan2( sin Δλ ⋅ cos φ2 , cos φ1 ⋅ sin φ2 − sin φ1 ⋅ cos φ2 ⋅ cos Δλ )
    φ1,λ1 is the start point
    φ2,λ2 the end point 
    
    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees
    :Returns:
      The bearing in degrees
    :Returns Type:
      float
    '''
    
    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = np.radians(pointA[0])
    lat2 = np.radians(pointB[0])
    
    diffLong = np.radians(pointB[1] - pointA[1])

    x = np.sin(diffLong) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (np.sin(lat1)
            * np.cos(lat2) * np.cos(diffLong))

    initial_bearing = np.arctan2(x, y)

    # Now we have the initial bearing but math.atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = np.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

# ------------------------------------------------------------------------------------------------

def create_grid_file(filename):
    ''' This function saves the coordinates for all
    points in the defined grid. This grid will then 
    be processed to calculate the scaling factors
    required by NEMO.
    
    :arg filename: Directory and name of netCDF4 file
    '''
    
    dataset = Dataset(filename, 'w')
    x = dataset.createDimension('x', nx)
    y = dataset.createDimension('y', ny)

    lons = dataset.createVariable('grid_lons', 'f8', ('x','y'))
    lats = dataset.createVariable('grid_lats', 'f8', ('x','y'))

    dataset.title = 'Mackenzie Canyon Coordinates Grid'
    dataset.author = 'Idalia A. Machuca'
    dataset.institution = 'Dept of Earth, Ocean & Atmospheric Sciences, University of British Columbia'
    dataset.source = 'bitbucket.org/CanyonsUBC/mackenzie_canyon/bathymetry/notebooks/'
    dataset.timeStamp = time.ctime(time.time())
    lons.standard_name = 'Longitude'
    lats.standard_name = 'Latitude'
    lons.units = 'degrees east'
    lons.units = 'degrees north'

    lons[:] = thelons[:]
    lats[:] = thelats[:]

    dataset.close()
    
# ------------------------------------------------------------------------------------------------

def compile_grid_coordinates(nx, ny, lon, lat):
    ''' Uses previously defined functions to calculate
    the longitudes and latitudes of all grid cells in
    the domain defined by its two southern corner
    points.
    
    :arg nx: Grid size (alongshore)
    :arg ny: Grid size (cross-shore)
    :arg lon: List of longitudes of corner points
    :arg lat: List of latitudes of corner points
    :returns: Longitudes and latitudes of all grid points
    '''

    lat2, lon2 = great_circle_points(lat[1], lon[1], lat[0], lon[0], nx)

    thelats = np.zeros((nx, ny))
    thelons = np.zeros_like(thelats)
    thelats[:,0] = lat2
    thelons[:,0] = lon2
    
    dx = geo_tools.haversine(thelons[-1,0], thelats[-1,0], thelons[-2,0], thelats[-2,0])

    angle = 0
    for j in range(1, ny):
        for i in range(0, nx-1):
            prevangle = angle
            bearing = calculate_initial_compass_bearing((thelats[i, j-1], thelons[i, j-1]), 
                                                        (thelats[i+1, j-1], thelons[i+1, j-1]))
            angle = bearing - 90.
            thelats[i,j], thelons[i,j] = step_forward(thelats[i,j-1], thelons[i,j-1], dx, angle)
        i = nx-1
        angle = 2*angle - prevangle
        thelats[i, j], thelons[i, j] = step_forward(thelats[i, j-1], thelons[i, j-1], dx, angle)

    return thelons, thelats

#----------------------------------------------------------------------------------
import numpy as np
from pyproj import Proj, transform

def transform_coords_geog_stere(lon_g, lat_g):
    ''' Convert geographic coordinates of all grid
    points in the domain to stereographic coordinates. 
    
    :arg lon_g: Geographic longitudes
    :arg lat_g: Geographic latitudes
    :returns: Stereographic longitudes and latitudes
    '''
    
    # Geographical coordinate system
    proj_geogr = Proj("+init=EPSG:4326")
    
    # IBCAO polar stereographic
    proj_stere = Proj("+init=EPSG:3996") 
    
    lon_s =  np.zeros_like(lon_g)
    lat_s =  np.zeros_like(lat_g)
    for i in np.arange(lon_g.shape[0]):
        lon_s[i,:], lat_s[i,:] = transform(proj_geogr, proj_stere, lon_g[i,:], lat_g[i,:])
    
    return lon_s, lat_s
