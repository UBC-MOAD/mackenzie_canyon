# This module creates a grid of coordinates
# for the canyon domain.

# Imports
import numpy as np
from netCDF4 import Dataset
import time
from salishsea_tools import geo_tools
from pyproj import Proj, transform

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

# ------------------------------------------------------------------------------------------------
 