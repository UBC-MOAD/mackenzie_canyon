# This module creates a grid of coordinates
# for the canyon domain.

# Imports
import numpy as np
from netCDF4 import Dataset
import time
from salishsea_tools import geo_tools
from pyproj import Proj, transform

def find_distance(p_one, p_two):
    ''' Finds the distance between two points.
    This is used to find the various canyon
    dimensions, such as widths and cross and along
    shore distances.
    '''
    xmax = abs(p_one[0])
    xmin = abs(p_two[0])
    ymax = p_one[1]
    ymin = p_two[1]
    dist = np.sqrt((xmax - xmin)**2 + (ymax - ymin)**2)
    return dist

# ------------------------------------------------------------------------------------------------

def match_distance(p1_x, p1_y, p2_x, m_slope, ideal):
    ''' Use a reference point [p1_x, p1_y] and the 
    slope of the line to find the second point
    [p2_x, p2_y]. dist gives the length of the line. 
    The goal is to realistic distances found here
    with those set in the idealized canyon domain.
    The ideal case would be for dist = dim leading
    to diff_dist ~ 0.
    '''
    p2_y = p1_y + (m_slope * (p2_x - p1_x))
    dist = find_distance([p1_x, p1_y], [p2_x, p2_y])
    diff = ideal - dist
    return p2_y, dist, diff

# ------------------------------------------------------------------------------------------------

def match_lines(p1_x, p1_y, p2_x_iterate, m_slope, ideal):
    '''This function searches for the stereographic coordinate
    which would create the line of the length required in 
    order to match with the corresponding side of the 
    idealized domain. p1_x and p1_y are the coordinates
    for the the point p1 where the line begines. m_slope is
    the slope of the line, p2_x_iterate is the coordinates
    through which the function iterates to find the one where
    the resulting line length is closes to the the value of 
    ideal which is the dimension of the idealized domain.
    '''
    for x_i in p2_x_iterate:
        y_i, dist_i, diff_i = match_distance(p1_x, p1_y, x_i, m_slope, ideal)
        if abs(diff_i) < 300.0:
            p2_x = x_i
            p2_y = y_i
            p2_dist = dist_i
            p2_diff = diff_i
        else:
            pass
    return p2_x, p2_y, p2_dist, p2_diff

# ------------------------------------------------------------------------------------------------

def match_domain(x_wall, y_wall, search_x, ax, slope=1.1, p_x0 = -1457500.0, p_y0 = 1348000.0):
    '''This is the main function used to create the realistic domain.
    It uses match_lines to create each section of the rectangle.
    It makes the right angle corners and the sides are as close
    as possible to the size of the idealized domain. The key here
    is to define the slope of the domain (the angle of the sides),
    the length and width of the idealized domain, and the initial
    reference point p2_x0 and p2_y0. They were called p2 because
    they were originally the second point in the line running
    along the axis of the canyon. This function also plots the
    lines. Search_x is the range of x coordinates that will be
    iterated through in order to find the one that produces the
    line length closest to the ideal size. The -500 is used
    when iterating leftward and 500 for rightward.'''

    # All the slopes needed to make the perpendicular corners
    m_slope = np.zeros(4)
    m_slope[0] = slope; m_slope[1] = -1/slope; m_slope[2] = slope; m_slope[3] = -1/slope

    # The values of the sides of the idealized domain
    ideal = np.zeros(4)
    ideal[0] = x_wall; ideal[1] = y_wall; ideal[2] = x_wall; ideal[3] = y_wall

    # Bottom Side
    lw = 1.5; 
    n=0; p1_x = p_x0; p1_y = p_y0 
    p2_xBR, p2_yBR, distBR, diffBR = match_lines(p1_x, p1_y, np.arange(p1_x, search_x[n], 500), m_slope[n], ideal[n])
    ax.plot([p1_x, p2_xBR], [p1_y, p2_yBR], 'k', linewidth=lw) 
    
    # Right Side
    n=1; p1_x = p2_xBR; p1_y = p2_yBR
    p2_xTR, p2_yTR, distTR, diffTR = match_lines(p1_x, p1_y, np.arange(p1_x, search_x[n], -500), m_slope[n], ideal[n])
    ax.plot([p1_x, p2_xTR], [p1_y, p2_yTR], 'k', linewidth=lw) 
    
    # Top Side
    n=2; p1_x = p2_xTR; p1_y = p2_yTR 
    p2_xTL, p2_yTL, distTL, diffTL = match_lines(p1_x, p1_y, np.arange(p1_x, search_x[n], -500), m_slope[n], ideal[n])
    ax.plot([p1_x, p2_xTL], [p1_y, p2_yTL], 'k', linewidth=lw) 

    # Left Side
    n=3; p1_x = p2_xTL; p1_y = p2_yTL
    p2_xBL, p2_yBL, distBL, diffBL = match_lines(p1_x, p1_y, np.arange(p1_x, search_x[n], 500), m_slope[n], ideal[n])
    ax.plot([p1_x, p2_xBL], [p1_y, p2_yBL], 'k', linewidth=lw)
    
    p2_BR = [round(p2_xBR/500.0)*500.0, round(p2_yBR/500.0)*500.0]
    p2_TR = [round(p2_xTR/500.0)*500.0, round(p2_yTR/500.0)*500.0]
    p2_TL = [round(p2_xTL/500.0)*500.0, round(p2_yTL/500.0)*500.0]
    p2_BL = [round(p2_xBL/500.0)*500.0, round(p2_yBL/500.0)*500.0]
    
    corner_lons = []
    corner_lats = []
    corner_all = [p2_BR, p2_TR, p2_TL, p2_BL]
    for i in np.arange(len(corner_all)):
        corner_lons.append(corner_all[i][0])
        corner_lats.append(corner_all[i][1])
    corner_lons, corner_lats
    
    return corner_lons, corner_lats, ax

# ------------------------------------------------------------------------------------------------

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

def create_grid(nx, ny, lon, lat):
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

def find_dx(lon_g, lat_g):
    ''' Finds the average distance between 
    two consecutive points in the grid of
    coordinates created by the functions above.
    
    :arg lon_g: Geographic longitudes
    :arg lat_g: Geographic latitudes
    :returns: Distance between points
    '''
    
    dx_start = geo_tools.haversine(lon_g[0,0], lat_g[0,0], lon_g[0,1], lat_g[0,1])
    dx_end = geo_tools.haversine(lon_g[-1,-2], lat_g[-1,-2], lon_g[-1,-1], lat_g[-1,-1])
    dx = np.mean([dx_start, dx_end])
    
    return dx
# ------------------------------------------------------------------------------------------------

def transform_coords(lon_orig, lat_orig, transformation):
    
    ''' Convert geographic coordinates of all grid
    points in the domain to stereographic coordinates. 
    
    :arg lon: Original longitudes
    :arg lat: Original latitudes
    :returns: Transformed longitudes and latitudes
    '''
    
    # Geographical coordinate system
    proj_geogr = Proj("+init=EPSG:4326")
    
    # IBCAO polar stereographic
    proj_stere = Proj("+init=EPSG:3996")
    
    lon_tran =  np.zeros_like(lon_orig)
    lat_tran =  np.zeros_like(lat_orig)
    for i in np.arange(lon_orig.shape[0]):
        if transformation == 'GS':
            lon_tran[i,:], lat_tran[i,:] = transform(proj_geogr, proj_stere, lon_orig[i,:], lat_orig[i,:])
        if transformation == 'SG':
            lon_tran[i,:], lat_tran[i,:] = transform(proj_stere, proj_geogr, lon_orig[i,:], lat_orig[i,:])
    
    return lon_tran, lat_tran

# ------------------------------------------------------------------------------------------------
