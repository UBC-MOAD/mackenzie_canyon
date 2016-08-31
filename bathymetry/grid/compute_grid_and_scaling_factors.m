% NAME: compute_grid_and_scaling_factors.m
%
% AUTHOR: J.-P. Paquin 
%
% DATE: July 2013
%
% REVISIONS: August 2016 - Idalia Machuca
%
% DESCRIPTION: Creates a coordinates.nc file for NEMO from outputs of the 
%              seagrid package (WHOI) available at 
%                 http://woodshole.er.usgs.gov/operations/modeling/seagrid/
%              Computes the latitudes and longitudes of u, v, and f grids
%              as simple averages 
%              and their respective scaling factors (based on ROMS code)
%
% HYPOTHESES:  Seagrid output file has been generated
%
% INPUTS: Seagrid output file 
%
% OUTPUT: Coordinates.nc file in compatible format for NEMO
%
% CALLED PGM & SCRIPTS: earthdist.m : compute distances
%
% NOTES: The final grid size is 2 grid points less in each direction to
%        avoid NaNs to be present in teh scaling factors. 
%
%--------------------------------------------------------------------------
clear all
clc
display('***** compute_grid_and_scaling_factors.m *****')

% Refer to salishsea_setup.ipynb for more information.
addpath /ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/grid/mexcdf_all/mexnc
addpath /ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/grid/mexcdf_all/netcdf_toolbox/netcdf
addpath /ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/grid/mexcdf_all/netcdf_toolbox/netcdf/nctype
addpath /ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/grid/mexcdf_all/netcdf_toolbox/netcdf/ncutility

%--- Outfile
fileout=('test_coordinates_seagrid_WestCoast2.nc')

%--- LOAD SEAGRID FILE 
%load(['seagrid_west_coast_1km_900x400_rot_new.mat']) % IM
load(['seagrid_west_coast_100x100_testAnchor.mat']) % IM
lon_T=s.geographic_grids{1,1}(:,:);
lat_T=s.geographic_grids{1,2}(:,:);

[dimx,dimy]=size(lon_T);

lat_u=NaN(dimx,dimy)    ;   lon_u=NaN(dimx,dimy);
lat_v=NaN(dimx,dimy)    ;   lon_v=NaN(dimx,dimy);
lat_f=NaN(dimx,dimy)    ;   lon_f=NaN(dimx,dimy);

display('Compute u, v points')
for jj=1:dimy-1
for ii=1:dimx-1

   lat_u(ii,jj) = 0.5 * ( lat_T(ii,jj) + lat_T(ii+1,jj) ) ;
   lon_u(ii,jj) = 0.5 * ( lon_T(ii,jj) + lon_T(ii+1,jj) ) ;

   lat_v(ii,jj) = 0.5 * ( lat_T(ii,jj) + lat_T(ii,jj+1) ) ;
   lon_v(ii,jj) = 0.5 * ( lon_T(ii,jj) + lon_T(ii,jj+1) ) ;

end
end

display('Compute f points')
for jj=1:dimy-1
for ii=1:dimx-1
   lat_f(ii,jj) = 0.5 * ( lat_u(ii,jj) + lat_u(ii  ,jj+1) ) ;
   lon_f(ii,jj) = 0.5 * ( lon_v(ii,jj) + lon_v(ii+1,jj  ) ) ;  
end
end

display('Put all lats/lons in 3D array')
lats=NaN(dimx,dimy,4);
lons=NaN(dimx,dimy,4);
for vv=1:4 % T,u,v,f in that order
  switch vv
      case 1 
        lats(:,:,vv)=lat_T(:,:);
        lons(:,:,vv)=lon_T(:,:);
      case 2
        lats(:,:,vv)=lat_u(:,:);
        lons(:,:,vv)=lon_u(:,:);        
      case 3
        lats(:,:,vv)=lat_v(:,:);
        lons(:,:,vv)=lon_v(:,:);
      case 4 
        lats(:,:,vv)=lat_f(:,:);
        lons(:,:,vv)=lon_f(:,:);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       COMPUTE THE SCALING FACTORS FOR T, u, v and f                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Compute distances between point for each grid... T,u,v,f')
radius = 6371*1000;

dist_lat=NaN(dimx,dimy,4);
dist_lon=NaN(dimx,dimy,4);
for vv=1:4 % T,u,v,f
for jj=1:dimy-1
for ii=1:dimx-1   
  % distance in lat  
  alat=lats(ii  ,jj  ,vv) ;
  alon=lons(ii  ,jj  ,vv) ;
  
  blat=lats(ii+1,jj,vv) ;
  blon=lons(ii+1,jj,vv) ;
  dist_lat(ii,jj,vv) = earthdist(alon, alat, blon, blat, radius); % eq to gy

  % distance in lon
  blat=lats(ii  ,jj+1,vv) ; 
  blon=lons(ii  ,jj+1,vv) ;
  dist_lon(ii,jj,vv) = earthdist(alon, alat, blon, blat, radius); % eq to gx

end
end
end

display('Compute scaling factors (inspired by seagrid2roms.m)')
sx=NaN(dimx,dimy,4);
sy=NaN(dimx,dimy,4);
for vv=1:4 % T,u,v,f
  sx(1:dimx-1,1:dimy  ,vv) = 0.5*(dist_lon(1:end-1,   :   ,vv) + dist_lon(2:end,   :  ,vv));
  sy(1:dimx  ,1:dimy-1,vv) = 0.5*(dist_lat(   :   ,1:end-1,vv) + dist_lat(  :  , 2:end,vv));
end
% sx and sy cannot be Inf, even if on land, so if values
% are Inf, set to an arbitrary non-zero value
sx(isinf(sx))=1.e+20;
sy(isinf(sy))=1.e+20;
sx(isnan(sx))=1.e+20;
sy(isnan(sy))=1.e+20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                     Write coordinates file NC file                    %%%%%
f=netcdf(fileout,'clobber');

xx = length(1:dimx-3) ;
yy = length(1:dimy-3) ;

f('x')=xx;
f('y')=yy;
f('time')='UNLIMITED';

f{'nav_lon'}=ncfloat('y','x');
f{'nav_lon'}.units='degrees_east';
f{'nav_lon'}.comment='at t points';

f{'nav_lat'}=ncfloat('y','x');
f{'nav_lat'}.units='degrees_north';
f{'nav_lat'}.comment='at t points';

f{'time'}=ncfloat('time');
f{'time'}.units='seconds since 0001-01-01 00:00:00';
f{'time'}.time_origin='0000-JAN-01 00:00:00';
f{'time'}.calendar='gregorian';

f{'time_steps'}=ncint('time');
f{'time_steps'}.units='seconds since 0001-01-01 00:00:00';
f{'time_steps'}.time_origin='0000-JAN-01 00:00:00';

f{'glamt'}=ncdouble('time', 'y', 'x');
f{'glamt'}.missing_value= 1.e+20;

f{'glamu'}=ncdouble('time', 'y', 'x');
f{'glamu'}.missing_value= 1.e+20;

f{'glamv'}=ncdouble('time', 'y', 'x');
f{'glamv'}.missing_value= 1.e+20;

f{'glamf'}=ncdouble('time', 'y', 'x');
f{'glamf'}.missing_value= 1.e+20;

f{'gphit'}=ncdouble('time', 'y', 'x');
f{'gphit'}.missing_value= 1.e+20;

f{'gphiu'}=ncdouble('time', 'y', 'x');
f{'gphiu'}.missing_value= 1.e+20;

f{'gphiv'}=ncdouble('time', 'y', 'x');
f{'gphiv'}.missing_value= 1.e+20;

f{'gphif'}=ncdouble('time', 'y', 'x');
f{'gphif'}.missing_value= 1.e+20;

f{'e1t'}=ncdouble('time', 'y', 'x');
f{'e1t'}.missing_value= 1.e+20;

f{'e1u'}=ncdouble('time', 'y', 'x');
f{'e1u'}.missing_value= 1.e+20;

f{'e1v'}=ncdouble('time', 'y', 'x');
f{'e1v'}.missing_value= 1.e+20;

f{'e1f'}=ncdouble('time', 'y', 'x');
f{'e1f'}.missing_value= 1.e+20;

f{'e2t'}=ncdouble('time', 'y', 'x');
f{'e2t'}.missing_value= 1.e+20;

f{'e2u'}=ncdouble('time', 'y', 'x');
f{'e2u'}.missing_value= 1.e+20;

f{'e2v'}=ncdouble('time', 'y', 'x');
f{'e2v'}.missing_value= 1.e+20;

f{'e2f'}=ncdouble('time', 'y', 'x');
f{'e2f'}.missing_value= 1.e+20;

f{'nav_lon'}(:,:) = permute(lons(1:dimx-3,1:dimy-3,1),[2,1]);
f{'nav_lat'}(:,:) = permute(lats(1:dimx-3,1:dimy-3,1),[2,1]);

f{'time'}(1:1)=0;
f{'time_step'}(1:1)= 0;

f{'glamt'}(1:1,:,:)= permute(lons(1:dimx-3,1:dimy-3,1)  ,[2,1]);
f{'glamu'}(1:1,:,:)= permute(lons(1:dimx-3,1:dimy-3,2)  ,[2,1]);
f{'glamv'}(1:1,:,:)= permute(lons(1:dimx-3,1:dimy-3,3)  ,[2,1]);
f{'glamf'}(1:1,:,:)= permute(lons(1:dimx-3,1:dimy-3,4)  ,[2,1]);
f{'gphit'}(1:1,:,:)= permute(lats(1:dimx-3,1:dimy-3,1)  ,[2,1]);
f{'gphiu'}(1:1,:,:)= permute(lats(1:dimx-3,1:dimy-3,2)  ,[2,1]);
f{'gphiv'}(1:1,:,:)= permute(lats(1:dimx-3,1:dimy-3,3)  ,[2,1]);
f{'gphif'}(1:1,:,:)= permute(lats(1:dimx-3,1:dimy-3,4)  ,[2,1]);
f{'e1t'}(1:1,:,:)= permute(dist_lon(1:dimx-3,1:dimy-3,1),[2,1]); 
f{'e1u'}(1:1,:,:)= permute(dist_lon(1:dimx-3,1:dimy-3,2),[2,1]);
f{'e1v'}(1:1,:,:)= permute(dist_lon(1:dimx-3,1:dimy-3,3),[2,1]);
f{'e1f'}(1:1,:,:)= permute(dist_lon(1:dimx-3,1:dimy-3,4),[2,1]);
f{'e2t'}(1:1,:,:)= permute(dist_lat(1:dimx-3,1:dimy-3,1),[2,1]);
f{'e2u'}(1:1,:,:)= permute(dist_lat(1:dimx-3,1:dimy-3,2),[2,1]);
f{'e2v'}(1:1,:,:)= permute(dist_lat(1:dimx-3,1:dimy-3,3),[2,1]);
f{'e2f'}(1:1,:,:)= permute(dist_lat(1:dimx-3,1:dimy-3,4),[2,1]);

close(f)
display('  END WRITING COORDINATES')
