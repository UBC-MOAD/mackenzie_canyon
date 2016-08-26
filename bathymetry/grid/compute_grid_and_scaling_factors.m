% NAME: compute_grid_and_scaling_factors.m
%
% AUTHOR: J.-P. Paquin 
%
% DATE: July 2013
%
% REVISIONS: 
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
display('***** compute_grid_and_scaling_factors.m *****')

addpath /users/staff/jppaquin/matlab/mexnc
addpath /users/staff/jppaquin/matlab/netcdf_toolbox/netcdf
addpath /users/staff/jppaquin/matlab/netcdf_toolbox/netcdf/nctype
addpath /users/staff/jppaquin/matlab/netcdf_toolbox/netcdf/ncutility

addpath ../seagrid

path_seagrid='/users/staff/jppaquin/NEMO_PREPARATION/1_Seagrid_generator';
%- Infile
file_seagrid='grid_mat/seagrid_west_coast_1km_900x400_rot_new.mat';
%file_seagrid='grid_mat/seagrid_west_coast_100x100_testAnchor.mat';

%- Outfile
fileout=([path_seagrid '/scaling_factors/coordinates_seagrid_WestCoast.nc']);


%--- LOAD SEAGRID FILE 
load([path_seagrid '/' file_seagrid])
lon_T=s.geographic_grids{1,1}(:,:);
lat_T=s.geographic_grids{1,2}(:,:);

[dimx,dimy]=size(lon_T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           COMPUTE LAT/LON FOR THE U, V AND f GRIDS                    %%%%%
%
% Arakawa-C grid
%
%           f---v---f---v---f  
%           |       |       |
%           u   T   u   T   u
%           |       |       |
%           f---v---f---v---f
%
%
% We define the latitudes longitudes in the output of seagrid to be the
% location of the T grid. Therefore the location of the u grid is defined
% as:
%       lat_u(i,j) = 0.5 * ( lat_T(i,j) + lat_T(i+1,j) )
%       lon_u(i,j) = 0.5 * ( lon_T(i,j) + lon_T(i+1,j) ) 
%
% similarly for the v grid
%       lat_v(i,j) = 0.5 * ( lat_T(i,j) + lat_T(i,j+1) )
%       lon_v(i,j) = 0.5 * ( lon_T(i,j) + lon_T(i,j+1) ) 
%
% and for the f grid:
%       lat_f(i,j) = 0.5 * ( lat_u(i,j) + lat_u(i,j+1) )
%       lon_f(i,j) = 0.5 * ( lon_v(i,j) + lon_v(i+1,j) ) 
%

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

% construire un tableau contenant toutes les lat-lons
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
%clear lat_u lon_u lat_v lon_v lat_f lon_f 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       COMPUTE THE SCALING FACTORS FOR T, u, v and f                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Compute distances between point for each grid... T,u,v,f')
%  earthdist -- Distance in meters between two lon/lats.
%    earthdist(alon, aloat, blon, blat, radius) returns the
%    distance in meters between locations (alon, alat)
%    and (blon, blat).  The default earth radius is
%    assumed to be 6371*1000 meters, the radius for
%    a sphere of equal-volume.
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

f{'nav_lon'}(:) = permute(lons(1:dimx-3,1:dimy-3,1),[2,1]);
f{'nav_lat'}(:) = permute(lats(1:dimx-3,1:dimy-3,1),[2,1]);

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ORIGINAL CODE from seagrid2roms.m taht computes the scaling %%%%%
%%%%%  required for ROMSfactors                                    %%%%%
% % Use "geometry" from seagrid file for computation of "pm", "pn".
% 
% % "geometry" contains distances computed from lon and lat in the
% % Seagrid routine "dosave.m" using the "earthdist.m" routine, which
% % assumes a spherical globe with a radius of 6371 km.
% 
% gx = geometry{1};   % Spherical distances in meters.
% gy = geometry{2};
% 
% 
% sx = 0.5*(gx(1:end-1, :) + gx(2:end, :));
% sy = 0.5*(gy(:, 1:end-1) + gy(:, 2:end));
% 
% pm = 1 ./ sx;
% pn = 1 ./ sy;
% 
% % pm and pn cannot be Inf, even if on land, so if values
% % are Inf, set to an arbitrary non-zero value
% pm(isinf(pm))=0.999e-3;
% pn(isinf(pn))=0.999e-3;
% pm(isnan(pm))=0.999e-3;
% pn(isnan(pn))=0.999e-3;
% nc{'pm'}(:) = pm;
% nc{'pn'}(:) = pn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%   FUNCTION earthdist.m   %%%%%%%%%%%%%%%%%%%%
% function theResult = earthdist(alon, alat, blon, blat, radius)
% 
% % earthdist -- Distance in meters between two lon/lats.
% %  earthdist(alon, aloat, blon, blat, radius) returns the
% %   distance in maters between locations (alon, alat)
% %   and (blon, blat).  The default earth radius is
% %   assumed to be 6371*1000 meters, the radius for
% %   a sphere of equal-volume.
% 
% if nargin < 4, help(mfilename), return, end
% if nargin < 5, radius = 6371*1000; end   % meters.
% 
% RCF = 180 / pi;
% 
% alon = alon / RCF;
% alat = alat / RCF;
% blon = blon / RCF;
% blat = blat / RCF;
% 
% c = cos(alat);
% ax = c .* cos(alon);
% ay = c .* sin(alon);
% az = sin(alat);
% 
% c = cos(blat);
% bx = c .* cos(blon);
% by = c .* sin(blon);
% bz = sin(blat);
% 
% result = acos(ax.*bx + ay.*by + az.*bz) .* radius;
% 
% if nargout > 0
% 	theResult = result;
% else
% 	disp(result)
% end
