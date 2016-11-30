% NAME: Compute_ScalingFactors.m
%
% AUTHOR: J.-P. Paquin 
%
% DATE: July 2013
%
% REVISIONS: October 2016 - Idalia A Machuca
%
% DESCRIPTION: Creates a coordinates.nc file for NEMO using a grid created 
%              by the great circle functions. Computes the latitudes and
%              longitudes of u, v, and f grids as simple averages and their
%              respective scaling factors (based on ROMS code).
%
% INPUTS: Grid file 
%         Original: Seagrid output file
%         Revised: Great Circle grid file
%
% OUTPUT: Coordinates.nc file in compatible format for NEMO
%
% CALLED PGM & SCRIPTS: earthdist.m : compute distances
%
% NOTES: The final grid size is 2 grid points less in each direction to
%        avoid NaNs to be present in the scaling factors. 
%
%--------------------------------------------------------------------------

clear all; clc

display('Compute_ScalingFactors.m is a revision of compute_grid_and_scaling_factors.m')

fileout = ('/ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/initial_sets/coordinates/coords_01.nc');
infile = ('/ocean/imachuca/Canyons/mackenzie_canyon/bathymetry/initial_sets/grid/grid_01.nc')

lon_T = ncread(infile, 'grid_lons');
lat_T = ncread(infile, 'grid_lats');

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
for vv=1:4
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

display('Compute distances between point for each grid T,u,v,f')
radius = 6371*1000;

dist_lat=NaN(dimx,dimy,4);
dist_lon=NaN(dimx,dimy,4);
for vv=1:4
for jj=1:dimy-1
for ii=1:dimx-1   
  alat=lats(ii  ,jj  ,vv) ;
  alon=lons(ii  ,jj  ,vv) ;
  
  blat=lats(ii+1,jj,vv) ;
  blon=lons(ii+1,jj,vv) ;
  dist_lat(ii,jj,vv) = earthdist(alon, alat, blon, blat, radius); % eq to gy

  blat=lats(ii  ,jj+1,vv) ; 
  blon=lons(ii  ,jj+1,vv) ;
  dist_lon(ii,jj,vv) = earthdist(alon, alat, blon, blat, radius); % eq to gx

end
end
end

display('Compute scaling factors (inspired by seagrid2roms.m)')
sx=NaN(dimx,dimy,4);
sy=NaN(dimx,dimy,4);
for vv=1:4 
  sx(1:dimx-1,1:dimy  ,vv) = 0.5*(dist_lon(1:end-1,   :   ,vv) + dist_lon(2:end,   :  ,vv));
  sy(1:dimx  ,1:dimy-1,vv) = 0.5*(dist_lat(   :   ,1:end-1,vv) + dist_lat(  :  , 2:end,vv));
end

sx(isinf(sx))=1.e+20;
sy(isinf(sy))=1.e+20;
sx(isnan(sx))=1.e+20;
sy(isnan(sy))=1.e+20;


display('Write netcdf file')
ncid = netcdf.create(fileout,'CLOBBER');

xx = length(1:dimx-3) ;
yy = length(1:dimy-3) ;

x = netcdf.defDim(ncid,'x',xx);
y = netcdf.defDim(ncid,'y',yy);
time = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));

nlon = netcdf.defVar(ncid,'nav_lon','float',[x,y]);
netcdf.putAtt(ncid,nlon,'units','degrees_east');
netcdf.putAtt(ncid,nlon,'comment','at t points');

nlat = netcdf.defVar(ncid,'nav_lat','float',[x,y]);
netcdf.putAtt(ncid,nlat,'units','degrees_north');
netcdf.putAtt(ncid,nlat,'comment','at t points');

ntim = netcdf.defVar(ncid,'time','float',time);
netcdf.putAtt(ncid,ntim,'units','seconds since 0001-01-01 00:00:00');
netcdf.putAtt(ncid,ntim,'time_origin','0000-JAN-01 00:00:00');
netcdf.putAtt(ncid,ntim,'calendar','gregorian');

nstp = netcdf.defVar(ncid,'time_steps','int',time);
netcdf.putAtt(ncid,nstp,'units','seconds since 0001-01-01 00:00:00');
netcdf.putAtt(ncid,nstp,'time_origin','0000-JAN-01 00:00:00');

glamt = netcdf.defVar(ncid,'glamt','double',[x,y,time]);
netcdf.putAtt(ncid,glamt,'missing_value',1.e+20);

glamu = netcdf.defVar(ncid,'glamu','double',[x,y,time]);
netcdf.putAtt(ncid,glamu,'missing_value',1.e+20);

glamv = netcdf.defVar(ncid,'glamv','double',[x,y,time]);
netcdf.putAtt(ncid,glamv,'missing_value',1.e+20);

glamf = netcdf.defVar(ncid,'glamf','double',[x,y,time]);
netcdf.putAtt(ncid,glamf,'missing_value',1.e+20);

gphit = netcdf.defVar(ncid,'gphit','double',[x,y,time]);
netcdf.putAtt(ncid,gphit,'missing_value',1.e+20);

gphiu = netcdf.defVar(ncid,'gphiu','double',[x,y,time]);
netcdf.putAtt(ncid,gphiu,'missing_value',1.e+20);

gphiv = netcdf.defVar(ncid,'gphiv','double',[x,y,time]);
netcdf.putAtt(ncid,gphiv,'missing_value',1.e+20);

gphif = netcdf.defVar(ncid,'gphif','double',[x,y,time]);
netcdf.putAtt(ncid,gphif,'missing_value',1.e+20);

e1t = netcdf.defVar(ncid,'e1t','double',[x,y,time]);
netcdf.putAtt(ncid,e1t,'missing_value',1.e+20);

e1u = netcdf.defVar(ncid,'e1u','double',[x,y,time]);
netcdf.putAtt(ncid,e1u,'missing_value',1.e+20);

e1v = netcdf.defVar(ncid,'e1v','double',[x,y,time]);
netcdf.putAtt(ncid,e1v,'missing_value',1.e+20);

e1f = netcdf.defVar(ncid,'e1f','double',[x,y,time]);
netcdf.putAtt(ncid,e1f,'missing_value',1.e+20);

e2t = netcdf.defVar(ncid,'e2t','double',[x,y,time]);
netcdf.putAtt(ncid,e2t,'missing_value',1.e+20);

e2u = netcdf.defVar(ncid,'e2u','double',[x,y,time]);
netcdf.putAtt(ncid,e2u,'missing_value',1.e+20);

e2v = netcdf.defVar(ncid,'e2v','double',[x,y,time]);
netcdf.putAtt(ncid,e2v,'missing_value',1.e+20);

e2f = netcdf.defVar(ncid,'e2f','double',[x,y,time]);
netcdf.putAtt(ncid,e2f,'missing_value',1.e+20);

netcdf.endDef(ncid)
netcdf.putVar(ncid,nlon,permute(lons(1:dimx-3,1:dimy-3,1),[2,1]))
netcdf.putVar(ncid,nlat,permute(lats(1:dimx-3,1:dimy-3,1),[2,1]))
netcdf.putVar(ncid,ntim,0,1,[0])
netcdf.putVar(ncid,nstp,0,1,[0]) 
netcdf.putVar(ncid,glamt,permute(lons(1:dimx-3,1:dimy-3,1)  ,[2,1]))
netcdf.putVar(ncid,glamu,permute(lons(1:dimx-3,1:dimy-3,2)  ,[2,1]))
netcdf.putVar(ncid,glamv,permute(lons(1:dimx-3,1:dimy-3,3)  ,[2,1]))
netcdf.putVar(ncid,glamf,permute(lons(1:dimx-3,1:dimy-3,4)  ,[2,1]))
netcdf.putVar(ncid,gphit,permute(lats(1:dimx-3,1:dimy-3,1)  ,[2,1]))
netcdf.putVar(ncid,gphiu,permute(lats(1:dimx-3,1:dimy-3,2)  ,[2,1]))
netcdf.putVar(ncid,gphiv,permute(lats(1:dimx-3,1:dimy-3,3)  ,[2,1]))
netcdf.putVar(ncid,gphif,permute(lats(1:dimx-3,1:dimy-3,4)  ,[2,1]))
netcdf.putVar(ncid,e1t,permute(dist_lon(1:dimx-3,1:dimy-3,1),[2,1]))
netcdf.putVar(ncid,e1u,permute(dist_lon(1:dimx-3,1:dimy-3,2),[2,1]))
netcdf.putVar(ncid,e1v,permute(dist_lon(1:dimx-3,1:dimy-3,3),[2,1]))
netcdf.putVar(ncid,e1f,permute(dist_lon(1:dimx-3,1:dimy-3,4),[2,1]))
netcdf.putVar(ncid,e2t,permute(dist_lat(1:dimx-3,1:dimy-3,1),[2,1]))
netcdf.putVar(ncid,e2u,permute(dist_lat(1:dimx-3,1:dimy-3,2),[2,1]))
netcdf.putVar(ncid,e2v,permute(dist_lat(1:dimx-3,1:dimy-3,3),[2,1]))
netcdf.putVar(ncid,e2f,permute(dist_lat(1:dimx-3,1:dimy-3,4),[2,1]))

netcdf.close(ncid)

display('END WRITING COORDINATES')
