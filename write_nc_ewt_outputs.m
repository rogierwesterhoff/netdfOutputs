% function write_nc_ewt_outputs(z,xaxis,yaxis,nc_filename)

%% made with the help of:
% http://code.guillaumemaze.org/tips/howtocreateacleannetcdffilefrommatlabusingthebuiltintoolbox
cc

performdeflate=false;

nc_filename='ewt_outputs_Mataura.nc';
if exist(nc_filename)
    warning('previous version existed. Deleting previous version...')
    delete(nc_filename)
end

% These parameters are loaded:
load('i:\GroundWater\Research\NIWA_NationalHydrologyProgram\Data\Task2_Data\2018-19_work\Matlab scripts\20190502_EWT_Mataura_250m_100yrs.mat');
plotmask = EWT.Plotmask;
transmissivity = EWT.transmissivity; transmissivity(~plotmask)=nan;
wtd = EWT.WTD; wtd(~plotmask)=nan;
head = EWT.Head; head(~plotmask)=nan;
K = EWT.K_m_s; K(~plotmask)=nan;
xaxis = EWT.xaxes;
yaxis = EWT.yaxes;
PotRech = EWT.Potential_recharge; PotRech(~plotmask)=nan;
ActRech = EWT.Actual_recharge; ActRech(~plotmask)=nan;
clear EWT

% Define axis, for example:
% nzseg_v3 = extractfield(S_DN3,'nzseg_v3');
% T = datenum(2007,1,1:365,0,0,0); % Time

% load('I:\GroundWater\Research\Smart\Satellite RS\Data_and_Research\New_Zealand_case\K_NZ_WGS84.mat')
% These parameters are loaded:
% lat, lon, K_WGS_30arcsec,sigmaK_WGS_30arcsec

% Define axis, for example:
% X = 0:2:360; % Longitude
% Y = -90:2:90; % Latitude
% T = datenum(2007:2007,7,1,0,0,0); % Time

% Create a random field to record in the cdf file dimensions are (T,Y,X):
% D=ones(1,length(yaxis),length(xaxis));
% D(1,:,:)=z;
% E=ones(1,length(lat),length(lon));
% E(1,:,:)=K_WGS_30arcsec;

ncid = netcdf.create(nc_filename,'netcdf4');

% Define useful constants:
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
fillValue = -9999;

% Define dimensions (this is the order they should be in!):
% dimidT = netcdf.defDim(ncid,'time',length(T));
% dimidZ = netcdf.defDim(ncid,'time',length(Z));
dimidY = netcdf.defDim(ncid,'y',length(yaxis));
dimidX = netcdf.defDim(ncid,'x',length(xaxis));

% Define axis variables:
% varid = netcdf.defVar(ncid,'time','double',[dimidT]);          
% netcdf.putAtt(ncid,varid,'standard_name','time');
% netcdf.putAtt(ncid,varid,'long_name','Time of measurements');
% netcdf.putAtt(ncid,varid,'units','days since 1900-01-01 00:00:00');
% netcdf.defVarFill(ncid,varid,false,fillValue);
% netcdf.putVar(ncid,varid,T-datenum(1900,1,1,0,0,0));

varid = netcdf.defVar(ncid,'y','float',[dimidY]);
netcdf.putAtt(ncid,varid,'standard_name','projection_y_coordinate');
netcdf.putAtt(ncid,varid,'long_name','Northing New Zealand Transverse Mercator');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.defVarFill(ncid,varid,false,fillValue);
netcdf.putVar(ncid,varid,yaxis);

varid = netcdf.defVar(ncid,'x','float',[dimidX]);
netcdf.putAtt(ncid,varid,'standard_name','projection_x_coordinate');
netcdf.putAtt(ncid,varid,'long_name','Easting New Zealand Transverse Mercator');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.defVarFill(ncid,varid,false,fillValue);
netcdf.putVar(ncid,varid,xaxis);

% check https://au.mathworks.com/help/matlab/ref/netcdf.defvar.html
% and https://publicwiki.deltares.nl/display/NETCDF/Coordinates
varid = netcdf.defVar(ncid,'transverse_mercator','char',[dimidY,dimidX]); 
netcdf.putAtt(ncid,varid,'grid_mapping_name','transverse_mercator');
netcdf.putAtt(ncid,varid,'projection_mapping_name','New Zealand Transverse Mercator (NZTM)');
netcdf.putAtt(ncid,varid,'semi_major_axis',6378137.0);
netcdf.putAtt(ncid,varid,'inverse_flattening',298.257222101);
netcdf.putAtt(ncid,varid,'latitude_of_projection_origin',0.0);
netcdf.putAtt(ncid,varid,'longitude_of_central_meridian',173.0);
netcdf.putAtt(ncid,varid,'scale_factor_at_central_meridian',0.9996);
netcdf.putAtt(ncid,varid,'false_easting',1.6e+06);
netcdf.putAtt(ncid,varid,'false_northing',1.0e+07);
netcdf.putAtt(ncid,varid,'spatial_ref','PROJCS[\"NZGD2000 / New Zealand Transverse Mercator 2000\",GEOGCS[\"NZGD2000\",DATUM[\"New_Zealand_Geodetic_Datum_2000\",SPHEROID[\"GRS 1980\",6378137,298.257222101,AUTHORITY[\"EPSG\",\"7019\"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY[\"EPSG\",\"6167\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4167\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",173],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",1600000],PARAMETER[\"false_northing\",10000000],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AUTHORITY[\"EPSG\",\"2193\"]]');
% netcdf.putAtt(ncid,varid,'proj4_params','+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs');
% netcdf.putAtt(ncid,varid,'EPSG_code','EPSG:2193');

%% define variables: 1
varname = 'hydraulic_head';
long_name = 'Hydraulic head in metres above mean sea level';
grid_mapping = 'transverse_mercator';
units = 'm';

% varid = netcdf.defVar(ncid,varname,'double',[dimidT,dimidY,dimidX]);  
varid = netcdf.defVar(ncid,varname,'float',[dimidY,dimidX]); 
netcdf.putAtt(ncid,varid,'long_name',long_name);
netcdf.putAtt(ncid,varid,'grid_mapping',grid_mapping);
netcdf.putAtt(ncid,varid,'units',units);
netcdf.defVarFill(ncid,varid,false,fillValue);
v = head;
    v(isnan(v)) = fillValue;
    netcdf.putVar(ncid,varid,v);  

%% define variables: 2
varname = 'water_table_depth';
long_name = 'water table depth in metres below ground level';
grid_mapping = 'transverse_mercator';
units = 'm';

% varid = netcdf.defVar(ncid,varname,'double',[dimidT,dimidY,dimidX]);  
varid = netcdf.defVar(ncid,varname,'float',[dimidY,dimidX]); 
netcdf.putAtt(ncid,varid,'long_name',long_name);
netcdf.putAtt(ncid,varid,'grid_mapping',grid_mapping);
netcdf.putAtt(ncid,varid,'units',units);
netcdf.defVarFill(ncid,varid,false,fillValue);
v = wtd;
    v(isnan(v)) = fillValue;
    netcdf.putVar(ncid,varid,v);  

%% define variables: 3
varname = 'transmissivity';
long_name = 'transmissivity';
grid_mapping = 'transverse_mercator';
units = 'm2 s-1';

% varid = netcdf.defVar(ncid,varname,'double',[dimidT,dimidY,dimidX]);  
varid = netcdf.defVar(ncid,varname,'float',[dimidY,dimidX]); 
netcdf.putAtt(ncid,varid,'long_name',long_name);
netcdf.putAtt(ncid,varid,'grid_mapping',grid_mapping);
netcdf.putAtt(ncid,varid,'units',units);
netcdf.defVarFill(ncid,varid,false,fillValue);
v = transmissivity;
    v(isnan(v)) = fillValue;
    netcdf.putVar(ncid,varid,v);  

%% define variables: 4
varname = 'hydraulic_conductivity';
long_name = 'estimate of saturated hydraulic conductivity (isotropic)';
grid_mapping = 'transverse_mercator';
units = 'm2 s-1';

% varid = netcdf.defVar(ncid,varname,'double',[dimidT,dimidY,dimidX]);  
varid = netcdf.defVar(ncid,varname,'float',[dimidY,dimidX]); 
netcdf.putAtt(ncid,varid,'long_name',long_name);
netcdf.putAtt(ncid,varid,'grid_mapping',grid_mapping);
netcdf.putAtt(ncid,varid,'units',units);
netcdf.defVarFill(ncid,varid,false,fillValue);
v = K;
    v(isnan(v)) = fillValue;
    netcdf.putVar(ncid,varid,v);

%% define variables: 5
varname = 'original_recharge';
long_name = 'groundwater recharge from the NGRM model (doi:10.3390/rs10010058)';
grid_mapping = 'transverse_mercator';
units = 'mm yr-1';

% varid = netcdf.defVar(ncid,varname,'double',[dimidT,dimidY,dimidX]);  
varid = netcdf.defVar(ncid,varname,'float',[dimidY,dimidX]); 
netcdf.putAtt(ncid,varid,'long_name',long_name);
netcdf.putAtt(ncid,varid,'grid_mapping',grid_mapping);
netcdf.putAtt(ncid,varid,'units',units);
netcdf.defVarFill(ncid,varid,false,fillValue);
v = PotRech;
    v(isnan(v)) = fillValue;
    netcdf.putVar(ncid,varid,v);

%% define variables: 6
varname = 'corrected_recharge';
long_name = 'groundwater recharge after correction from the groundwater model';
grid_mapping = 'transverse_mercator';
units = 'mm yr-1';

% varid = netcdf.defVar(ncid,varname,'double',[dimidT,dimidY,dimidX]);  
varid = netcdf.defVar(ncid,varname,'double',[dimidY,dimidX]); 
netcdf.putAtt(ncid,varid,'long_name',long_name);
netcdf.putAtt(ncid,varid,'grid_mapping',grid_mapping);
netcdf.putAtt(ncid,varid,'units',units);
netcdf.defVarFill(ncid,varid,false,fillValue);
v = ActRech;
    v(isnan(v)) = fillValue;
    netcdf.putVar(ncid,varid,v);

%% define global attributes
netcdf.putAtt(ncid,NC_GLOBAL,'title','EWT outputs of the Mataura Catchment')
netcdf.putAtt(ncid,NC_GLOBAL,'long_title','EWT outputs of the Mataura catchment.')
netcdf.putAtt(ncid,NC_GLOBAL,'comments','Made in the NZWaM-Hydro project, funded by MBIE, NZ. Citation to be used: Westerhoff, R., White, P., and Miguez-Macho, G.: Application of an improved global-scale groundwater model for water table estimation across New Zealand, Hydrol. Earth Syst. Sci., 22, 6449-6472, https://doi.org/10.5194/hess-22-6449-2018, 2018.')
netcdf.putAtt(ncid,NC_GLOBAL,'institution','GNS Science')
netcdf.putAtt(ncid,NC_GLOBAL,'source','Matlab script')

netcdf.putAtt(ncid,NC_GLOBAL,'Conventions','CF-1.6')   
netcdf.putAtt(ncid,NC_GLOBAL,'Conventions_help','http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html')

netcdf.putAtt(ncid,NC_GLOBAL,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'))
netcdf.putAtt(ncid,NC_GLOBAL,'CreatedBy',getenv('LOGNAME'))
netcdf.putAtt(ncid,NC_GLOBAL,'MatlabSource',which('write_nc_ewt_outputs.m'))

if performdeflate
    netcdf.defVarDeflate(ncid,varid,true,true,5);
end
netcdf.close(ncid)

ncdisp(nc_filename)