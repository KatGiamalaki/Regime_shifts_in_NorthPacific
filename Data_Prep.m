%%  Written for daily SLP  by NCEP/NCAR reanalysis I (https://psl.noaa.gov/cgi-bin/db_search/DBListFiles.pl?did=198&tid=94405&vid=676) 
% % From 01/01/1948 until 19/11/2016 - change accordingly 

lon = ncread('slp.1948.nc','lon');
lat = ncread('slp.1948.nc','lat');
total_ndays = 25143;

slp = NaN * ones(size(lon,1),size(lat,1),total_ndays);

%% Upload data into one array dims: lon * lat * total_ndays

a = 1;
for ii = 1948:2016
    temp_slp = ncread(['slp.' num2str(ii) '.nc'],'slp');
    if size(temp_slp,3)==365
        slp(:,:,a:a+364) = temp_slp;
        a  = a+365;
    elseif size(temp_slp,3)==366
        slp(:,:,a:a+365) = temp_slp;
        a = a+366;
    end
end

%% North Pacific area selection

ji = find(lon>100 & lon<245); 
ij = find(lat<65 & lat>15);
NP_slp=slp(ji,ij,:);

%% 

clearvars -except NP_slp

