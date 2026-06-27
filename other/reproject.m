addpath('C:\Users\aashr\ISSM-Windows-MATLAB\bin');

load('data_pinns_Amery.mat'); 

epsgCode = 32742;  
proj = projcrs(epsgCode);

sz = size(x_gd);
x_vec = x_gd(:);
y_vec = y_gd(:);

[lat_vec, lon_vec] = projinv(proj, x_vec, y_vec);

lat = reshape(lat_vec, sz);
lon = reshape(lon_vec, sz);

figure('Name','Scatter Lon/Lat');
scatter(lon_vec, lat_vec, 6, 'filled');
axis equal;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title(sprintf('Reprojected Points (EPSG:%d) - Scatter', epsgCode));