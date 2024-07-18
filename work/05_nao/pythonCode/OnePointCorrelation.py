#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: pamelatueam
"""

import cartopy.crs as ccrs
import netCDF4 as nc
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
import numpy as np
import cftime

# Define Earth's gravitational acceleration in m/s^2
g = 9.80665

# read the data
path_in = "work/05_nao/data/geopotential_data.nc"
ncfile = nc.Dataset(path_in)

geopotential = ncfile.variables['z'][:]  # shape (time, expver, latitude, longitude)
geopotential = geopotential/g
lons = ncfile.variables['longitude'][:]
lats = ncfile.variables['latitude'][:]

time = ncfile.variables['time']
time_data = time[:]

# Convert time data to datetime objects
time_data_cftime = cftime.num2date(time_data, units=time.units, calendar=time.calendar)

# Extract year-month format from cftime objects
year_month_list = [(time.year, time.month) for time in time_data_cftime]


djfm_indices = []

# Iterate over each tuple in year_month_list
for i, (year, month) in enumerate(year_month_list):
    # Check if the month is December, January, or February
    if month in [12, 1, 2, 3]:
        djfm_indices.append(i)


# Select the DJFM months
geopotential = geopotential[djfm_indices]

# Compute anomalies by removing the time-mean
geopotential_mean = geopotential.mean(axis=0)
geopotential_anomalie = geopotential - geopotential_mean

for var in geopotential_anomalie:
    print(var)

# Select a grid point (latitude and longitude indices)
lat_index = lat_idx = np.argmin(np.abs(lats - 65))  # Replace with the index of your chosen latitude
lon_index = np.argmin(np.abs(lons - (-18)))  # Replace with the index of your chosen longitude
lat_index
lon_index
# Extract the time series for the selected grid point
point_time_series = geopotential_anomalie[:, lat_index, lon_index]

    
# Initialize an array to store the correlation coefficients
correlation_map = np.empty((geopotential_anomalie.shape[1], geopotential_anomalie.shape[2]))

# Compute the correlation of the point_time_series with all other time series
for i in range(geopotential_anomalie.shape[1]):
    for j in range(geopotential_anomalie.shape[2]):
        grid_point_time_series = geopotential_anomalie[:,  i, j]
        correlation_map[i, j] = np.corrcoef(point_time_series, grid_point_time_series)[0, 1]
        
selected_point_corr = correlation_map[lat_index, lon_index]
grid_points = [(50, -200), (-60, -250)] 
# Plot the one-point correlation map
clevs = np.arange(-0.7, 1.0, 0.2)
proj = ccrs.Orthographic(central_longitude=-0, central_latitude=75)
ax = plt.axes(projection=proj)
ax.set_global()
ax.coastlines()
ax.gridlines()
ax.contourf(lons, lats, correlation_map, levels=clevs,
                           cmap=plt.cm.RdYlBu_r, transform=ccrs.PlateCarree(), extend='both')
# Plot positive correlation contours with solid lines
ax.contour(lons, lats, correlation_map, levels=clevs[clevs > 0],
           colors='k', linewidths=0.8, linestyles='solid', transform=ccrs.PlateCarree())
# Plot negative correlation contours with dashed lines
ax.contour(lons, lats, correlation_map, levels=clevs[clevs < 0],
           colors='k', linewidths=0.8, linestyles='dashed', transform=ccrs.PlateCarree())

# Add the correlation coefficients for the selected grid points
for lat_index, lon_index in grid_points:
    corr_coeff = correlation_map[lat_index, lon_index]
    ax.text(lons[lon_index], lats[lat_index], f'{corr_coeff:.2f}',
            transform=ccrs.PlateCarree(), horizontalalignment='center',
            verticalalignment='center', fontsize=12, color='black')

plt.title('One-Point Correlation Map (DJFM)', fontsize=16)
plt.show()
