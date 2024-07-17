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
import pandas as pd
import cftime
from eofs.standard import Eof
import cartopy.feature as cfeature

# load the data
path_in = "work/05_nao/data/mean_sea_level_data.nc"
ncfile = nc.Dataset(path_in)

# extract the mean sea level pressure, latitude, logitude and time
mean_sea_level = ncfile.variables['msl'][:]
lons = ncfile.variables['longitude'][:]
lats = ncfile.variables['latitude'][:]
time = ncfile.variables['time']
time_data = time[:]

# change the date format
time_data_cftime = cftime.num2date(time_data, units=time.units, calendar=time.calendar)

# Extract year-month format from cftime objects
year_month_list = [(time.year, time.month) for time in time_data_cftime]

# Initialize an empty list to store the indices of DJFM months
djfm_indices = []

# Iterate over each tuple in year_month_list
for i, (year, month) in enumerate(year_month_list):
    # Check if the month is December, January, or February
    if month in [12, 1, 2, 3]:
        djfm_indices.append(i)

# Convert the list of indices to a numpy array
djfm_indices_array = np.array(djfm_indices)

# Select the DJFM months
mean_sea_level_djfm = mean_sea_level[djfm_indices]


# Compute anomalies by removing the time-mean.
mean_sea_level_mean = mean_sea_level_djfm.mean(axis = 0)
mean_sea_level_anomalie = mean_sea_level_djfm - mean_sea_level_mean

# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(lats)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(mean_sea_level_anomalie, weights=wgts)

# Retrieve the leading EOF, expressed as the covariance between the leading PC
# time series and the input MSLP anomalies at each grid point.
eof1 = solver.eofsAsCovariance(neofs=1)
eof1_2d = eof1[:, 0, :]

# Get the explained variance associated with eof1
explained_variance = solver.varianceFraction()
eof1_explained_variance = explained_variance[0]*100

# Plot the leading EOF expressed as covariance in the European/Atlantic domain.

clevs = np.linspace(-900, 1000, 20)  #generates contour levels for the plot ranging from -70 to 80 with 11 intervals


proj = ccrs.Orthographic(central_longitude=-0, central_latitude=75)
ax = plt.axes(projection=proj)
ax.coastlines()
ax.set_global()
ax.gridlines()
fill=ax.contourf(lons, lats, eof1.squeeze(), levels=clevs,
cmap=plt.cm.RdYlBu_r, transform=ccrs.PlateCarree())

ax.contour(lons, lats, eof1.squeeze(), levels=clevs[clevs > 0],
           colors='grey', linewidths=1, linestyles='solid', transform=ccrs.PlateCarree())

# Plot negative correlation contours with dashed lines (dashed line for blue)
ax.contour(lons, lats, eof1.squeeze(), levels=clevs[clevs < 0],
           colors='gray', linewidths=1, linestyles='dashed', transform=ccrs.PlateCarree())

plt.text(0.5, 1.05, f'EOF1 MSLP DJFM 1940-2024: {eof1_explained_variance:.1f}', transform=ax.transAxes,
       horizontalalignment='center', fontsize=20)
plt.show()
