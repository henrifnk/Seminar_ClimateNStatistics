#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: pamelatueam
"""
# Import the tools we are going to need 

import cartopy.crs as ccrs
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cftime
from sklearn.cluster import KMeans
from eofs.standard import Eof
import pandas as pd

# Load the data
path_in = "work/05_nao/data/mean_sea_level_data.nc"
ncfile = nc.Dataset(path_in)

mean_sea_level = ncfile.variables['msl'][:]  # shape (time, expver, latitude, longitude)
lons = ncfile.variables['longitude'][:]
lats = ncfile.variables['latitude'][:]

time = ncfile.variables['time']
time_data = time[:]


# Convert time data to datetime objects
time_data_cftime = cftime.num2date(time_data, units=time.units, calendar=time.calendar)

# Select the DJFM (December-January-February-March) months
djf_indices = [i for i, t in enumerate(time_data_cftime) if t.month in [12, 1, 2,3]]

mean_sea_level = mean_sea_level[djf_indices]

# Compute anomalies by removing the time-mean
mean_sea_level_mean = mean_sea_level.mean(axis=0)
mean_sea_level_anomalie = mean_sea_level - mean_sea_level_mean

# Flatten the spatial dimensions
n_samples,  n_lat, n_lon = mean_sea_level_anomalie.shape

# Compute the EOFs
coslat = np.cos(np.deg2rad(lats)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(mean_sea_level_anomalie, weights=wgts)
eofs = solver.eofs(neofs=4)  # Shape (4, n_lat, n_lon)
pcs = solver.pcs(npcs=4)     # Shape (n_samples, 4)

# Perform k-means clustering on the PCs
n_clusters = 4
kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(pcs)

# Get the cluster labels
labels = kmeans.labels_

# Calculate the mean PCs for each cluster
mean_pcs = np.array([pcs[labels == i].mean(axis=0) for i in range(n_clusters)])

# Project the mean PCs back to the spatial domain
nao_regimes = np.dot(mean_pcs, eofs.reshape(4, -1))

# Reshape nao_regimes to the original spatial dimensions
nao_regimes = nao_regimes.reshape(n_clusters, n_lat, n_lon)

explained_variance = solver.varianceFraction()

for i, var in enumerate(explained_variance[:4]):
    print(f"Explained variance of EOF{i+1}: {var:.2%}")

# Plot the 4 NAO regimes
proj = ccrs.Orthographic(central_longitude=-20, central_latitude=60)
fig, axes = plt.subplots(2, 2, subplot_kw={'projection': proj}, figsize=(15, 10))

regime_indices = [0,1,2, 3]
clevs = np.linspace(-800, 900, 20)
for i, ax in zip(regime_indices, axes.flat):
    ax.set_global()
    ax.coastlines()
    ax.gridlines()
    regime_plot = ax.contourf(lons, lats, nao_regimes[i], levels=clevs, cmap=plt.cm.RdYlBu_r, transform=ccrs.PlateCarree())
    ax.contour(lons, lats, nao_regimes[i], levels=clevs[clevs > 0],
               colors='grey', linewidths=1, linestyles='solid', transform=ccrs.PlateCarree())
    ax.contour(lons, lats, nao_regimes[i], levels=clevs[clevs < 0],
               colors='gray', linewidths=1, linestyles='dashed', transform=ccrs.PlateCarree())
    ax.set_title(f'{"NAO-" if i == 0 else "NAO+" if i==1 else "Blocking" if i == 2 else "Atl.Ridge"}\nExplained Variance: {explained_variance[i]:.2%}', fontsize=16)

fig.suptitle('4 MSLP Weather Regimes DJFM 1940-2024', fontsize=20)
#fig.colorbar(regime_plot, ax=axes, orientation='horizontal', fraction=0.046, pad=0.04)

plt.show()
