import xarray as xr
import pandas as pd
import numpy as np
import geopandas as gpd
import os
import matplotlib.pyplot as plt

path = 'data/raw/meteo'
# ds = xr.open_dataset(os.path.join('..',path,'2023_ERA5_Regen.nc'))
ds = xr.open_dataset('/Users/lemarx/Documents/01_projects/Seminar_ClimateNStatistics/work/07-hydroLSTM/data/raw/meteo/2023_ERA5_Regen.nc')

catch_path = os.path.join('..','data/raw/catchments')
# polygons = gpd.read_file(os.path.join(catch_path,'regen_catchments_right.shp'))
polygons = gpd.read_file('/Users/lemarx/Documents/01_projects/Seminar_ClimateNStatistics/work/07-hydroLSTM/data/raw/catchments/regen_catchments_right.shp')

def calc_zonal_stats(geom,ds,time):
    #some weird geopandas things
    gdf = gpd.GeoDataFrame(geometry=[geom])
    geom = gdf['geometry']
    area = geom.area
    lon_list = np.arange(12,14.25,0.25)
    lat_list = np.arange(48,50.25,0.25)
    splits = np.empty(0)
    totals = {}
    for var in list(ds.keys()):
        for lat in lat_list:
            for lon in lon_list:
                clip_area = geom.clip_by_rect(xmin = lon, xmax = lon + 0.25, ymin = lat, ymax = lat + 0.25).area.values[0]
                if clip_area == 0:
                    continue
                var_value = ds[var].sel(latitude=lat, longitude=lon,time = time, method='nearest',tolerance = 0.001).values
                np.append(splits,(clip_area/area) * var_value)
        totals[var] = splits.sum()
    return totals


ds_daily = xr.Dataset({
    'tp': ds['tp'].resample(time='D').mean(dim='time'),
    'ssr': ds['ssr'].resample(time='D').mean(dim='time'),
    'd2m': ds['d2m'].resample(time='D').mean(dim='time'),
    't2m': ds['t2m'].resample(time='D').mean(dim='time'),
    'sp': ds['sp'].resample(time='D').mean(dim='time')
})

polygons.apply(lambda row: calc_zonal_stats(geom = row['geometry'],ds = ds_daily,time = '2023-01-01'),axis = 1)
