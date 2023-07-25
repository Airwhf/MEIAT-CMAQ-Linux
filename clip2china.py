import xarray as xr


dataset = xr.open_dataset('Input/SA/landscan-global-2017_nodata.nc')
dataset = dataset.sel(lon=slice(70, 150), lat=slice(70, 0))

# print(dataset)
dataset.to_netcdf('Input/SA/china_population_2017.nc')
