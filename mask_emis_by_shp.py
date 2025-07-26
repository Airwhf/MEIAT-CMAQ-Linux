import os
import glob
import re
import geopandas as gpd
import regionmask
import xarray as xr
import numpy as np

if __name__ == '__main__':
    # Define the output directory.
    output_dir = '/Volumes/project_new/data_pkucn/MIX_MEIAT_masked'
    os.makedirs(output_dir, exist_ok=True)
    
    # Method
    method = 'remove' # [keep or remove]
                    # Keep the value in the shapefile extent.
                    # Remove the value in the shapefile extent.
    
    # Define the shpefile path
    shp_for_reg_emis = 'shapefiles/CN-Province_WGS84.shp'
    
    # Read the shpefile
    gdf = gpd.read_file(shp_for_reg_emis)
    
    # Regional emission file paths.
    reg_emis_files = glob.glob('/Volumes/project_new/data_pkucn/mix_emis_2017/mix_emis_2017_meiat/*.nc')
    
    for reg_emis_file in reg_emis_files:
        reg_emis_ds = xr.open_dataset(reg_emis_file)
        mask = regionmask.mask_geopandas(gdf, reg_emis_ds.coords['longitude'], reg_emis_ds.coords['latitude'])
        if method == 'keep':
            mask = np.isnan(mask.values) # outside of the shapefile is True
        elif method == 'remove':
            mask = ~np.isnan(mask.values) # inside of the shapefile is True
        else:
            raise ValueError(f"Unknown method: {method}. Please use 'keep' or 'remove'.")
        
        # Create new emission files.
        output_name = f'{output_dir}/{os.path.basename(reg_emis_file)}'
        
        # If mask == True, set to 0; if mask == False, keep original value
        masked_emission = xr.where(mask, 0, reg_emis_ds['emission'])
        
        masked_emission.to_netcdf(output_name)
        print(f"Masked emission file saved to: {output_name} (method: {method})")
    
    
    