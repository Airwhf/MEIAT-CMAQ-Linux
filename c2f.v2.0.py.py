import os
import re
import pyioapi
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon 
import pyproj
import xarray as xr
import numpy as np
import PseudoNetCDF as pnc
import datetime
import glob
import f90nml

os.environ['IOAPI_ISPH'] = '6370000.'

def user_control():
    import datetime
    specified_time = datetime.datetime(2025, 6, 30, 23, 59)
    # 获取当前时间
    current_time = datetime.datetime.now()
    # 检查当前时间是否已经超过指定时间
    if current_time > specified_time:
        return False
    else:
        return True

def zonal(file, coarse_shpf, griddesc, gridname, fid):
    """_summary_

    Args:
        file (_type_): This is the path of NetCDF file. Please note that the netCDF file must with lon, lat and z variables. 
        coarse_shpf (_type_): This is the path of shapefile.
        griddesc (_type_): This is the file output from MCIP.
        gridname (_type_): This is the grid name of model.
        fid (_type_): This is the label of shapefile.

    Returns:
        _type_: pandas.DataFrame.
    """    
    dataset = xr.open_dataset(file)
    
    print("Original dataset:")
    print(dataset)
    print()
    
    print('Input dataset range:')
    lon_min, lon_max = dataset.lon.values.min(), dataset.lon.values.max()
    lat_min, lat_max = dataset.lat.values.min(), dataset.lat.values.max()
    print(f"lon_min, lon_max = {lon_min}, {lon_max}")
    print(f"lat_min, lat_max = {lat_min}, {lat_max}")
    print()
    
    # TODO Finish You can update for limit the range of dataset by longitude and latitude.
    # Create a Template
    gf = pnc.pncopen(griddesc, format='griddesc', GDNAM=gridname, SDATE=1970001)
    lons, lats = gf.variables['longitude'], gf.variables['latitude']
    lon_min, lon_max = lons.min(), lons.max()
    lat_min, lat_max = lats.min(), lats.max()
    print("Model range:")
    print(f"lon_min, lon_max = {lon_min}, {lon_max}")
    print(f"lat_min, lat_max = {lat_min}, {lat_max}")
    print()
    dataset = dataset.sel(lon=slice(lon_min, lon_max), lat=slice(lat_max, lat_min))
    print("Limited dataset:")
    print(dataset)
    print()
    
    # TODO interpolate the data and convert units to t/yr.
    # Calculate the new latitude and longitude coordinates with 0.01 degree interval
    lat_new = np.arange(lat_min, lat_max + 0.01, 0.01)
    lon_new = np.arange(lon_min, lon_max + 0.01, 0.01)
    dataset['z'] = dataset['z'].interp(lat=lat_new, lon=lon_new, method='linear')
    # Convert kg/km2/s to t/yr.
    dataset['z'] = dataset['z'] * 0.001 * 31536000 * 1.1131944444444444 ** 2
    
    
    lons, lats = np.meshgrid(dataset.lon.values, dataset.lat.values)
    # z = dataset.z.values
    # Create a DataFrame from the 'z' variable in the NetCDF dataset
    df_z = dataset['z'].to_dataframe().reset_index()
    
    gdf = gpd.read_file(coarse_shpf)
    gdf = gdf.to_crs(epsg=4326)
    # Flatten the grid to 1D and pair up the latitudes and longitudes
    latlon = np.dstack([lons.flatten(), lats.flatten()])[0]

    # Create a GeoDataFrame from the latitude and longitude pairs
    gdf_points = gpd.GeoDataFrame(pd.DataFrame(latlon, columns=['lon', 'lat']), 
                                geometry=gpd.points_from_xy(latlon[:,0], latlon[:,1]))
    # Set the same CRS for the points as the polygons
    gdf_points.crs = gdf.crs
    # Perform the spatial join
    gdf_sjoin = gpd.sjoin(gdf_points, gdf, how='left', op='within')
    # Merge the 'z' values into the spatial join result
    gdf_sjoin_z = gdf_sjoin.merge(df_z, how='left', on=['lon', 'lat'])

    # Group by the polygon index and sum the 'z' values
    df_sum = gdf_sjoin_z.groupby('index_right')['z'].sum().reset_index()
    df_sum[fid] = gdf.loc[df_sum.index_right.values][fid].values
    return df_sum


def draw_model_grid(griddesc_file):
    """_summary_

    This function will create a shapefile for model grid.
    
    Args:
        griddesc_file (_type_): This is the file path of GRIDDESC file from MCIP.
    """
    
    # Read the GRIDDESC file.
    griddesc = pyioapi.GRIDDESC(griddesc_file)
    # Get grid parameters
    grid = griddesc.get_grid(list(griddesc.grids)[0])
    # Get coordinate system parameters
    coord = griddesc.get_coord(list(griddesc.coords)[0])
    # Set output names based on grid and coordinate system names
    outshp_name = list(griddesc.grids)[0]

    # Initialise the starting grid quadrant range
    ringXleftOrigin = grid.XORIG
    ringXrightOrigin = ringXleftOrigin + grid.XCELL
    ringYbottomOrigin = grid.YORIG
    ringYtopOrigin = ringYbottomOrigin + grid.YCELL
    
    # Traverse the columns and write each column to the grid
    col = 1
    idd_num = 1
    grids_list = []  # Used to store the individual grids generated
    while col<=grid.NCOLS:
        # Initialise, initialise the upper and lower ranges for each column write completion
        ringYtop = ringYtopOrigin
        ringYbottom = ringYbottomOrigin
        # Iterate through the rows, creating and writing for each grid row in this column
        row = 1
        while row<=grid.NROWS:
            # Create the first grid in the upper left corner
            ring = Polygon([(ringXleftOrigin, ringYbottom),
                             (ringXleftOrigin, ringYtop ),
                             (ringXrightOrigin, ringYtop),
                             (ringXrightOrigin, ringYbottom)])
            # Write geometric polygons
            geo_df = gpd.GeoDataFrame(data=[[idd_num,row,col]],
                          geometry=[ring],
                          columns=['ID','rownum', 'colnum'])
            grids_list.append(geo_df)
    
            #Next polygon, update top and bottom ranges
            row+=1
            ringYtop = ringYtop + grid.YCELL
            ringYbottom = ringYbottom + grid.YCELL
            
            idd_num+=1
        # After a column write is complete, next column, update left and right ranges
        col+=1
        ringXleftOrigin = ringXleftOrigin + grid.XCELL
        ringXrightOrigin = ringXrightOrigin + grid.XCELL
    # Merge grids in the list
    # gdf_fishgrid = pd.concat(grids_list)
    gdf_fishgrid = gpd.GeoDataFrame(pd.concat(grids_list), geometry='geometry')


    # Setting up a mesh projection using GRIDDESC's parameters
    wrf_proj = pyproj.Proj('+proj=lcc '+'+lon_0='+str(coord.XCENT)+' lat_0='+str(coord.YCENT)
                    +' +lat_1='+str(coord.P_ALP)+' +lat_2='+str(coord.P_BET))
    # Setting up the grid projection
    gdf_fishgrid.crs = wrf_proj.crs
    # Output model mesh shp file
    output_name = 'Output/' + outshp_name + '_grid.shp'
    gdf_fishgrid.to_file(output_name,
                         driver='ESRI Shapefile',
                         encoding='utf-8')
    return output_name, grid.XCELL, grid.YCELL


def define_parent_grid(model_grid, parent_grid, gdnm):
    # Load the shapefiles into GeoDataFrames
    grid_gdf = gpd.read_file(model_grid)
    MEIC_gdf = gpd.read_file(parent_grid)
    output_name = f'Output/{gdnm}_ParentGrid.csv'

    # Check the data
    # grid_gdf.head(), MEIC_gdf.head()

    # Make sure the GeoDataFrames are using the same Coordinate Reference System (CRS)
    MEIC_gdf = MEIC_gdf.to_crs(grid_gdf.crs)

    # Create a new column in the grid GeoDataFrame to hold the name of the MEIC region each grid cell belongs to
    grid_gdf["PR"] = None  # * PR == Parent Region

    # For each grid cell, find the MEIC region it belongs to
    for i, row in grid_gdf.iterrows():
        # Get the geometry of the grid cell
        grid_geometry = row.geometry

        # Find the MEIC regions that intersect with the grid cell
        intersecting_regions = MEIC_gdf[MEIC_gdf.geometry.intersects(grid_geometry)]

        # If there is one intersecting region, set the name of the MEIC region in the grid GeoDataFrame
        if len(intersecting_regions) == 1:
            grid_gdf.at[i, "PR"] = intersecting_regions.iloc[0].NAME
        # If there are multiple intersecting regions, find the one with the largest intersection area
        elif len(intersecting_regions) > 1:
            intersection_areas = intersecting_regions.geometry.intersection(grid_geometry).area
            largest_intersection_index = intersection_areas.idxmax()
            grid_gdf.at[i, "PR"] = intersecting_regions.loc[largest_intersection_index].NAME

    grid_gdf[['ID','rownum','colnum','PR']].to_csv(output_name)
    return grid_gdf[['ID','rownum','colnum','PR']]


def allocator_factor(allocator_file, sector, grid_info, allocator_type):
    
    if allocator_type == 'raster':
        #* Raster allocator.
        # Statistics coarse values of allocator.
        parent_df = zonal(allocator_file, zone_file, gridf, gdnm, 'NAME')
        new_column_names = {'z': 'TotalValue'}
        parent_df.rename(columns=new_column_names, inplace=True)
        
        # Statistics fine values of allocator.
        son_df = zonal(allocator_file, model_grid, gridf, gdnm, 'ID')
        new_column_names = {'z': 'Finevalue'}
        son_df.rename(columns=new_column_names, inplace=True)
        
        # Add Total value to fine value table.
        grid_info = grid_info.set_index(grid_info.ID.values)
        parent_grids = []
        for i in son_df.ID.values:
            parent_grid = grid_info.loc[i, 'PR']
            parent_grids.append(str(parent_grid))
        son_df['NAME'] = parent_grids
        factor_df = pd.merge(son_df, parent_df, on='NAME', how='left')
        factor_df['AF'] = factor_df['Finevalue'] / np.where(factor_df['TotalValue']==0, 1, factor_df['TotalValue'])
        factor_df = pd.merge(factor_df, grid_info, on='ID', how='left')
        factor_df.to_csv(f'{factor_dir}/{sector}.csv')
    # elif allocator_type == 'line':     
    #      #* Line allocator.  
    return factor_df


def encode_title(file_name):
    # Get the species name from file name.
    condition = f"(.*?)_(.*?)_(.*?)__(.*?)__(.*?).csv"
    encode_name = re.findall(condition, file_name)[0]
    label = encode_name[0]
    year = encode_name[1]
    month = encode_name[2]
    sector = encode_name[3]
    pollutant = encode_name[4]
    dict = {"label": label,
            "year": year,
            "month": month,
            "sector": sector,
            "pollutant": pollutant}
    return dict


def create_source(label, year, mm, development_list, emission_factor_dir, emission_data_dir, out_path):
    for development in development_list:
        emission_factor_file = f'{emission_factor_dir}/{development}.csv'
        # if os.path.exists(f"{out_path}/source-{label}-{development}-{year}-{mm}.csv"):
        #     continue

        ef_data = pd.read_csv(emission_factor_file)
        ef_data['NAME'] = np.array(ef_data['NAME'], dtype='str')
        
        result = pd.DataFrame(columns=["LON", "LAT"])
        result['LON'] = ef_data['LON']
        result['LAT'] = ef_data['LAT']

        files = glob.glob(f'{emission_data_dir}/{label}_{year}_{mm}__{development}__*.csv')

        for file in files:
            specie = encode_title(os.path.basename(file))["pollutant"]
            data = pd.read_csv(file)   
            data['NAME'] = np.array(data['NAME'], dtype='str')
            
            #* This is the new version.
            _result = pd.merge(ef_data, data, on='NAME', how='left')
            _result['FineEmission'] = _result['AF'] * _result['z']
            _result['FineEmission'].values[np.isnan(_result['FineEmission'])] = 0
            # print(ef_data, data)
            result[specie] = _result['FineEmission']
            
            #* This is the old version.
            # values = []
            # for i in range(ef_data.shape[0]):
            #     temp_area = ef_data['NAME'].values[i]
            #     temp_ef = ef_data['AF'].values[i]
            #     temp_big_emis = data[data['NAME'].isin([temp_area])]['z'].values
            #     temp_small_emis = temp_ef * temp_big_emis
            #     try:
            #         values.append(temp_small_emis[0])
            #     except IndexError:
            #         values.append(0.0)
            # result[specie] = values

        result.to_csv(f"{out_path}/source-{label}-{development}-{year}-{mm}.csv", index=False)


def create_emission_file(emission_date, grid_desc, grid_name, label, sector, inventory_year, inventory_mechanism, target_mechanism):
    # Convert the emission date to other format.
    datetime_emission = pd.to_datetime(emission_date)
    yyyymmdd = datetime.datetime.strftime(datetime_emission, "%Y%m%d")
    # yyyy = datetime.datetime.strftime(datetime_emission, "%Y")
    mm = datetime.datetime.strftime(datetime_emission, "%m")
    # dd = datetime.datetime.strftime(datetime_emission, "%d")
    yyyyjjj = datetime.datetime.strftime(datetime_emission, "%Y%j")
    w = datetime.datetime.strftime(datetime_emission, "%w")

    # Create template file.
    gf = pnc.pncopen(grid_desc, GDNAM=grid_name, format="griddesc", SDATE=int(yyyyjjj), TSTEP=10000, withcf=False)
    gf.updatetflag(overwrite=True)
    tmpf = gf.sliceDimensions(TSTEP=[0] * 25)
    max_col_index = getattr(tmpf, "NCOLS") - 1
    max_row_index = getattr(tmpf, "NROWS") - 1

    # Create the source file and read it.
    source_file = f"output/source/source-{label}-{sector}-{inventory_year}-{mm}.csv"
    data = pd.read_csv(source_file)

    # Add I and J coordinate and calculate the total emission.
    I, J = gf.ll2ij(data["LON"].values, data["LAT"].values)
    # data["I"], data["J"] = data["rownum"], data["colnum"]
    data["I"], data["J"] = I, J
    celltotal = data.groupby(["I", "J"], as_index=False).sum()
    celltotal = celltotal[
        (celltotal["I"] >= 0)
        & (celltotal["J"] >= 0)
        & (celltotal["I"] <= max_col_index)
        & (celltotal["J"] <= max_row_index)
        ]
    # Read species file.
    species_file = (
        f"species/{inventory_mechanism}_{target_mechanism}_speciate_{sector}.csv"
    )
    species_info = pd.read_csv(species_file)
    fname_list = species_info.pollutant.values
    var_list = species_info.emission_species.values
    factor_list = species_info.split_factor.values
    divisor_list = species_info.divisor.values
    origin_units = species_info.inv_unit.values
    target_units = species_info.emi_unit.values

    # Read the temporal file.
    # _monthly_factor = pd.read_csv("temporal/monthly.csv")
    _weekly_factor = pd.read_csv("temporal/weekly.csv")
    _hourly_factor = pd.read_csv("temporal/hourly.csv")
    # monthly_factor = _monthly_factor[sector].values
    weekly_factor = _weekly_factor[sector].values
    hourly_factor = _hourly_factor[sector].values

    # Loop the species and create the variable to IOAPI file.
    items = zip(
        fname_list, var_list, factor_list, divisor_list, origin_units, target_units
    )
    for fname, var, split_factor, divisor, origin_unit, target_unit in items:
        try:
            # Extract the current pollutant.
            df = celltotal[["I", "J", fname]]

            # Convert monthly emission to weekly emission.
            weekly_values = df[fname].values * 0.25  # Version 1.0 and Version 1.1
            # weekly_values = converted_values * 0.25

            # Convert weekly emission to daily emission.
            daily_values = weekly_values * weekly_factor[int(w)]

            # Convert daily emission to hourly emission.
            df_list = []
            for hour_i in range(24):
                _df = pd.DataFrame(columns=["J", "I", "hour", "values"])
                _df["J"] = df.J.values
                _df["I"] = df.I.values
                _df["hour"] = np.zeros(df.shape[0]) + hour_i
                _df["values"] = daily_values * hourly_factor[hour_i]
                df_list.append(_df)
            result = pd.concat(df_list)
            # print(result)

            # Convert original units to target units and input the split_factor.
            if origin_unit == "Mmol" and target_unit == "mol/s":
                result["values"] = result["values"] * 1000000.0 / 3600.0 * split_factor
            elif origin_unit == "Mg" and target_unit == "g/s":
                result["values"] = result["values"] * 1000000.0 / 3600.0 * split_factor
            elif origin_unit == "Mg" and target_unit == "mol/s":
                result["values"] = (
                        result["values"] * 1000000.0 / 3600.0 / divisor * split_factor
                )

            # Convert the I, J, hour to int.
            result[["hour", "J", "I"]] = result[["hour", "J", "I"]].astype("int")
            h = result.hour
            i = result.I
            j = result.J

            # Create the variable of emission.
            evar = tmpf.createVariable(var, "f", ("TSTEP", "LAY", "ROW", "COL"))
            if target_unit == "mol/s":
                evar.setncatts(dict(units="moles/s", long_name=var, var_desc=var))
            elif target_unit == "g/s":
                evar.setncatts(dict(units="g/s", long_name=var, var_desc=var))

            evar[h, h * 0, j, i] = result["values"].values

        except KeyError:
            # If don not have this pollutant in GeoTIFF, skip it.
            print(f"Warning: Do not have the pollutant named {fname}.")
            continue
    # Get rid of initial DUMMY variable
    del tmpf.variables["DUMMY"]

    # Update TFLAG to be consistent with variables
    tmpf.updatetflag(tstep=10000, overwrite=True)

    # Remove VAR-LIST so that it can be inferred
    delattr(tmpf, "VAR-LIST")
    tmpf.updatemeta()

    # Save out.
    os.makedirs('Output/ModelReady', exist_ok=True)
    output_name = f"Output/ModelReady/{target_mechanism}_{sector}_{grid_name}_{yyyymmdd}.nc"
    tmpf.save(output_name, format="NETCDF3_CLASSIC")
            
            
if __name__ == "__main__":
    # --------------------------------------------------------------------------------------------------------
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!Author Information!!!!!!!!!!!!!!!!!!!!!!!!!!')
    if user_control() is True:
        print("### This system is developed by Haofan Wang.           ###")
        print("### Email: wanghf58@mail2.sysu.edu.cn                  ###")
        print("### -----------------Acknowledgements----------------- ###")
        print("1. Yang Zhang from College of Resources and Environment, Chengdu University of Information Technology")
        print("2. Yuhong Qiao from Chongqing Insitute of Eco-Environmental Science")
        print("3. Kai Wu from Department of Civil and Environmental Engineering, University of California, Irvine.")
        print("### -------------------------------------------------- ###")
    else:
        print("### This system is developed by Haofan Wang.           ###")
        print("### You can contact me for any suggestions.            ###")
        print("### Email: wanghf58@mail2.sysu.edu.cn                  ###")
        print("### ************************************************** ###")
        print("### The current version has expired.                   ###")
        print("### Please contact me to request the latest version.   ###")
        print("### ************************************************** ###")
        exit()
    # --------------------------------------------------------------------------------------------------------
    print()
    #TODO Read the namelist.input
    input_file = f90nml.read("namelist.input")
    
    #* This is the global input parameters.
    gridf = input_file["global"]['griddesc_file']
    gdnm = input_file["global"]['griddesc_name']
    start_date = input_file["global"]['start_date']
    end_date = input_file["global"]['end_date']
    
    #* In the 'model grid' module, it only use the `gridf` and `gdnm`.
    
    #* In the `Define the parent grid` module, it use the below parameter and the model grid output.
    zone_file = input_file["global"]['coarse_grid_file']
    
    #* In the `Coarse emission` module, it use the below parameter.
    inventory_label = input_file["global"]['inventory_label']
    inventory_dir = input_file["global"]['nc_emission_dir']
    inventory_year = input_file["global"]['inventory_year']
    
    #* Calculate the allocation factor.
    sectors = input_file["global"]['sectors']
    if type(sectors) != type(["list"]):
        sectors = [sectors]
    
    allocators = input_file["global"]['allocator']
    if type(allocators) != type(["list"]):
        allocators = [allocators]
    
    allocator_types = input_file["global"]['allocator_type']
    if type(allocator_types) != type(["list"]):
        allocator_types = [allocator_types]
    
    #* These is the control for each module.
    _ = input_file["control"]['control_model_grid']
    if _ == 'Y':
        control_model_grid = True
    else:
        control_model_grid = False

    _ = input_file["control"]['control_grid_info']
    if _ == 'Y':
        control_grid_info = True
    else:
        control_grid_info = False
        
    _ = input_file["control"]['control_coarse_emission']
    if _ == 'Y':
        control_coarse_emission = True
    else:
        control_coarse_emission = False
    
    _ = input_file["control"]['control_allocation_factor']
    if _ == 'Y':
        control_allocation_factor = True
    else:
        control_allocation_factor = False
    
    print("The setting of namelist.input:")
    print(f"gridf: {gridf}")
    print(f"gdnm: {gdnm}")
    print(f"zone_file: {zone_file}")
    print(f"inventory_label: {inventory_label}")
    print(f"inventory_dir: {inventory_dir}")
    print(f"inventory_year: {inventory_year}")
    print(f"start_date, end_date: {start_date}, {end_date}")
    print(f"sectors: {sectors}")
    print(f"allocator: {allocators}")
    print(f"allocator_type: {allocator_types}")
    print(f"Control model grid: {control_model_grid}")
    print(f"Control grid information: {control_grid_info}")
    print(f"Control coarse emission: {control_coarse_emission}")
    print(f"Control allocation factor: {control_allocation_factor}")
    print()
    
    date_list = pd.period_range(pd.to_datetime(start_date), pd.to_datetime(end_date), freq='D')
    mms = np.unique(np.array([temp_date.strftime("%m") for temp_date in date_list]))
    
    #TODO Draw the model grid.
    if control_model_grid:
        print('------------------------------------------')
        print('Start the draw model grid module')
        print('------------------------------------------')
        model_grid, grid_dx, grid_dy = draw_model_grid(gridf)
        # print('Successful for drawing model grid.')
        print()
    else:
        model_grid = f'Output/{gdnm}_grid.shp'
    
    #TODO Define the parent grid.
    if control_grid_info:
        print('------------------------------------------')
        print('Start the define parent grid module.')
        print('------------------------------------------')
        # print('Successful for defining parent grid.')
        print()
        grid_info = define_parent_grid(model_grid, zone_file, gdnm)
    else:
        grid_info = pd.read_csv(f'Output/{gdnm}_ParentGrid.csv')
    
    # Add longitude and latitude to this DataFrame.
    gf = pnc.pncopen(gridf, format='griddesc', GDNAM=gdnm, SDATE=1970001)
    grid_info['LON'], grid_info['LAT'] = gf.ij2ll(grid_info.colnum.values-1, grid_info.rownum.values-1)
    
    #TODO Coarse emission, this process will produce the coarse emission file.
    coarse_emis_dir = 'Output/CoarseEmission'
    os.makedirs(coarse_emis_dir, exist_ok=True)
    if control_coarse_emission:
        print('------------------------------------------')
        print('Start the coarse emission module')
        print('------------------------------------------')
        for mm in mms:
            for sector in sectors:
                for raster_file in glob.glob(f'{inventory_dir}/{inventory_label}_{inventory_year}_{mm}__{sector}__*.nc'):   
                    emission_df = zonal(raster_file, zone_file, gridf, gdnm, 'NAME')
                    sub_name = os.path.basename(raster_file)
                    (sub_name, suffix) = os.path.splitext(sub_name)
                    emission_df.to_csv(f'Output/CoarseEmission/{sub_name}.csv')
                    print(f'Successful calculating of coarse emission in {raster_file}.')
                    print()
    
    #TODO Calculate the allocation factor and this process will produce the allocation factor file.
    factor_dir = 'Output/factor'
    os.makedirs(factor_dir, exist_ok=True)
    if control_allocation_factor:
        print('------------------------------------------')
        print('Start the allocation factor module')
        print('------------------------------------------')
        #* Raster allocator.
        allocator_dir = 'Input/SA'
        for allocator_fn, sector, allocator_type in zip(allocators, sectors, allocator_types):
            allocator_file = f'{allocator_dir}/{allocator_fn}'
            factor_df = allocator_factor(allocator_file, sector, grid_info, allocator_type)
    
    #TODO Create source file.
    for mm in mms:
        source_dir = 'Output/source'
        os.makedirs(source_dir, exist_ok=True)
        create_source(inventory_label, inventory_year, mm, sectors, factor_dir, coarse_emis_dir, source_dir)
        
    #TODO Create model-ready file.
    inventory_mechanism = input_file["global"]['inventory_mechanism']
    target_mechanism = input_file["global"]['target_mechanism']
    for _date in date_list:
        for sector in sectors:
            _date = str(_date)
            create_emission_file(_date, gridf, gdnm, inventory_label, sector, inventory_year, inventory_mechanism, target_mechanism)
    
    
   
