import datetime
import glob
import os.path
import numpy as np
import re
import xarray as xr
import shutil

def time_format():
    return f'{datetime.datetime.now()}|> '


def emis_meic(input_dir, output_dir):
    # output_dir = input_dir
    os.makedirs(output_dir, exist_ok=True)
    files = glob.glob(f"{input_dir}/*.asc")
    for file in files:
        sub_name = os.path.basename(file)
        condition = f"(.*?)_(.*?)_(.*?)_(.*?).asc"
        # condition = f"(.*?)_(.*?)__(.*?)__(.*?).asc"  # For 2019 and 2020.
        encode_name = re.findall(condition, sub_name)[0]
        year = r"%.4d" % int(encode_name[0])
        mm = r"%.2d" % int(encode_name[1])
        sector = encode_name[2]
        pollutant = encode_name[3].replace(".", "")
        output_name = f"{output_dir}/MEIC_{year}_{mm}__{sector}__{pollutant}.nc"

        # 以只读模式打开文件
        with open(file, 'r', encoding='utf-8') as file:
            # 逐行读取文件内容，并以空格为分隔符分割每行，最后将其添加到数组（列表）中
            lines = [line.strip().split() for index, line in enumerate(file) if index >= 6]

        # 打印读取到的数组（列表）
        _ = np.array(lines, dtype="float")
        z = _
        # z = np.where(_ == -9999.0, 0.0, _)

        # 最大最小经纬度
        min_long, min_lat, max_long, max_lat = 70.0, 10.0, 150.0, 60.0

        # 分辨率
        x_resolution = 0.25
        y_resolution = 0.25

        # 计算栅格的行和列
        width = int((max_long - min_long) / x_resolution)
        height = int((max_lat - min_lat) / y_resolution)

        # Create netcdf file via xarray.
        # 创建一个数据数组 (DataArray) 从你的 z 数组
        data_array = xr.DataArray(z, dims=['latitude', 'longitude'])

        # 添加坐标信息
        data_array = data_array.assign_coords(longitude=np.linspace(min_long, max_long, width), latitude=np.linspace(max_lat, min_lat, height))

        # 保存数据数组为一个数据集 (Dataset)
        dataset = data_array.to_dataset(name='emission')

        # 添加坐标信息
        dataset.longitude.attrs['units'] = 'degrees_east'
        dataset.latitude.attrs['units'] = 'degrees_north'

        # 保存数据集为 NetCDF 文件
        dataset.to_netcdf(output_name)
        
        print(f"{time_format()} {output_name} has been created.")



def emis_mix(input_dir, output_dir, year):
    os.makedirs(output_dir, exist_ok=True)

    # Search the files.
    files = glob.glob(f"{input_dir}/*.nc")

    # Get the name of pollutants.
    for file in files:
        ds = xr.open_dataset(file)
        lats = ds.coords["lat"].__array__()
        lons = ds.coords["lon"].__array__()
        lonmin, latmax, lonmax, latmin = lons.min(), lats.max(), lons.max(), lats.min()
        num_lon = lons.shape[0]
        num_lat = lats.shape[0]
        res = 0.25

        # Set sectors.
        sectors = ["POWER", "INDUSTRY", "RESIDENTIAL", "TRANSPORT", "AGRICULTURE"]
        file_name = os.path.basename(file)
        condition = fr"MICS_Asia_(.*?)_{year}_0.25x0.25.nc"
        pollutant = re.findall(condition, file_name)[0].replace(".", "")

        for var_name in list(ds.keys()):
            var = ds[var_name]
            months = var.__getattr__("time").values
            for i in range(12):
                month = "%.2d" % months[i]
                temp_var = var[i, ...].values
                sector = var_name.split("_")[-1]

                if sector == "TRANSPORT":
                    sector_label = "transportation"
                else:
                    sector_label = sector.lower()

                output_file = f"{output_dir}/MIX_{year}_{month}__{sector_label}__{pollutant}.nc"
                # 创建一个数据数组 (DataArray) 从你的 z 数组
                data_array = xr.DataArray(temp_var, dims=['latitude', 'longitude'])
                # 添加坐标信息
                data_array = data_array.assign_coords(longitude=lons, latitude=lats)

                # 保存数据数组为一个数据集 (Dataset)
                dataset = data_array.to_dataset(name='emission')

                # 添加坐标信息
                dataset.longitude.attrs['units'] = 'degrees_east'
                dataset.latitude.attrs['units'] = 'degrees_north'

                # 保存数据集为 NetCDF 文件
                dataset.to_netcdf(output_file)
                
                print(f"{time_format()} {output_file} has been created.")

def mosaic(mix_file, meic_file, output_file, missing_value=-9999):
    # 打开东亚地区的排放文件和中国地区的排放文件
    east_asia_file = xr.open_dataset(mix_file)
    china_file = xr.open_dataset(meic_file)
    
    # 确保坐标信息匹配，对中国地区的数据进行坐标对齐
    china_file = china_file.reindex({'longitude': east_asia_file['longitude'], 'latitude': east_asia_file['latitude']}, method='nearest', fill_value=missing_value)
    # print(china_file['emission'].values)
    
    # 识别中国地区的部分（假设中国部分的无效值为missing_value）
    china_mask = (china_file['emission'] == missing_value)
    # print(china_mask)

    # 从中国地区的排放文件中提取中国部分的数据
    china_emissions = china_file['emission']

    # 使用中国地区数据填充东亚地区的数据
    east_asia_file['emission'].values = np.where(~china_mask, china_emissions, east_asia_file['emission'].values)
    # print(east_asia_file['emission'])

    # 保存填充后的数据为新的 netCDF 文件
    east_asia_file.to_netcdf(output_file)

    # 关闭文件
    east_asia_file.close()
    china_file.close()
    print(f"{time_format()} {output_file} has been created.")

def merge_mixmeic(meic_dir, meic_year, mix_dir, mix_year, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    meic_files = glob.glob(f"{meic_dir}/*_{meic_year}_*.nc")
    template_mix_file = glob.glob(f"{mix_dir}/*_{mix_year}_*.nc")[0]
    template_mix_dataset = xr.open_dataset(template_mix_file)
    for meic_file in meic_files:
        basename = os.path.basename(meic_file)
        condition = f"MEIC_{meic_year}_(.*?)__(.*?)__(.*?).nc"
        encode_name = re.findall(condition, basename)[0]
        mm = encode_name[0]
        sector = encode_name[1]
        pollutant = encode_name[2]
        
        mix_file = f"{mix_dir}/MIX_{mix_year}_{mm}__{sector}__{pollutant}.nc"
        if os.path.exists(mix_file) is False:
            # copy a mix file.
            template_mix_dataset['emission'].values = template_mix_dataset['emission'].values * 0.0
            template_mix_dataset.to_netcdf(mix_file)
            template_mix_dataset.close()
        
        output_file = f'{output_dir}/MEICMIX_{meic_year}_{mm}__{sector}__{pollutant}.nc'
        mosaic(mix_file, meic_file, output_file)
        # print(f"{time_format()} {output_file} has been merged.")
            
    pass
    
if __name__ == '__main__':
    mix_file = r'F:\data\Emission\MIX\2010\MEIAT\MIX_2010_01__agriculture__NH3.nc'
    meic_file = r'F:\data\Emission\MEICv1.4\2020\MEIAT\MEIC_2020_01__agriculture__NH3.nc'
    output_file = 'MEICMIX_2010_01__agriculture__NH3.nc'
    # mosaic(mix_file, meic_file, output_file)
    mix_dir, mix_year = r'F:\data\Emission\MIX\2010\MEIAT', 2010
    meic_dir, meic_year = r'F:\data\Emission\MEICv1.4\2020\MEIAT', 2020
    output_dir = r'F:\data\Emission\MEICMIX\2020\MEIAT'
    merge_mixmeic(meic_dir, meic_year, mix_dir, mix_year, output_dir)
