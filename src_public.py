import datetime
import glob
import os.path
import numpy as np
import re
import xarray as xr

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



def emis_mix(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Search the files.
    files = glob.glob(f"{input_dir}/*.nc")

    # Get the name of pollutants.
    for file in tqdm.tqdm(files):
        ds = xr.open_dataset(file)
        lats = ds.coords["lat"].__array__()
        lons = ds.coords["lon"].__array__()
        lonmin, latmax, lonmax, latmin = lons.min(), lats.max(), lons.max(), lats.min()
        num_lon = lons.shape[0]
        num_lat = lats.shape[0]
        res = 0.25
        transform = Affine.translation(lonmin - res / 2, latmin - res / 2) * Affine.scale(res, res)

        # Set sectors.
        sectors = ["POWER", "INDUSTRY", "RESIDENTIAL", "TRANSPORT", "AGRICULTURE"]
        file_name = os.path.basename(file)
        condition = fr"MICS_Asia_(.*?)_{year}_0.25x0.25.nc"
        pollutant = re.findall(condition, file_name)[0]

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

                tiffile = f"MIX_{year}_{month}__{sector_label}__{pollutant}.tiff"
                with rio.open(f"{output_dir}/{tiffile}",
                              'w',
                              driver='GTiff',
                              height=num_lat,
                              width=num_lon,
                              count=1,
                              dtype=temp_var.dtype,
                              crs='+proj=latlong',
                              transform=transform, ) as dst:
                    dst.write(temp_var, 1)
                # print(f"Finish and output {output_dir}/{tiffile}.")

