import datetime
import glob
import os.path
import numpy as np
import re
import xarray as xr

def time_format():
    return f'{datetime.datetime.now()}|> '


def emis_meic(input_dir):
    output_dir = input_dir
    # os.makedirs(output_dir, exist_ok=True)
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

# if __name__ == '__main__':
#     input_dir = r'F:\data\Emission\MEICv1.4\2020'
#     emis_meic(input_dir)

