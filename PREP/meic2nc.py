#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Time      :2023/7/23 10:48
# @Author    :Haofan Wang

import glob
import os

import numpy as np
import xarray as xr
import re

if __name__ == "__main__":
    print("This script is written by Haofan Wang.")
    # ------------------------------------------
    input_dir = r"/mnt/f/data/MEICv1.4/2020"
    output_dir = r"/mnt/e/GitHub/MEIAT-CMAQ-Linux/Input/MEICv1.4-2020"
    # ------------------------------------------

    os.makedirs(output_dir, exist_ok=True)

    files = glob.glob(f"{input_dir}/*.asc")

    for file in files:
        print(file)
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
        z = np.array(np.where(_ == -9999.0, 0.0, _), dtype='float')

        # 最大最小经纬度
        min_long, min_lat, max_long, max_lat = 70.0, 10.0, 150.0, 60.0

        # 分辨率
        x_resolution = 0.25
        y_resolution = 0.25

        lons = np.arange(min_long, max_long, x_resolution) + 0.25
        lats = np.flip(np.arange(min_lat, max_lat, y_resolution)) - 0.25

        
        # 使用已知的经纬度和数值创建xarray数据数组
        # coords=[('lat', lats, {'units': 'degrees_north'}),
        #         ('lon', lons, {'units': 'degrees_east'})]
        
        # da = xr.DataArray(z, coords=coords, dims=['lat', 'lon'])

        da = xr.DataArray(
            z,
            coords=[('lat', lats, {'units': 'degrees_north'}),
                    ('lon', lons, {'units': 'degrees_east'})],
            attrs={'description': 'Example data'}
        )
        
        # 创建xarray数据集并将数据数组添加进去
        ds = xr.Dataset({'z': da})
        # 将 DataArray 转换为 xarray.Dataset
        dataset = xr.Dataset({'z': da})
        
        dataset.to_netcdf(output_name, format='NETCDF4')




