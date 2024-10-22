#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName  :reformat_edgar.py
# @Time      :2024/10/21 23:15:03
# @Author    :Haofan Wang

# 此数据用于将EDGAR数据库中的排放清单数据转换为MEIAT-CMAQ所需的格式：
#
# 原始的排放清单下载地址如下：
#   https://edgar.jrc.ec.europa.eu/dataset_ap81#p4
# 在此链接中下载`bkl_{sector}_emi_nc.zip`数据，该数据包含气体和气溶胶的排放清单。
# 下载 `Monthly sector-specific gridmaps (2000-2022)`

# 该数据一共包含以下几个部门：
#   AGRICULTURE
#   BUILDINGS
#   FUEL_EXPLOITATION
#   IND_COMBUSTION
#   IND_PROCESSES
#   POWER_INDUSTRY
#   TRANSPORT
#   WASTE

import os
import xarray as xr

if __name__ == "__main__":
    
    # 设置输入目录
    input_dir = r'/Volumes/project/Download/EDGAR'
    
    # 设置输出目录
    output_dir = r'/Volumes/project/Emissions/EDGAR_reformat'
    os.makedirs(output_dir, exist_ok=True)
    
    # 设置部门列表以及月分配系数（月分配系用于处理VOC）
    sector_dict = {
        'AGRICULTURE': [0.06, 0.06, 0.07, 0.08, 0.10, 0.12, 0.12, 0.12, 0.10, 0.08, 0.06, 0.03],  # Higher activity in summer months
        'BUILDINGS': [0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.05, 0.05, 0.06, 0.07, 0.09, 0.13],  # Higher in winter due to heating
        'FUEL_EXPLOITATION': [0.08, 0.08, 0.08, 0.08, 0.08, 0.09, 0.10, 0.10, 0.09, 0.08, 0.07, 0.07],  # Slight seasonal variation
        'IND_COMBUSTION': [0.08, 0.08, 0.08, 0.08, 0.08, 0.09, 0.10, 0.10, 0.09, 0.08, 0.07, 0.07],  # Similar to fuel exploitation
        'IND_PROCESSES': [0.08, 0.08, 0.08, 0.08, 0.08, 0.09, 0.10, 0.10, 0.09, 0.08, 0.07, 0.07],  # Consistent throughout the year
        'POWER_INDUSTRY': [0.12, 0.11, 0.10, 0.09, 0.08, 0.06, 0.05, 0.05, 0.06, 0.08, 0.10, 0.10],  # Higher demand in winter
        'TRANSPORT': [0.09, 0.08, 0.08, 0.08, 0.08, 0.09, 0.09, 0.09, 0.09, 0.09, 0.07, 0.07],  # Fairly consistent throughout the year
        'WASTE': [0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.09, 0.09, 0.09, 0.08, 0.08, 0.07]  # Small seasonal variation
    }
        
    # 设置输入物种名称
    species_in = [
        'AP_BC', 'AP_CO', 'AP_NH3', 'AP_NOx', 'AP_OC',
        'AP_PM2.5', 'AP_PM10', 'AP_SO2', 'AP_NMVOC', 
        'VOC_spec_voc1', 'VOC_spec_voc2', 'VOC_spec_voc3',
        'VOC_spec_voc4', 'VOC_spec_voc5', 'VOC_spec_voc6',
        'VOC_spec_voc7', 'VOC_spec_voc8', 'VOC_spec_voc9',
        'VOC_spec_voc10', 'VOC_spec_voc11', 'VOC_spec_voc12',
        'VOC_spec_voc13', 'VOC_spec_voc14', 'VOC_spec_voc15',
        'VOC_spec_voc16', 'VOC_spec_voc17', 'VOC_spec_voc18',
        'VOC_spec_voc19', 'VOC_spec_voc20', 'VOC_spec_voc21',
        'VOC_spec_voc22', 'VOC_spec_voc23', 'VOC_spec_voc24',
        'VOC_spec_voc25', 
    ]
    
    # EDGAR VOC 转 CB06 物种分配表
    voc_species_map = {
        'voc1': 'CB06_ETOH', # 46.2
        'voc2': 'CB06_ETHA', # 30
        'voc3': 'CB06_PRPA', # 44
        'voc4': '-',
        'voc5': '-',
        'voc6': 'CB06_SOAALK', # 112
        'voc7': 'CB06_ETH', # 28
        'voc8': '-',
        'voc9': 'CB06_ETHY', # 26
        'voc10': 'CB06_ISOP', # 118.1
        'voc11': 'CB06_TERP', # 136.2
        'voc12': '-',
        'voc13': 'CB06_BENZ', # 78.1
        'voc14': 'CB06_TOL', # 92.1
        'voc15': 'CB06_XYLMN', # 106.2
        'voc16': '-',
        'voc17': '-',
        'voc18': '-',
        'voc19': '-',
        'voc20': '-',
        'voc21': 'CB06_ECH4', # 16
        'voc22': 'CB06_ALD2', # 44
        'voc23': 'CB06_KET', # 72.1
        'voc24': 'CB06_FACD', # 46
        'voc25': '-',
    }
    
    # 设置需要处理的年份
    year = 2019
    
    # 遍历部门和物种，将其转换为MEIAT-CMAQ格式
    for sector in sector_dict.keys():
        for species in species_in:
            
            # 处理VOC排放
            if 'VOC_spec_' in species:
                voc_number = species.split('_')[-1]
                
                # 读取对应的物种分配表
                voc_speice = voc_species_map[voc_number]
                if voc_speice == '-':
                    continue
                
                src = f'{input_dir}/v8.1_FT2022_{species}_{year}_bkl_{sector}_emi.nc'
                if not os.path.exists(src):
                    print('文件不存在：', src)
                    continue
                
                # 获取排放数据和经纬度信息
                ds = xr.open_dataset(src)['emissions']
                lons = ds.coords['lon'].values
                lats = ds.coords['lat'].values
                
                for month_i in range(12):
                    # 获取对应月份的数据
                    month = f'{month_i + 1:02d}'
                    temp_var = ds.values * sector_dict[sector][month_i]  # 进行月分配 
                    
                    # 设置输出文件名
                    dst = f'{output_dir}/FT2022_{year}_{month}__{sector}__{voc_speice}.nc'
                    
                    if os.path.exists(dst):
                        continue
                    
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
                    dataset.to_netcdf(dst)
                    
                    # print(f'{dst} saved!')
            
            # 处理其他污染物排放
            if 'AP_' in species:
                
                if 'PM10' in species: # 计算PMC
                    pm10_src = f'{input_dir}/v8.1_FT2022_AP_PM10_{year}_bkl_{sector}_emi.nc'
                    pm25_src = f'{input_dir}/v8.1_FT2022_AP_PM2.5_{year}_bkl_{sector}_emi.nc'
                    
                    if not os.path.exists(pm10_src) or not os.path.exists(pm25_src):
                        continue
                    
                    # 获取排放数据和经纬度信息
                    ds = xr.open_dataset(pm10_src)['emissions']
                    lons = ds.coords['lon'].values
                    lats = ds.coords['lat'].values
                    
                    for month_i in range(12):
                        
                        # 获取对应月份的数据
                        month = f'{month_i + 1:02d}'
                        
                        pm10_ds = xr.open_dataset(pm10_src)['emissions'][month_i, ...].values
                        pm25_ds = xr.open_dataset(pm25_src)['emissions'][month_i, ...].values
                        temp_var = pm10_ds - pm25_ds
                        
                        # 设置输出文件名
                        dst = f'{output_dir}/FT2022_{year}_{month}__{sector}__PMC.nc'
                        
                        # 如果文件已经存在 则直接跳过
                        if os.path.exists(dst):
                            continue
                        
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
                        dataset.to_netcdf(dst)
                    
                    
                else:
                    # continue  # 循环调试
                    src = f'{input_dir}/v8.1_FT2022_{species}_{year}_bkl_{sector}_emi.nc'
                    if not os.path.exists(src):
                        print('文件不存在：', src)
                        continue
                    # print(src)
                    # MEIC_{year}_{mm}__{sector}__{pollutant}.nc
                    
                    # 获取排放数据和经纬度信息
                    ds = xr.open_dataset(src)['emissions']
                    lons = ds.coords['lon'].values
                    lats = ds.coords['lat'].values
                    
                    for month_i in range(12):
                        
                        # 获取对应月份的数据
                        month = f'{month_i + 1:02d}'
                        temp_var = ds[month_i, ...].values
                        
                        # 设置输出文件名
                        dst = f'{output_dir}/FT2022_{year}_{month}__{sector}__{species[3::]}.nc'
                        
                        if os.path.exists(dst):
                            continue
                        
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
                        dataset.to_netcdf(dst)
                        
                        # print(f'{dst} saved!')
                
                
                
                
        
        