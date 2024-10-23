import os
import time
import pandas as pd
from multiprocessing import Pool
from src_f2c import *

# 定义函数来处理每个任务
def process_emission(args):
    sector, emission_date, griddesc_file, griddesc_name, emission_dir, inventory_mechanism, target_mechanism, output_dir, shapefactor = args
    source2cmaq(str(emission_date), griddesc_file, griddesc_name, sector, emission_dir, inventory_mechanism, target_mechanism, output_dir, shapefactor=shapefactor)

if __name__ == "__main__":
    # ========================================================================================
    # Set GRIDDESC configuration.
    griddesc_file = 'input/GRIDDESC.UK80x80'
    griddesc_name = 'UK80x80'    
    
    # Set the inventory with Geotiff format.
    emission_dir = r'/Volumes/project/Emissions/EDGAR_reformat'  
    sectors = ['AGRICULTURE', 'BUILDINGS', 'FUEL_EXPLOITATION',
               'IND_COMBUSTION', 'IND_PROCESSES', 'POWER_INDUSTRY', 'TRANSPORT', 'WASTE']
    
    # Set the inventory.
    inventory_label = 'FT2022'
    inventory_year = 2019
    
    # Set the inventory period.
    start_date = '2018-12-19'
    end_date = '2020-01-01'
    
    # Species allocation.
    inventory_mechanism = 'EDGAR'
    target_mechanism = 'CB06'
    
    # shape factor
    shapefactor = 2
    
    # 根据系统的CPU数量选择合适的进程数
    num_processes = 4  
    # ========================================================================================
    
    start_time = time.time()
    output_dir = f'model_emission_{griddesc_name}'
    os.makedirs(output_dir, exist_ok=True)
    
    # 生成日期周期
    periods = pd.period_range(pd.to_datetime(start_date), pd.to_datetime(end_date), freq='D')
    
    # 创建参数列表
    task_list = []
    for sector in sectors:
        for emission_date in periods:
            task_list.append((sector, emission_date, griddesc_file, griddesc_name, emission_dir, inventory_mechanism, target_mechanism, output_dir, shapefactor))
    
    # 使用多进程池进行并行处理
    with Pool(processes=num_processes) as pool:
        pool.map(process_emission, task_list)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print()
    print(f"### Time consuming: {elapsed_time} s ###")
