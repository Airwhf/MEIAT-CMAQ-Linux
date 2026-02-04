import os
import sys
import time
import pandas as pd
from multiprocessing import Pool, cpu_count
from functools import partial
from src_f2c import *

sys.stdout = open(os.devnull, 'w')  # 可选：关闭标准输出

# 单个任务执行函数
def run_source2cmaq(sector, emission_date, griddesc_file, griddesc_name, emission_dir,
                    inventory_mechanism, target_mechanism, output_dir, shapefactor):
    source2cmaq(
        str(emission_date), griddesc_file, griddesc_name,
        sector, emission_dir, inventory_mechanism,
        target_mechanism, output_dir, shapefactor=shapefactor
    )

if __name__ == "__main__":
    # ========================================================================================
    griddesc_file = 'input/GRIDDESC.CN30KM'
    griddesc_name = 'CN30KM'

    emission_dir = r'/work/home/wanghf58/MEIAT-CMAQ-Linux-main/emissions/MIX_MEIAT_masked_remove'
    sectors = ['residential', 'power', 'industry', 'agriculture', 'transportation']
    inventory_label = 'MIX'
    inventory_year = 2017
    start_date = '2022-12-19'
    end_date = '2024-01-02'
    inventory_mechanism = 'MEIC-CB05'
    target_mechanism = 'CB06'
    shapefactor = 2
    core_nums = 64
    # ========================================================================================

    start_time = time.time()
    output_dir = f'model_emission_{griddesc_name}'
    os.makedirs(output_dir, exist_ok=True)

    periods = pd.period_range(pd.to_datetime(start_date), pd.to_datetime(end_date), freq='D')

    # 构建任务列表
    tasks = [
        (sector, emission_date)
        for sector in sectors for emission_date in periods
    ]

    # 定义参数函数（使用 partial 绑定不变参数）
    task_func = partial(
        run_source2cmaq,
        griddesc_file=griddesc_file,
        griddesc_name=griddesc_name,
        emission_dir=emission_dir,
        inventory_mechanism=inventory_mechanism,
        target_mechanism=target_mechanism,
        output_dir=output_dir,
        shapefactor=shapefactor
    )

    # 启动进程池
    with Pool(processes=core_nums) as pool:
        pool.starmap(task_func, tasks)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"### Time consuming (multiprocessing): {elapsed_time:.2f} s ###")

