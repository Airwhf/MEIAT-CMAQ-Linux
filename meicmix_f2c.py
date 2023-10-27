# @Time    : 2023/02/21 19:17
# @Author  : Haofan Wang
# @Version : python3.9
# @Email   : wanghf58@mail2.sysu.edu.cn
import os
import time
from src_f2c import *
from src_public import *


if __name__ == "__main__":
    # ========================================================================================
    # Set GRIDDESC configuration.
    griddesc_file = 'input/GRIDDESC.CN27km'
    griddesc_name = 'CN27km'    
    
    # Set the inventory with Geotiff format.
    meic_asc_dir = r'F:\data\Emission\MEICv1.4\2020'  
    mix_dir = r'F:\data\Emission\MIX\2010'
    mix_year = 2010
    meic_year = 2020
    merge_output_dir = r'F:\data\Emission\MEICMIX\2020'
    sectors = ['residential', 'power', 'industry', 'agriculture', 'transportation']
    
    # Set the inventory.
    inventory_label = 'MEICMIX'
    inventory_year = meic_year 
    
    # Set the inventory period.
    start_date = '2017-01-01'
    end_date = '2017-01-02'
    
    # Species allocation.
    inventory_mechanism = 'MEIC-CB05'
    target_mechanism = 'CB06'
    
    # shape factor
    shapefactor = 4
    # ========================================================================================
    
    start_time = time.time()
    # process meic.
    meic_emission_dir = f'{meic_asc_dir}/MEIAT'
    os.makedirs(meic_emission_dir, exist_ok=True)
    emis_meic(meic_asc_dir, meic_emission_dir)  
    
    # process mix.
    mix_emission_dir = f'{mix_dir}/MEIAT'
    os.makedirs(mix_emission_dir, exist_ok=True)
    emis_mix(mix_dir, mix_emission_dir, mix_year)
    
    # merge meic and mix.
    emission_dir = merge_output_dir
    merge_mixmeic(meic_emission_dir, meic_year, mix_emission_dir, mix_year, emission_dir)
    
    output_dir = f'model_emission_{griddesc_name}'
    os.makedirs(output_dir, exist_ok=True)
    periods = pd.period_range(pd.to_datetime(start_date), pd.to_datetime(end_date), freq='D')
    for sector in sectors:
        for emission_date in periods:
            source2cmaq(str(emission_date), griddesc_file, griddesc_name, sector, emission_dir, inventory_mechanism, target_mechanism, output_dir, shapefactor=shapefactor)
            
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"### Time consuming: {elapsed_time} s ###")
