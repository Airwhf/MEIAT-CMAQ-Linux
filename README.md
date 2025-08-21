# MEIAT-CMAQ Linux Version 5.0

1. The species profile from MEIC are archived `https://github.com/Airwhf/My_WRF_CMAQ_Configuration/tree/main/share`.

## 使用说明

### 1. 使用`f2c.py`将其他排放数据转换为CMAQ所需的排放文件格式。

在此目录下执行`python f2c.py`，即可将其他排放数据转换为CMAQ所需的排放文件格式。

运行之前请先设置：

```python
# ========================================================================================
# Set GRIDDESC configuration.
griddesc_file = 'input/GRIDDESC.CN27km'
griddesc_name = 'CN27km'    

# Set the inventory with Geotiff format.
emission_dir = r'/Volumes/project/Emissions/EDGAR_reformat'  
sectors = ['AGRICULTURE', 'BUILDINGS', 'FUEL_EXPLOITATION',
            'IND_COMBUSTION', 'IND_PROCESSES', 'POWER_INDUSTRY', 'TRANSPORT', 'WASTE']

# Set the inventory.
inventory_label = 'FT2022'
inventory_year = 2019

# Set the inventory period.
start_date = '2017-01-01'
end_date = '2017-01-02'

# Species allocation.
inventory_mechanism = 'EDGAR'
target_mechanism = 'CB06'

# shape factor
shapefactor = 2
# ========================================================================================
```




