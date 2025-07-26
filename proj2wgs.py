file = 'shapefiles/CN-Province_WGS84.shp'


import geopandas as gpd

# Read the shapefile
gdf = gpd.read_file(file)

# Convert Chinese to English in the 'NAME' field using if/elif
def convert_name(name):
    mapping = {
        '黑龙江': 'Heilongjiang',
        '新疆': 'Xinjiang',
        '山西': 'Shanxi',
        '宁夏': 'Ningxia',
        '西藏': 'Tibet',
        '山东': 'Shandong',
        '河南': 'Henan',
        '江苏': 'Jiangsu',
        '安徽': 'Anhui',
        '湖北': 'Hubei',
        '浙江': 'Zhejiang',
        '江西': 'Jiangxi',
        '湖南': 'Hunan',
        '云南': 'Yunnan',
        '贵州': 'Guizhou',
        '福建': 'Fujian',
        '广西': 'Guangxi',
        '广东': 'Guangdong',
        '海南': 'Hainan',
        '吉林': 'Jilin',
        '辽宁': 'Liaoning',
        '天津': 'Tianjin',
        '青海': 'Qinghai',
        '甘肃': 'Gansu',
        '陕西': 'Shaanxi',
        '内蒙古': 'Inner Mongolia',
        '重庆': 'Chongqing',
        '河北': 'Hebei',
        '上海': 'Shanghai',
        '北京': 'Beijing',
        '台湾': 'Taiwan',
        '香港': 'Hong Kong',
        '澳门': 'Macau',
        '四川': 'Sichuan'
    }
    return mapping.get(name, name)

if 'NAME' in gdf.columns:
    gdf['NAME'] = gdf['NAME'].apply(convert_name)
else:
    print("No 'NAME' field found in the shapefile.")

# Convert to WGS84 projection (EPSG:4326)
gdf_wgs84 = gdf.to_crs(epsg=4326)

# Save the reprojected shapefile
gdf_wgs84.to_file('data/shapefiles/CN-Province_WGS84.shp')



