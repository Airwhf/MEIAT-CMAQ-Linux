#!/bin/bash

# 基础URL
base_url="https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v81_FT2022_VOC_spec"

# 部门列表
departments=(
  "bkl_POWER_INDUSTRY"
  "bkl_FUEL_EXPLOITATION"
  "bkl_IND_COMBUSTION"
  "bkl_TRANSPORT"
  "bkl_BUILDINGS"
  "bkl_IND_PROCESSES"
  "bkl_AGRICULTURE"
  "bkl_WASTE"
)

# 循环下载 voc1 到 voc25 的文件
for i in {1..25}
do
  for dept in "${departments[@]}"
  do
    # 生成文件下载链接
    file_url="${base_url}/voc${i}/${dept}/${dept}_emi_nc.zip"
    
    # 生成保存的文件名（vocX_部门名称_emi_nc.zip）
    output_file="voc${i}_${dept}_emi_nc.zip"
    
    # 使用 wget 下载文件并保存为指定的文件名
    echo "正在下载: $file_url 保存为: $output_file"
    wget -O "$output_file" "$file_url"
  done
done

