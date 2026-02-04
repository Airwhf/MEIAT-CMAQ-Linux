# MEIAT-CMAQ Linux Version 5.0

1. The species profile from MEIC are archived `https://github.com/Airwhf/My_WRF_CMAQ_Configuration/tree/main/share`.

## 使用说明

### 1. 融合处理MEIC和MIX排放清单

1. Run `./meic_f2c.py` for meic emissions.
2. Run `./reformat_mix.py` for original format mix files.
3. Run `./mask_emis_by_shp.py` for post-mix files. **Keep: method = 'remove'**
4. Run `./mix_f2c.py` for CMAQ-ready files.





