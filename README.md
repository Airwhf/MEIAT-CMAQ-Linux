# MEIAT-CMAQ For Linux



创建环境：

```shell
conda create -n meiat python=3.9
```

进入环境：

```shell
conda activate meiat
```

安装第三方库文件：

```shell
pip install -r requirements.txt
pip install PseudoNetCDF
```

以MEIC数据为例进行演示：

1. 进入`PREP`目录，打开`meic2nc.py`文件，修改输入和输出目录，并运行此程序。

```shell
python ./meic2nc.py
```

配置完成`nameilist.input`以后运行，测试运行过程中，只需要将`nc_emission_dir`修改为上述步骤中的输出目录即可。

```shell
python c2f.py
```


