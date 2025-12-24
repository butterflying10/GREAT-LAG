# =================================================================
# PPP画图主程序
# =================================================================
import logging
import os
from tqdm import tqdm

from PythonScripts.dependency.gnss_flt_io import readGreatFltFile_v2
from PythonScripts.snx_to_crd import readCrdFile
from PythonScripts.dependency.gnss_great_draw_flt import gnss_great_draw_compare_flt

siteList = [
    "JFNG"
]

# 数据的年份
year = 2024
# 数据的年积日
day = 122
# 设置数据类型：浮点解(ppp-float)或者固定解(ppp-fixed)
type1 = "JFNG GPS"
type2 = "JFNG BDS"

# 对应 type1 的测站数据
fltPath1 = rf"..\data_ppp\float_2024122"
# 对应 type2 的测站数据
fltPath2 = rf"..\data_ppp\float_2024122"
# 保存的结果目录
savePath = rf"..\output"
# gnss_crd_io.py 脚本输出的 crd 文件路径
crdPath = rf"..\data_ppp\crd\snx_igs_2024_122.crd"

crdData = readCrdFile(crdPath)
for doy in tqdm(range(day, day + 1)):
    yeardoy_savePath = os.path.join(savePath, f"{year:04d}_{doy:03d}_sta_test")
    for site in siteList:
        logging.info(f"{site}")
        for i in range(1):
            fltFilePath1 = os.path.join(fltPath1, f"{site}-SIMU_LEO_G.flt")
            fltFilePath2 = os.path.join(fltPath2, f"{site}-SIMU_LEO_C.flt")
            if not os.path.exists(fltFilePath2):
                continue
            fltData2 = readGreatFltFile_v2(fltFilePath2)
            if not os.path.exists(fltFilePath1):
                continue
            fltData1 = readGreatFltFile_v2(fltFilePath1)
            gnss_great_draw_compare_flt(site, type1, type2, fltData1, fltData2, crdData, yeardoy_savePath)
