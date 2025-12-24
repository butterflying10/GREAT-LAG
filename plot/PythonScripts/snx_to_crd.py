# =====================================================================================================
# sinex文件转crd文件，读写crd文件，并将其中的测站名称，测站天线类型，测站精密坐标存储下来
# =====================================================================================================

import logging
import os
import sys


pyDirName, pyFileName = os.path.split(os.path.abspath(__file__))
sys.path.append(f"{pyDirName}/../")

from PythonScripts.dependency.gnss_file_tool import isFileExist
from PythonScripts.dependency.gnss_sinex_io import readSinexFile, t_sinex
from PythonScripts.dependency.gnss_timestran_tool import doy2gpswd

# ==================================================================
# 读取并存储crd文件
# ==================================================================
def readCrdFile(crd_path):
    """
    读取crd文件，并返回一个索引dict：posData[prn][gpst]
    """
    if not isFileExist(crd_path):
        logging.error(f"Check your path : {crd_path}")
    fr = open(crd_path)
    lines = fr.readlines()
    fr.close()

    # 数据类
    dataSinex = t_sinex()

    # 数据行
    for line in lines:
        # 跳过备注行
        if line[0] == "*" or line[0] == "<":
            continue

        # 解码
        site = str.strip(line[0:5])
        site_stax = str.strip(line[6:27])
        site_stay = str.strip(line[28:49])
        site_staz = str.strip(line[50:71])
        site_stadx = str.strip(line[72:83])
        site_stady = str.strip(line[84:95])
        site_stadz = str.strip(line[96:107])
        site_longitude = str.strip(line[108:119])
        site_latitude  = str.strip(line[120:131])
        site_height = str.strip(line[132:139])
        site_antenna    = str.strip(line[140:160])
        site_receiver    = str.strip(line[161:181])

        # 存储
        if site not in dataSinex.site_list:
            dataSinex.site_list.add(site)
            dataSinex.site_stax[site] = site_stax
            dataSinex.site_stay[site] = site_stay
            dataSinex.site_staz[site] = site_staz
            dataSinex.site_stadx[site] = site_stadx
            dataSinex.site_stady[site] = site_stady
            dataSinex.site_stadz[site] = site_stadz
            dataSinex.site_longitude[site] = site_longitude
            dataSinex.site_latitude[site] = site_latitude
            dataSinex.site_height[site] = site_height
            dataSinex.site_antenna[site] = site_antenna
            dataSinex.site_receiver[site] = site_receiver

    return dataSinex

# ==================================================================
# 编码并存储crd文件
# crd_data : t_sinex 数据存储类
# crd_path : str()   文件存放路径
# ==================================================================
def WriteCrdFile(crd_data, crd_path):
    """
    编写crd文件
    """

    # 检查路径
    file_name = os.path.basename(crd_path)
    file_path = os.path.dirname(crd_path)
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    # 考虑到数据量不是很大，全部存储为list
    writeLines = list()

    # 文件头 
    headLine = "*CODE ____APRIORI_VALUE____ ____APRIORI_VALUE____ ____APRIORI_VALUE____ __STD_DEV__ __STD_DEV__ __STD_DEV__ _LONGITUDE_ _LATITUDE__ HEIGHT_ ___RECEIVER_TYPE____ ____ANTENNA_TYPE____\n"
    writeLines.append(headLine)

    # 循环测站写入
    # *CODE ____APRIORI_VALUE____ ____APRIORI_VALUE____ ____APRIORI_VALUE____ __STD_DEV__ __STD_DEV__ __STD_DEV__ _LONGITUDE_ _LATITUDE__ HEIGHT_ ___RECEIVER_TYPE____ ____ANTENNA_TYPE____ 
    writeData = crd_data
    emptyData = " "
    site_list = list(writeData.site_list)
    site_list.sort()
    for site in site_list:
        # 跳过这种无效测站
        if site == "----":
            continue
        # CODE
        writeLine = str(f"{site:>5s} ")
        # ____APRIORI_VALUE____ X
        if site in writeData.site_stax:
            writeLine = writeLine + f"{writeData.site_stax[site]:>21.6f} "
        else:
            writeLine = writeLine + f"{emptyData:>21s} "
        # ____APRIORI_VALUE____ Y
        if site in writeData.site_stay:
            writeLine = writeLine + f"{writeData.site_stay[site]:>21.6f} "
        else:
            writeLine = writeLine + f"{emptyData:>21s} "
        # ____APRIORI_VALUE____ Z
        if site in writeData.site_staz:
            writeLine = writeLine + f"{writeData.site_staz[site]:>21.6f} "
        else:
            writeLine = writeLine + f"{emptyData:>21s} "
        # __STD_DEV__ X
        if site in writeData.site_stadx:
            writeLine = writeLine + f"{writeData.site_stadx[site]:>11.4f} "
        else:
            writeLine = writeLine + f"{emptyData:>11s} "
        # __STD_DEV__ Y
        if site in writeData.site_stady:
            writeLine = writeLine + f"{writeData.site_stady[site]:>11.4f} "
        else:
            writeLine = writeLine + f"{emptyData:>11s} "
        # __STD_DEV__ Z
        if site in writeData.site_stadz:
            writeLine = writeLine + f"{writeData.site_stadz[site]:>11.4f} "
        else:
            writeLine = writeLine + f"{emptyData:>11s} "
        # _LONGITUDE_
        if site in writeData.site_longitude:
            writeLine = writeLine + f"{writeData.site_longitude[site]:>11.6f} "
        else:
            writeLine = writeLine + f"{emptyData:>11s} "
        # _LATITUDE__
        if site in writeData.site_latitude:
            writeLine = writeLine + f"{writeData.site_latitude[site]:>11.6f} "
        else:
            writeLine = writeLine + f"{emptyData:>11s} "
        # HEIGHT_
        if site in writeData.site_height:
            writeLine = writeLine + f"{writeData.site_height[site]:>7.1f} "
        else:
            writeLine = writeLine + f"{emptyData:>7s} "
        # ___RECEIVER_TYPE____
        if site in writeData.site_receiver:
            writeLine = writeLine + f"{writeData.site_receiver[site]:>20s} "
        else:
            writeLine = writeLine + f"{emptyData:>20s} "
        # ____ANTENNA_TYPE____
        if site in writeData.site_antenna:
            writeLine = writeLine + f"{writeData.site_antenna[site]:>20s} "
        else:
            writeLine = writeLine + f"{emptyData:>20s} "            
        writeLine = writeLine + "\n"
        writeLines.append(writeLine)

    # XML 文件部分
    for site in site_list:
        # 跳过这种无效测站
        if site == "----":
            continue

        # CODE
        writeLine = str(f"<rec id=\"{site}\" obj=\"SNX\" ")
        # ____APRIORI_VALUE____ X
        if site in writeData.site_stax:
            writeLine = writeLine + f"X=\"{writeData.site_stax[site]:>21.6f}\" "

        # ____APRIORI_VALUE____ Y
        if site in writeData.site_stay:
            writeLine = writeLine + f"Y=\"{writeData.site_stay[site]:>21.6f}\" "

        # ____APRIORI_VALUE____ Z
        if site in writeData.site_staz:
            writeLine = writeLine + f"Z=\"{writeData.site_staz[site]:>21.6f}\" "

        # __STD_DEV__ X
        if site in writeData.site_stadx:
            writeLine = writeLine + f"dX=\"{writeData.site_stadx[site]:>11.4f}\" "

        # __STD_DEV__ Y
        if site in writeData.site_stady:
            writeLine = writeLine + f"dY=\"{writeData.site_stady[site]:>11.4f}\" "

        # __STD_DEV__ Z
        if site in writeData.site_stadz:
            writeLine = writeLine + f"dZ=\"{writeData.site_stadz[site]:>11.4f}\" "

        # ___RECEIVER_TYPE____
        if site in writeData.site_receiver:
            writeLine = writeLine + f"rec=\"{writeData.site_receiver[site]:>20s}\" "

        # ____ANTENNA_TYPE____
        if site in writeData.site_antenna:
            writeLine = writeLine + f"ant=\"{writeData.site_antenna[site]:>20s}\" "

        writeLine = writeLine + f"/> \n"
        writeLines.append(writeLine)

    # 写入文件
    fw = open(crd_path, "w+")
    fw.writelines(writeLines)
    fw.close()

if __name__ == '__main__':
    
    logging.basicConfig(level=logging.DEBUG)
    sinex_data = None
    # 年
    year = 2024
    # 年积日
    day = 122
    # snx 文件夹路径
    sinex_path = rf"..\data_ppp\snx"
    # crd 文件夹路径（输出目录）
    crd_path = rf"..\data_ppp\crd"

    for doy in range(day, day+1):
        sinex_file = os.path.join(sinex_path, f"IGS0OPSSNX_{year:04d}{doy:03d}0000_01D_01D_SOL.SNX")
        crd_file = os.path.join(crd_path, f"snx_igs_{year:04d}_{doy:03d}.crd")

        if not os.path.exists(sinex_file):
            gpsw, gpswk = doy2gpswd(2024, doy)
            sinex_file = os.path.join(sinex_path, f"igs{str(year)[2:4]}P{gpsw:04d}{gpswk:01d}.snx")
        if not os.path.exists(sinex_file):
            continue

        logging.info(f"Read Sinex File {sinex_file}")
        sinex_data = readSinexFile(sinex_file, sinex_data)

        WriteCrdFile(sinex_data, crd_file)
        sinex_data = None
    exit(0)
