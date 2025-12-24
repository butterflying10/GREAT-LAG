# -*- coding: utf-8 -*-
# ===================================================================
# 读取GREAT的flt格式
# ===================================================================

import logging
import os
import sys

pyDirName, pyFileName = os.path.split(os.path.abspath(__file__))
sys.path.append(f"{pyDirName}/../")

from PythonScripts.dependency.gnss_file_tool import isFileExist
from PythonScripts.dependency.gnss_math_tool import get_round
from PythonScripts.dependency.gnss_time_tool import epoch2gpst
from PythonScripts.dependency.gnss_pos_io import t_great_pos_data


# ==================================================================
# 读取并存储flt文件
# ==================================================================
def readGreatFltFile(posFilePath):
    """
    读取Pos文件，并返回一个索引dict：posData[gpsweek][gpssecd]
    """
    if not isFileExist(posFilePath):
        logging.error(f"Check your path : {posFilePath}")
        return False

    fr = open(posFilePath)
    lines = fr.readlines()
    fr.close()

    intv = 0
    posData = dict()
    epoNum = 0
    # 用于计算采样间隔
    preSec = 0
    nowSec = 0
    for line in lines:
        if "%" in line:
            continue

        data = t_great_pos_data()
        tmpData = ' '.join(line.split()).split()
        gpsweek = str()
        gpssecd = str()

# =EPO  2023-03-25 01:48:30.000      524910.0000       6510.0000  -2279829.1214   5004706.3738   3219777.3305         0.0077         0.0029         0.0022        -9.2136         0.0415              8         0.8837         0.0479         0.0096          Float
        # utc time
        if "-" in tmpData[1]:
            ymd = tmpData[1].split("-")
            hms = tmpData[2].split(":")
            igpsweek,igpssecd =  epoch2gpst([int(ymd[0]), int(ymd[1]), int(ymd[2]), int(hms[0]), int(hms[1]), float(hms[2])])
            gpsweek = f"{igpsweek:04d}"
            gpssecd = get_round(float(igpssecd),"0.01")
            gpssecd = f"{gpssecd:.3f}"
        else:
            continue

        intGsod = int(float(gpssecd))
        intGgpw = int(gpsweek)
        if intGgpw not in posData:
            posData[intGgpw] = dict()

        # 计算采样间隔
        preSec = float(nowSec)
        nowSec = float(gpssecd)
        if intv == 0:
            intv = nowSec - preSec
        else:
            if intv > (nowSec - preSec):
                intv = nowSec - preSec

# =EPO  2023-03-25 01:48:30.000      524910.0000       6510.0000  -2279829.1214   5004706.3738   3219777.3305         0.0077         0.0029         0.0022        -9.2136         0.0415              8         0.8837         0.0479         0.0096          Float
        data.gpst = gpsweek + " " + gpssecd
        data.x = float(tmpData[5])
        data.y = float(tmpData[6])
        data.z = float(tmpData[7])
        data.Q = int(0)
        data.ns = int(tmpData[13])
        data.sdx = float(tmpData[8])
        data.sdy = float(tmpData[9])
        data.sdz = float(tmpData[10])
        data.e = float(0)
        data.n = float(0)
        data.u = float(0)
        #data.dataline = line

        if "Float" == tmpData[17]:
            data.Q = 2
        if "Fixed" == tmpData[17]:
            data.Q = 1

        if intGsod not in posData[intGgpw]:
            posData[intGgpw][intGsod] = data
            epoNum = epoNum + 1
        #print(line)
    logging.info(f"Read flt file : {posFilePath}, epo number is {epoNum}")
    return intv, posData


# ==================================================================
# 读取并存储flt文件, great-2.0版本的格式
# ==================================================================
def readGreatFltFile_v2(posFilePath):
    """
    读取Pos文件，并返回一个索引dict：posData[gpsweek][gpssecd]
    """
    if not isFileExist(posFilePath):
        logging.error(f"Check your path : {posFilePath}")
        return False

    fr = open(posFilePath)
    lines = fr.readlines()
    fr.close()

    posData = dict()
    epoNum = 0

    for line in lines:
        if "#" in line:
            continue

        data = t_great_pos_data()
        tmpData = ' '.join(line.split()).split()

        intGsod = int(float(tmpData[0]))
        intGgpw = int(0)
        if intGgpw not in posData:
            posData[intGgpw] = dict()

        data.gpst = f"{0:04d}" + " " + "{intGsod:06d}"
        data.x = float(tmpData[1])
        data.y = float(tmpData[2])
        data.z = float(tmpData[3])
        data.Q = int(0)
        data.ns = int(tmpData[13])
        data.sdx = float(tmpData[7])
        data.sdy = float(tmpData[8])
        data.sdz = float(tmpData[9])
        data.e = float(0)
        data.n = float(0)
        data.u = float(0)
        #data.dataline = line

        if "Float" == tmpData[16]:
            data.Q = 2
        if "Fixed" == tmpData[16]:
            data.Q = 1

        if intGsod not in posData[intGgpw]:
            posData[intGgpw][intGsod] = data
            epoNum = epoNum + 1
        #print(line)
    logging.info(f"Read flt file : {posFilePath}, epo number is {epoNum}")
    return posData

# if __name__ == '__main__':
#     TleBlocks = {}
#     floPosData = readGreatFltFile("/home/jdhuang/data/projects/run_ppplsq_bdb/PosFile/G/JFNG_2023084_G_2_IF.flt")
