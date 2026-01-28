# -*- coding: utf-8 -*-
# ===================================================================
# 读取GREAT的flt格式, 并绘图
# ===================================================================

import logging
import os
import sys
from matplotlib import pyplot as plt

pyDirName, pyFileName = os.path.split(os.path.abspath(__file__))
sys.path.append(f"{pyDirName}/../")
from PythonScripts.dependency.gnss_greatPos_plot import getPosTimeIndex
from PythonScripts.dependency.gnss_cord_tool import ecef2enu, ecef2pos
from PythonScripts.dependency.gnss_math_tool import get_rms, get_standard_deviation


def gnss_great_draw_compare_flt(site_name, type1, type2, fltData1, fltData2, crdData, savePath):
    """
    绘制两个flt文件的比较结果
    """
    # 检查文件夹是否存在
    picDir = os.path.dirname(savePath)
    if not os.path.exists(picDir):
        os.makedirs(picDir)
    if not os.path.exists(savePath):
        os.makedirs(savePath)

    # 检查精密坐标是否存在，没有的话直接退出了
    preciseCrd = [0.0, 0.0, 0.0]
    if site_name not in crdData.site_list:
        return False
    else:
        preciseCrd[0] = float(crdData.site_stax[site_name])
        preciseCrd[1] = float(crdData.site_stay[site_name])
        preciseCrd[2] = float(crdData.site_staz[site_name])

    # 创建保存文件
    picPath = os.path.join(savePath, f"{site_name}_{type1}_vs_{type2}_flt.png")
    sumPath = os.path.join(savePath, f"{site_name}_{type1}_vs_{type2}_flt.sum")
    picTitle = f"{site_name}_{type1}_vs_{type2}"

    # 创建汇总文件
    logging.info(f"Plot {picPath}")
    fwSumFile = open(sumPath, 'w')

    # plt.rcParams['agg.path.chunksize'] = 1000
    font1 = {'weight': 600, 'size': 15}
    font2 = {'weight': 550, 'size': 12}
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax1 = plt.subplot(3,1,1)
    plt.title(site_name, font1)
    ax2 = plt.subplot(3,1,2)
    ax3 = plt.subplot(3,1,3)

    ax1.set_ylabel('East/m', font1)
    ax1.grid(axis="y")
    # ax1.set_ylim(-0.1, 0.1)
    ax1.tick_params(axis='both', colors='black', direction='out', labelsize=11, width=1, length=1, pad=5)

    ax2.set_ylabel('North/m', font1)
    ax2.grid(axis="y")
    # ax2.set_ylim(-0.1, 0.1)
    ax2.tick_params(axis='both', colors='black', direction='out', labelsize=11, width=1, length=1, pad=5)
    
    ax3.set_ylabel('Up/m', font1)
    ax3.grid(axis="y")
    # ax3.set_ylim(-0.1, 0.1)
    ax3.set_xlabel('Period of processing',font1)
    ax3.tick_params(axis='both', colors='black', direction='out', labelsize=11, width=1, length=1, pad=5)


    # 小时格式的坐标轴
    xSod = list()
    xStr = list()
    for i in range(0, 24, 2):
        xSod.append(i * 3600)
        xStr.append(f"{i:02d}")

    ax1.set_xticks(xSod)
    ax1.set_xticklabels(xStr)
    ax2.set_xticks(xSod)
    ax2.set_xticklabels(xStr)
    ax3.set_xticks(xSod)
    ax3.set_xticklabels(xStr)

    # 存数据
    plotX1 = list()
    plotY1 = dict()
    plotY1[0] = list()
    plotY1[1] = list()
    plotY1[2] = list()

    plotX2 = list()
    plotY2 = dict()
    plotY2[0] = list()
    plotY2[1] = list()
    plotY2[2] = list()

    # 根据时间循环
    # 绘图
    ## 文件1
    timeIndex = getPosTimeIndex(fltData1, 1)
    allNum = len(timeIndex)
    fixNum = 0
    floNum = 0
    for i in range(0,allNum):
        strGpswSow = str(timeIndex[i])
        gpsw = int(strGpswSow[0:strGpswSow.find(".")])
        gsow = int(strGpswSow[strGpswSow.find(".")+1:])
        # logging.info(f"Plot {strGpswSow}")
        if gpsw in fltData1 and gsow in fltData1[gpsw]:
            # 计算e, n, u坐标
            posCrd = [fltData1[gpsw][gsow].x, fltData1[gpsw][gsow].y, fltData1[gpsw][gsow].z]
            difCrd = [0, 0, 0]
            pos = ecef2pos(preciseCrd)
            for x in range(0, 3):
                difCrd[x] = posCrd[x] - preciseCrd[x]
            enu = ecef2enu(pos, difCrd)

            if fltData1[gpsw][gsow].Q == 1:
                fixNum = fixNum + 1
            elif fltData1[gpsw][gsow].Q == 2:
                floNum = floNum + 1
            
            sod = gsow % 86400
            plotX1.append(sod)
            plotY1[0].append(enu[0])
            plotY1[1].append(enu[1])
            plotY1[2].append(enu[2])

            # 输出
            writeLine = f"{type1} {gsow} {enu[0]:16.8f} {enu[1]:16.8f} {enu[2]:16.8f}"
            fwSumFile.write(writeLine + "\n")

    ## 文件2
    timeIndex = getPosTimeIndex(fltData2, 1)
    allNum = len(timeIndex)
    fixNum = 0
    floNum = 0

    #ylx
    col=[]
    colors = ['black','r', 'g']


    for i in range(0,allNum):
        strGpswSow = str(timeIndex[i])
        gpsw = int(strGpswSow[0:strGpswSow.find(".")])
        gsow = int(strGpswSow[strGpswSow.find(".")+1:])
        # logging.info(f"Plot {strGpswSow}")
        if gpsw in fltData2 and gsow in fltData2[gpsw]:
            # 计算e, n, u坐标
            posCrd = [fltData2[gpsw][gsow].x, fltData2[gpsw][gsow].y, fltData2[gpsw][gsow].z]
            difCrd = [0, 0, 0]
            pos = ecef2pos(preciseCrd)
            for x in range(0, 3):
                difCrd[x] = posCrd[x] - preciseCrd[x]
            enu = ecef2enu(pos, difCrd)

            if fltData2[gpsw][gsow].Q == 1:
                fixNum = fixNum + 1
            elif fltData2[gpsw][gsow].Q == 2:
                floNum = floNum + 1    
            col.append(colors[fltData2[gpsw][gsow].Q])
            fixingrate=fixNum/(floNum+fixNum)

            sod = gsow % 86400
            plotX2.append(sod)
            plotY2[0].append(enu[0])
            plotY2[1].append(enu[1])
            plotY2[2].append(enu[2])

            # 输出
            writeLine = f"{type2} {gsow} {enu[0]:16.8f} {enu[1]:16.8f} {enu[2]:16.8f}"
            fwSumFile.write(writeLine + "\n")


    

    ax1.plot(plotX1, plotY1[0], c='b', marker='.', linestyle='None', label=f"{type1}",linewidth=1,markersize=1)
    ax1.plot(plotX2, plotY2[0], c='r', marker='.', linestyle='None', label=f"{type2}",linewidth=1,markersize=1)
    #ax1.scatter(plotX2, plotY2[0], c=col, marker='.', linestyle='None', label=f"{type2}")
    ax2.plot(plotX1, plotY1[1], c='b', marker='.', linestyle='None', label=f"{type1}",linewidth=1,markersize=1)
    ax2.plot(plotX2, plotY2[1], c='r', marker='.', linestyle='None', label=f"{type2}",linewidth=1,markersize=1)
    #ax2.scatter(plotX2, plotY2[1], c=col, marker='.', linestyle='None', label=f"{type2}")
    ax3.plot(plotX1, plotY1[2], c='b', marker='.', linestyle='None', label=f"{type1}",linewidth=1,markersize=1)
    ax3.plot(plotX2, plotY2[2], c='r', marker='.', linestyle='None', label=f"{type2}",linewidth=1,markersize=1)
    #ax3.scatter(plotX2, plotY2[2], c=col, marker='.', linestyle='None', label=f"{type2}")
    ax1.legend()
    ax2.legend()
    ax3.legend()

    # 写入统计文件
    ## 文件1
    rms_e1 = get_rms(plotY1[0])
    rms_n1 = get_rms(plotY1[1])
    rms_u1 = get_rms(plotY1[2])
    std_e1 = get_standard_deviation(plotY1[0])
    std_n1 = get_standard_deviation(plotY1[1])
    std_u1 = get_standard_deviation(plotY1[2])
    writeLine = f"{type1} SumMark ENU RMS1 = {rms_e1:16.8f} {rms_n1:16.8f} {rms_u1:16.8f}"
    fwSumFile.write(writeLine + "\n")
    writeLine = f"{type1} SumMark ENU STD1 = {std_e1:16.8f} {std_n1:16.8f} {std_u1:16.8f}"
    fwSumFile.write(writeLine + "\n")
    ax1.text(0.7, 1.20, f"{type1} RMS = {rms_e1:6.3f},{rms_n1:6.3f},{rms_u1:6.3f}", bbox=dict(facecolor='none', alpha=0.0, pad=6), fontdict=font2, transform=ax1.transAxes)
    ax1.text(0.7, 1.07, f"{type1} STD = {std_e1:6.3f},{std_n1:6.3f},{std_u1:6.3f}", bbox=dict(facecolor='none', alpha=0.0, pad=6), fontdict=font2, transform=ax1.transAxes)

    ## 文件2
    rms_e2 = get_rms(plotY2[0])
    rms_n2 = get_rms(plotY2[1])
    rms_u2 = get_rms(plotY2[2])
    std_e2 = get_standard_deviation(plotY2[0])
    std_n2 = get_standard_deviation(plotY2[1])
    std_u2 = get_standard_deviation(plotY2[2])
    writeLine = f"{type2} SumMark ENU RMS2 = {rms_e2:16.8f} {rms_n2:16.8f} {rms_u2:16.8f}"
    fwSumFile.write(writeLine + "\n")
    writeLine = f"{type2} SumMark ENU STD2 = {std_e2:16.8f} {std_n2:16.8f} {std_u2:16.8f}"
    fwSumFile.write(writeLine + "\n")

    ## 两个文件差异
    writeLine = f"{type2} SumMark RMS DIF = {rms_e1-rms_e2:16.8f} {rms_n1-rms_n2:16.8f} {rms_u1-rms_u2:16.8f}"
    fwSumFile.write(writeLine + "\n")
    writeLine = f"{type2} SumMark STD DIF = {std_e1-std_e2:16.8f} {std_n1-std_n2:16.8f} {std_u1-std_u2:16.8f}"
    fwSumFile.write(writeLine + "\n")

    writeLine= f"fixingrate = {fixingrate*100:6.8f}%"
    fwSumFile.write(writeLine + "\n")

    ax1.text(0.0, 1.20, f"{type2} RMS = {rms_e2:6.3f},{rms_n2:6.3f},{rms_u2:6.3f}", bbox=dict(facecolor='none', alpha=0.0, pad=6), fontdict=font2, transform=ax1.transAxes)
    ax1.text(0.0, 1.07, f"{type2} STD = {std_e2:6.3f},{std_n2:6.3f},{std_u2:6.3f}", bbox=dict(facecolor='none', alpha=0.0, pad=6), fontdict=font2, transform=ax1.transAxes)
    fig.text(0.9, 0.04, f"fixing rate = {fixingrate*100:6.2f}%",ha='right', bbox=dict(facecolor='none', alpha=0.0, pad=6), fontdict=font2)

    # 关闭
    fwSumFile.close()
    fig.subplots_adjust()
    fig.tight_layout()

    #plt.grid(axis="y")
    plt.savefig(picPath, dpi=300)
    plt.close()
    plt.clf()
    #plt.show()
    return True
