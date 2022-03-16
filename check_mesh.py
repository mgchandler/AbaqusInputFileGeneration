# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 13:05:47 2022

@author: mc16535
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import time

def secsToString(seconds):
    return "{:02d}:{:02d}:{:02d}".format(int(seconds/3600), int((seconds%3600)/60), int(seconds%60))

def formatTo3SF(itn):
    if np.log10(itn) >= 2:
        return "{:4d}".format(int(itn))
    elif np.log10(itn) >= 1:
        return "{:3.1f}".format(itn)
    else:
        return "{:3.2f}".format(itn)

def itnToString(itn):
    if np.log10(itn) >= 18:
        return "{}E".format(formatTo3SF(itn*10**-18))
    elif np.log10(itn) >= 15:
        return "{}P".format(formatTo3SF(itn*10**-15))
    elif np.log10(itn) >= 12:
        return "{}T".format(formatTo3SF(itn*10**-12))
    elif np.log10(itn) >= 9:
        return "{}G".format(formatTo3SF(itn*10**-9))
    elif np.log10(itn) >= 6:
        return "{}M".format(formatTo3SF(itn*10**-6))
    elif np.log10(itn) >= 3:
        return "{}k".format(formatTo3SF(itn*10**-3))
    else:
        return "{:4d} ".format(int(itn))

def updateProgress(startTime, currentItn, totalItn):
    currentTime = time.time()
    timeElapsed = secsToString(currentTime - startTime)
    timeRemains = secsToString((currentTime - startTime) * (float(totalItn) / float(currentItn+1) - 1))
    
    barChars = 50
    percPrint = int(100*float(currentItn+1)/float(totalItn))
    percBarFull = int(barChars*float(currentItn+1)/float(totalItn))
    bar = "".join(["#" for ii in range(percBarFull)]) + "".join([" " for ii in range(barChars-percBarFull)])
    
    print("\r{:3d}%|{}|{}/{} Elap:{} ETA:{}  ".format(percPrint, bar, itnToString(currentItn), itnToString(totalItn), timeElapsed, timeRemains), end='')
    sys.stdout.flush()

nodes = np.zeros((3, 3250000, 5))
elements = np.zeros((4, 6250000, 5), dtype=int)

for ii in [0,1]:
    if ii == 0:
        filename = "I_45npw_02mm_b.inp"
    else:
        filename = "I_45npw_02mm_{}.inp".format(ii)
        
    in_node = False
    in_ele  = False
    nodecount = 0
    elecount  = 0
    
    t1 = time.time_ns()*10**-9
    print("Starting {}".format(filename))
    with open(filename, 'r') as file:
        for count, line in enumerate(file):
            if "*Node" in line:
                in_node = True
                continue
            elif "*Element" in line:
                in_node = False
                in_ele  = True
                continue
            elif "*Nset" in line:
                break
            
            if in_node:
                nodes[:, nodecount, ii] = [float(val) for val in line.split(',')]
                nodecount += 1
            if in_ele:
                elements[:, elecount, ii] = [float(val) for val in line.split(',')]
                elecount += 1
            
            if (nodecount + elecount+1)%(int((nodes.shape[1] + elements.shape[1])/1000)+1) == 0:
                updateProgress(t1, nodecount + elecount, nodes.shape[1] + elements.shape[1])
    print("")

plotIdx1 = 0
plotIdx2 = 1
midX = 7.118e-3
radX = 1e-3
midY = -26.563e-3
radY = .6e-3
plt.figure(figsize=(10,6), dpi=100)
# for ii in range(elements.shape[1]):
#     if abs(nodes[1, elements[1, ii, plotIdx1]-1, plotIdx1] - midX) < radX and abs(nodes[2, elements[1, ii, plotIdx1]-1, plotIdx1] - midY) < radY:
#         plt.plot(
#             [nodes[1, elements[1, ii, plotIdx1]-1, plotIdx1], nodes[1, elements[2, ii, plotIdx1]-1, plotIdx1], nodes[1, elements[3, ii, plotIdx1]-1, plotIdx1], nodes[1, elements[1, ii, plotIdx1]-1, plotIdx1]],
#             [nodes[2, elements[1, ii, plotIdx1]-1, plotIdx1], nodes[2, elements[2, ii, plotIdx1]-1, plotIdx1], nodes[2, elements[3, ii, plotIdx1]-1, plotIdx1], nodes[2, elements[1, ii, plotIdx1]-1, plotIdx1]],
#             c='C0', linewidth=1.5)
# for ii in range(elements.shape[1]):
#     if abs(nodes[1, elements[1, ii, plotIdx2]-1, plotIdx2] - midX) < radX and abs(nodes[2, elements[1, ii, plotIdx2]-1, plotIdx2] - midY) < radY:
#         plt.plot(
#             [nodes[1, elements[1, ii, plotIdx2]-1, plotIdx2], nodes[1, elements[2, ii, plotIdx2]-1, plotIdx2], nodes[1, elements[3, ii, plotIdx2]-1, plotIdx2], nodes[1, elements[1, ii, plotIdx2]-1, plotIdx2]],
#             [nodes[2, elements[1, ii, plotIdx2]-1, plotIdx2], nodes[2, elements[2, ii, plotIdx2]-1, plotIdx2], nodes[2, elements[3, ii, plotIdx2]-1, plotIdx2], nodes[2, elements[1, ii, plotIdx2]-1, plotIdx2]],
#             c='C1', linewidth=.6)
plt.scatter(nodes[1, :, plotIdx1], nodes[2, :, plotIdx1], s=.5, c='C0')
plt.scatter(nodes[1, :, plotIdx2], nodes[2, :, plotIdx2], s=.5, c='C1')
plt.xlim(midX - radX, midX + radX)
plt.ylim(midY - radY, midY + radY)
plt.show()