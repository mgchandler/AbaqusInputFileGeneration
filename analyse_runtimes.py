# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:13:23 2021

@author: mc16535
"""

import os
import time

def parseTimeString(timeString):
    """ Converts time as printed by BP to struct_time object from time module"""
    timeString = timeString.replace('BST', 'GMT Daylight Time')
    return time.strptime(timeString, "%a %b %d %H:%M:%S %Z %Y")

os.chdir("C:\\Users\\mc16535\\OneDrive - University of Bristol\\Documents\\Postgrad\\Coding\\Abaqus\\FMC Generation\\v7\\Output\\Sensitivity Variation\\io")

fileList = os.listdir()

times = []

for filename in fileList:
    if ".sh.o" in filename and "convert_to_mat" not in filename:
        startTimes = []
        endTimes = []
        
        data = {}
        data['filename'] = filename
        
        with open(filename, 'r') as file:
            for line in file:
                if "Start Time: " in line:
                    startTimes.append(parseTimeString(line.split("Start Time: ")[1].replace('\n', '')))
                if "End Time: " in line:
                    endTimes.append(parseTimeString(line.split("End Time: ")[1].replace('\n', '')))
                    
            data['Start Times'] = startTimes
            data['End Times'] = endTimes
            
        times.append(data)
        
runtimes = []

avgRuntime = 0
avgCounter = 0
        
for file in times:
    start = time.mktime(file['Start Times'][0])
    end = time.mktime(file['End Times'][-1])
    runtime = end - start
    runtimes.append([file['filename'], runtime])
    if "RunArrayJob" in file['filename']:
        avgRuntime += runtime
        avgCounter += 1
avgRuntime /= avgCounter