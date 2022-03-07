# -*- coding: utf-8 -*-
"""
Created on Mon May 17 15:02:50 2021

@author: mc16535
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir("C:\\Users\\mc16535\\OneDrive - University of Bristol\\Documents\\Postgrad\\Coding\\Abaqus\\comparison with matlab")

with open("ray tracing test.csv", 'r') as file:
    csv_reader = csv.reader(file, delimiter=',')
    count = 0
    
    time = []
    real = []
    imag = []
    
    for row in file:
        row = row.split(',')
        time.append(row[0])
        comp = row[1].split(' ')
        real.append(float(comp[0]))
        imag.append(float('{}{}'.format(comp[1], comp[2][:-2])))
        
time = np.array(time,dtype=np.float)
real = np.array(real,dtype=np.float)
imag = np.array(imag,dtype=np.float)
data = real + imag*1j

start_idx = 400
end_idx = 700

plt.figure(figsize=(8,6))
plt.plot(time[start_idx:end_idx]*10**6, np.abs(data[start_idx:end_idx]), color='r', linestyle='--', label='abs')
plt.plot(time[start_idx:end_idx]*10**6, np.imag(data[start_idx:end_idx]), color='C1', label='imag')
plt.plot(time[start_idx:end_idx]*10**6, np.real(data[start_idx:end_idx]), color='C0', label='real')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,1,0]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.xlabel('Time (μs)')
plt.ylabel('Signal')

peak_1 = 450
peak_2 = 615

plt.plot([time[peak_1]*10**6, time[peak_1]*10**6], [-np.max(np.abs(data)), np.max(np.abs(data))], color='gray', linestyle='--')
plt.plot([time[peak_2]*10**6, time[peak_2]*10**6], [-np.max(np.abs(data)), np.max(np.abs(data))], color='gray', linestyle='--')

plt.text(time[peak_1-50]*10**6, 0.6, 'd = {:3.1f}mm'.format(1000 * 6320 * time[peak_1] / 2), color='gray')
plt.text(time[peak_2-50]*10**6, 0.6, 'd = {:3.1f}mm'.format(1000 * 6320 * time[peak_2] / 2), color='gray')