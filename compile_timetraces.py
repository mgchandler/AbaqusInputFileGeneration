# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:47:57 2021

@author: mc16535
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

num_els = 32

Jobname = sys.argv[1]
tx = int(Jobname.split('_')[-1])
filename = '{}.dat'.format(sys.argv[1])
All_Tts = np.array([])

with open(filename, 'r') as f:
    count = 0
    for line in f:
        
        if line[:4] == "Node":
            if count != 0:
                if len(All_Tts) == 0:
                    All_Tts = np.reshape(Tt, (Tt.shape[0], Tt.shape[1], 1))
                else:
                    All_Tts = np.append(All_Tts, np.reshape(Tt, (Tt.shape[0], Tt.shape[1], 1)), axis=2)
            Tt = np.zeros((2, 1))
        elif line != '\n':
            data = np.reshape(np.double(line.split(' ')), (2, 1))
            Tt = np.append(Tt, data, axis=1)
        else:
            continue
        
        count += 1
        
All_Tts = np.append(All_Tts, np.reshape(Tt, (Tt.shape[0], Tt.shape[1], 1)), axis=2)

nodes_per_element = int(All_Tts.shape[2] / num_els)

Tts = np.zeros((All_Tts.shape[0], All_Tts.shape[1], num_els))
for rx in range(num_els):
    Tts[0, :, rx] = All_Tts[0, :, rx]
    for node in range(nodes_per_element):
        Tts[1, :, rx] += All_Tts[1, :, node + nodes_per_element*rx]
        
for rx in range(num_els):
    with open('tx{}-rx{}.dat'.format(tx, rx+1), 'w') as file:
        for line in range(Tts.shape[1]):
            file.write('{} {}\n'.format(Tts[0, line, rx], Tts[1, line, rx]))
        
# for el in range(num_els):
#     plt.plot(Tts[0, :, el], Tts[1, :, el])
# plt.show()

# plt.plot(All_Tts[0, :, 0], All_Tts[1, :, 0])
# plt.show()

