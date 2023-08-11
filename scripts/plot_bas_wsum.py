#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 23:07:36 2023

@author: ashraya
"""

import matplotlib.pyplot as plt
import csv

filename = open('../out.txt', 'r')
file = csv.reader(filename)

data = []
for col in file:
    data.append(int(col[0]))

plt.figure()
plt.plot(data[0:int(len(data))])
x = [8]
x = x * int(len(data))
plt.plot(x)
plt.ylim(0, 30)
plt.show()