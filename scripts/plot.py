

from numpy import genfromtxt
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from PIL import Image
from io import BytesIO

# neuron_stats
#-|-------------------------------------------------------------------------------------------------------------------------------|
# |  0 |      1      |       2       |     3    |     4    |      5     |   6   |        7       |         8        |      9      |
#-|----|-------------|---------------|----------|----------|------------|-------|----------------|------------------|-------------|
# |time|basket_active|basket_inactive|basket_ref|olm_active|olm_inactive|olm_ref|pyramidal_active|pyramidal_inactive|pyramidal_ref|
#-|-------------------------------------------------------------------------------------------------------------------------------|

time_start = 0
time_end = 0

with open('results/ca_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    row = next(reader)

    data = genfromtxt('results/ca_stats.csv', delimiter='\t')
    data = np.delete(data, (0), axis=0)

    t = data[:, 0]
    time_start = t[0]
    time_end = t[-1]
    E = []
    I = []
    EXT = []

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "ex"):
                E = data[:, i]
            elif (name == "in"):
                I = data[:, i]
            else:
                EXT = data[:, i]

    # Plot active neuron stats
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))
    ax1.plot(t, E)
    ax1.set_xlabel('time, t(ms)')
    ax1.set_ylabel('E(t)')
    ax1.set_title('Excitatory population')
    ax1.grid()

    ax2.plot(t, I)
    ax2.set_xlabel('time, t(ms)')
    ax2.set_ylabel('I(t)')
    ax2.set_title('Inhibitory population')
    ax2.grid()

    ax3.plot(t, EXT)
    ax3.set_xlabel('time, t(ms)')
    ax3.set_ylabel('External input')
    ax3.set_title('External input')
    ax3.grid()

    plt.tight_layout()

'''
    plt.figure(figsize=(8, 3))
    dataR = genfromtxt('results/ex.csv', delimiter='\t')
    dataPRand = genfromtxt('results/ca_p_rand_stats.csv', delimiter='\t')
    dataT = dataR.T
    dataS = np.delete(dataT, 0, axis=0)
    dataS[ dataS ==-1 ] = np.nan
    dataS1 = dataS.copy()
    dataS2 = dataS.copy()
    dataS3 = dataS.copy()
    dataS4 = dataS.copy()
    dataS5 = dataS.copy()

    for i in range(len(dataS)):
        for j in range(len(dataS[i])):
            for p in range(len(dataPRand)):
                if (np.isnan(dataS[i][j])):
                    continue
                else:
                    if (dataS[i][j] not in dataPRand[p]):
                        if p == 0:
                            dataS1[i][j] = np.nan
                        elif p == 1:
                            dataS2[i][j] = np.nan
                        elif p == 2:
                            dataS3[i][j] = np.nan
                        elif p == 3:
                            dataS4[i][j] = np.nan
                        elif p == 4:
                            dataS5[i][j] = np.nan

        plt.scatter(t, dataS5[i], marker="*", s=20, color='red')
        plt.scatter(t, dataS1[i], marker="s", s=50, color='blue')
        plt.scatter(t, dataS2[i], marker="^", s=40, color='yellow')
        plt.scatter(t, dataS3[i], marker="o", s=30, color='green')
        plt.scatter(t, dataS4[i], marker="_", s=30, color='pink')

    plt.title('Spike raster plot')
    plt.xlabel('time (ms)')
    plt.ylabel('Excitatory population: neuron number')
    # plt.legend((dots, stars), ('Truely recalled neurons', 'Falsely recalled neurons'))
'''

# FFT
# Number of sample points
N = len(t)
# sample spacing
T = 1.0 / len(t)
yf = fft(E)
xf = fftfreq(N, T)[:N//2]
plt.figure()
plt.plot(xf[0:100], 1.0 / 10 * np.abs(yf[0:N//10]))
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.grid()
plt.show()
