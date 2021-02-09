from numpy import genfromtxt
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# neuron_stats
#-|-------------------------------------------------------------------------------------------------------------------------------|
# |  0 |      1      |       2       |     3    |     4    |      5     |   6   |        7       |         8        |      9      |
#-|----|-------------|---------------|----------|----------|------------|-------|----------------|------------------|-------------|
# |time|basket_active|basket_inactive|basket_ref|olm_active|olm_inactive|olm_ref|pyramidal_active|pyramidal_inactive|pyramidal_ref|
#-|-------------------------------------------------------------------------------------------------------------------------------|

with open('results/ca_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    row = next(reader)

    data = genfromtxt('results/ca_stats.csv', delimiter='\t')
    data = np.delete(data, (0), axis=0)

    t = data[:, 0]
    E = []
    I = []
    EXT = []

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "ex"):
                E = data[:, i]
            elif (name == "ext"):
                EXT = data[:, i]
            else:
                I = data[:, i]

    # Plot active neuron stats
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    ax1.plot(data[:, 0], E)
    ax1.set_xlabel('time, t(ms)')
    ax1.set_ylabel('E(t)')
    ax1.set_title('Excitatory population')

    ax2.plot(data[:, 0], I)
    ax2.set_xlabel('time, t(ms)')
    ax2.set_ylabel('I(t)')
    ax2.set_title('Inhibitory population')

    plt.tight_layout()

'''
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
with open('results/ex.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader: 
        if (len(row)-1 == 0):
            x = [row[0]] * len(row)
            row = 0
            ax1.scatter(x, row)
        else:
            x = [row[0]] * (len(row)-1)
            row = [int(i) for i in row]
            ax1.scatter(x, row[1:])

with open('results/in.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader: 
        if (len(row)-1 == 0):
            x = [row[0]] * len(row)
            row = 0
            ax2.scatter(x, row)
        else:
            x = [row[0]] * (len(row)-1)
            row = [int(i) for i in row]
            ax2.scatter(x, row[1:])

with open('results/ext.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    index = 0
    for row in reader: 
        if (len(row)-1 == 0):
            x = [row[0]] * len(row)
            row = 0
            ax3.scatter(x, row)
        else:
            x = [row[0]] * (len(row)-1)
            row = [int(i) for i in row]
            ax3.scatter(x, row[1:])

ax1.set_ylabel('Excitatory population: neuron number')
ax2.set_ylabel('Inhibitory population: neuron number')
ax3.set_ylabel('External pseudo population: neuron number')
ax3.set_xlabel('time (ms)')
'''
plt.show()

'''
    # FFT
    # Number of sample points
    N = len(t)
    # sample spacing
    T = 1.0 / len(t)
    yf = fft(E)
    xf = fftfreq(N, T)[:N//2]
    plt.figure()
    plt.plot(xf, 1.0 / 10 * np.abs(yf[0:N//2]))
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.grid()
    plt.show()

with open('results/ca_bin_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')

    dataR = genfromtxt('results/ca_bin_stats.csv', delimiter='\t')
    data = dataR.T

    data = np.delete(data, 0, axis=0)
    data = np.delete(data, (len(data)-1), axis=0)

    bins = []
    E = []
    I = []
    EXT = []

    for row in reader:
        if (row[0] == "ex"):
            E = row[1:-1]
            E = [int(i) for i in E]
        elif (row[0] == "in"):
            I = row[1:-1]
            I = [int(i) for i in I]
        elif (row[0] == "ext"):
            EXT = row[1:-1]
            EXT = [int(i) for i in EXT]
        elif (row[0] == "bins"):
            bins = row[1:-1]
            bins = [int(i) for i in bins]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    ax1.bar(bins, E)
    ax1.plot(bins, E)
    ax1.set_xlabel('time, t(ms)')
    ax1.set_ylabel('E(t)')
    ax1.set_title('Excitatory population')

    ax2.bar(bins, I)
    ax2.plot(bins, I)
    ax2.set_xlabel('time, t(ms)')
    ax2.set_ylabel('I(t)')
    ax2.set_title('Inhibitory population')

    plt.tight_layout()
    plt.show()
'''
