from numpy import genfromtxt
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize

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
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))
    ax1.plot(data[:, 0], E)
    ax1.set_xlabel('time, t(ms)')
    ax1.set_ylabel('E(t)')
    ax1.set_title('Excitatory population')

    ax2.plot(data[:, 0], I)
    ax2.set_xlabel('time, t(ms)')
    ax2.set_ylabel('I(t)')
    ax2.set_title('Inhibitory population')

    ax3.plot(data[:, 0], EXT)
    ax3.set_xlabel('time, t(ms)')
    ax3.set_ylabel('External input')
    ax3.set_title('External input')

    plt.tight_layout()
    plt.show()

'''
    width = 10

    bins = []
    for i in range(len(t) + 1):
        if ( i != 0 and i % width == 0):
            bins.append(i)
    values = []
    count = 0
    for i in range(len(E)):
        if (E[i] > 0):
            count = count + 1
        if (i != 0 and i % width == 0):
            values.append(count)
            count = 0
    values.append(count)
'''

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
            E = row[1:]
        elif (row[0] == "in"):
            I = row[1:]
        elif (row[0] == "ext"):
            EXT = row[1:]
        elif (row[0] == "bins"):
            bins = row[1:]

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))
    ax1.bar(bins, E)
    ax1.set_xlabel('time, t(ms)')
    ax1.set_ylabel('E(t)')
    ax1.set_title('Excitatory population')

    ax2.bar(bins, I)
    ax2.set_xlabel('time, t(ms)')
    ax2.set_ylabel('I(t)')
    ax2.set_title('Inhibitory population')

    ax3.bar(bins, EXT)
    ax3.set_xlabel('time, t(ms)')
    ax3.set_ylabel('External input')
    ax3.set_title('External input')

    plt.tight_layout()
    plt.show()
