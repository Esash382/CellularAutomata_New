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

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "ex"):
                E = data[:, i]
            else:
                I = data[:, i]

    # Plot active neuron stats
    fig = plt.figure(1, figsize=(8, 6)) 
    plt.plot(t, E)
    plt.xlabel('time, t(ms)')
    plt.ylabel('E(t)')
    plt.title('Excitatory population')
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

def func1(x, a, b, c):
    # a Gaussian distribution
    # return a * np.exp(-(x-b)**2/(2*c**2))
    return a + b*x + c*x*x

def func(x, a, b):
    return a * np.sin(b * x)

with open('results/ca_bin_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    row = next(reader)

    dataR = genfromtxt('results/ca_bin_stats.csv', delimiter='\t')
    data = dataR.T
    data = np.delete(data, (len(data)-1), axis=0)

    bins = data[:, 0]
    E = []

    for i in range(len(data[0])):
        if (i%2 != 0):
            E = data[:, i]
            break

    width = 1/1000
    # fit bar plot data using curve_fit
    # popt, pcov = curve_fit(func, bins, E)
    popt, pcov = optimize.curve_fit(func, bins, E, p0=[2, 2])

    x = np.linspace(bins[0], bins[-1], 100)
    y = func(x, *popt)

    plt.figure()
    plt.bar(bins, E, color='b', width=5)
    plt.xticks(bins)
    #plt.plot(bins, E, c='g')
    #plt.plot(bins, y + max(E), c='g')

    plt.show()
