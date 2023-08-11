# neuron_stats
#-|-------------------------------------------------------------------------------------------------------------------------------|
# |  0 |      1      |       2       |     3    |     4    |      5     |   6   |        7       |         8        |      9      |
#-|----|-------------|---------------|----------|----------|------------|-------|----------------|------------------|-------------|
# |time|basket_active|basket_inactive|basket_ref|olm_active|olm_inactive|olm_ref|pyramidal_active|pyramidal_inactive|pyramidal_ref|
#-|-------------------------------------------------------------------------------------------------------------------------------|

from numpy import genfromtxt
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

t = []

guess_freq = 1
guess_amplitude = 0.1
guess_phase = 0
guess_offset = 0
guess_fs = 1

def my_sin(x, freq, amp, phase, fs):
    return amp * np.sin(((2 * np.pi * freq * x) + phase) / fs)

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
    CA3 = []

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "ex"):
                E = data[:, i]
            elif (name == "in"):
                I = data[:, i]
            else:
                CA3 = data[:, i]

    # Plot active neuron stats
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    ax1.title.set_text('Cellular automata simulation of E-I model')
    if (len(E) > 0):
        ax1.plot(t, E)
        ax1.set_ylabel('E_CA1(t)')

    if (len(I) > 0):
        ax2.plot(t, I)
        ax2.set_ylabel('I_CA1I(t)')

    if (len(CA3) > 0):
        ax3.plot(t, I)
        ax3.set_ylabel('Ext(t)')

    plt.tight_layout()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
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

with open('results/ext.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    index = 0
    for row in reader: 
        if (len(row)-1 == 0):
            x = [row[0]] * len(row)
            row = 0
            ax2.scatter(x, row)
        else:
            x = [row[0]] * (len(row)-1)
            row = [int(i) for i in row]
            ax2.scatter(x, row[1:])

ax1.set_ylabel('Excitatory population: neuron number')
ax1.set_ylim(1, 10)
ax2.set_ylabel('External pseudo population: neuron number')
ax2.set_xlabel('time (ms)')
ax2.set_ylim(1, 10)
plt.xticks(np.arange(0, 1000, step=100))
plt.yticks(np.arange(1, 10, step=1))
plt.show()
