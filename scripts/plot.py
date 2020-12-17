from numpy import genfromtxt
import csv
import numpy as np
import matplotlib.pyplot as plt

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
    E_CA1 = []
    I_CA1P = []
    I_CA1I = []
    I_S = []
    CA3 = []
    PS = []

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "pyramidal"):
                E_CA1 = data[:, i]
            elif (name == "hippocamposeptal"):
                I_CA1P = data[:, i]
            elif (name == "interneurons"):
                I_CA1I = data[:, i]
            elif (name == "ca3"):
                CA3 = data[:, i]
            elif (name == "ps"):
                PS = data[:, i]
            else:
                I_S = data[:, i]

    # Plot active neuron stats
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, sharex = True, figsize = (8, 6))
    ax1.set_title('CA1 dynamics')
    ax1.plot(t, E_CA1)
    ax1.set_ylabel('Pyramidal')
    ax2.plot(t, I_CA1P)
    ax2.set_ylabel('Hippocampo-septal')
    ax3.plot(t, I_CA1I)
    ax3.set_ylabel('Interneurons')
    ax4.plot(t, I_S)
    ax4.set_ylabel('Septum')
    ax5.plot(t, CA3)
    ax5.set_ylabel('CA3')
    ax6.plot(t, PS)
    ax6.set_ylabel('PS')
    ax6.set_xlabel('time (ms)')
    plt.show()

with open('results/ca_bin_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')

    bins = []
    E = []
    I_CA1P = []
    I_CA1I = []
    I_S = []
    CA3 = []
    PS = []

    for row in reader:
        if (row[0] == "pyramidal"):
            E = row[1:-1]
            E = [int(i) for i in E]
        elif (row[0] == "septum"):
            I_S = row[1:-1]
            I_S = [int(i) for i in I_S]
        elif (row[0] == "hippocamposeptal"):
            I_CA1P = row[1:-1]
            I_CA1P = [int(i) for i in I_CA1P]
        elif (row[0] == "interneurons"):
            I_CA1I = row[1:-1]
            I_CA1I = [int(i) for i in I_CA1I]
        elif (row[0] == "ca3"):
            CA3 = row[1:-1]
            CA3 = [int(i) for i in CA3]
        elif (row[0] == "ps"):
            PS = row[1:-1]
            PS = [int(i) for i in PS]
        elif (row[0] == "bins"):
            bins = row[1:-1]
            bins = [int(i) for i in bins]

    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, sharex = True, figsize = (9, 9))
    ax1.set_title('CA1 dynamics')
    ax1.bar(bins, E, width = 5)
    ax1.set_ylabel('Pyramidal')

    ax2.bar(bins, I_CA1I, width = 5)
    ax2.set_ylabel('Interneurons')

    ax3.bar(bins, I_CA1P, width = 5)
    ax3.set_ylabel('Hippocampo-septal')

    ax4.bar(bins, I_S, width = 5)
    ax4.set_ylabel('Septum')

    ax5.bar(bins, CA3, width = 5)
    ax5.set_ylabel('CA3')

    ax6.bar(bins, PS, width = 5)
    ax6.set_ylabel('PS')
    ax6.set_xlabel('time (ms)')

    plt.tight_layout()
    plt.show()

    # Plot inactive neuron stats
    # fig = plt.figure(figsize=(8, 6))
    # for i in range(len(row)):
    #     if (row[i].find("_inactive") > 0):
    #         name = row[i][:row[i].find('_')]
    #         plt.plot(data[:, 0], data[:, i], label=name+' cells')
    # plt.title('Neuron stats')
    # plt.xlabel('time')
    # plt.ylabel('inactive neuron %')
    # plt.legend()
    # plt.show()

# fig = plt.figure(2, figsize=(8, 6))
# plt.plot(data[:, 0], data[:, 2], label='active inhibitory neurons')
# plt.plot(data[:, 0], data[:, 8], label='inactive inhibitory neurons')
# plt.plot(data[:, 0], data[:, 5], label='refractory inhibitory neurons')
# plt.title('Neuron stats')
# plt.xlabel('time')
# plt.ylabel('neuron %')
# plt.legend()
# plt.show()

# fig = plt.figure(1, figsize=(8, 6))
# ax1 = fig.add_subplot(211)
# ax1.plot(data[:, 0], data[:, 1], label='active excitatory neurons')
# ax1.plot(data[:, 0], data[:, 7], label='inactive excitatory neurons')
# ax1.plot(data[:, 0], data[:, 4], label='refractory excitatory neurons')
# ax1.set_title('Neuron stats')
# ax1.set_xlabel('time')
# ax1.set_ylabel('neuron %')
# plt.grid(True)
# plt.legend()

# ax2 = fig.add_subplot(212)
# ax2.plot(data[:, 0], data[:, 2], label='active inhibitory neurons')
# ax2.plot(data[:, 0], data[:, 8], label='inactive inhibitory neurons')
# ax2.plot(data[:, 0], data[:, 5], label='refractory inhibitory neurons')
# ax2.set_title('Neuron stats')
# ax2.set_xlabel('time')
# ax2.set_ylabel('neuron %')
# plt.grid(True)
# plt.legend()

# plt.tight_layout()
# plt.show()
