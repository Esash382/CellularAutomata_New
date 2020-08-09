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
    E_CA3 = []
    I_CA3P = []
    I_CA3I = []
    I_B3 = []
    I_BS3 = []

    I_S = []

    E_CA1 = []
    I_CA1P = []
    I_CA1I = []
    I_B1 = []
    I_BS1 = []
    I_BP = []

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "pyramidal1"):
                E_CA1 = data[:, i]
            elif (name == "hippocamposeptal1"):
                I_CA1P = data[:, i]
            elif (name == "interneurons1"):
                I_CA1I = data[:, i]
            elif (name == "basket1"):
                I_B1 = data[:, i]
            elif (name == "bistratified1"):
                I_BS1 = data[:, i]
            elif (name == "backprojection"):
                I_BP = data[:, i]
            if (name == "pyramidal3"):
                E_CA3 = data[:, i]
            elif (name == "hippocamposeptal3"):
                I_CA3P = data[:, i]
            elif (name == "interneurons3"):
                I_CA3I = data[:, i]
            elif (name == "basket3"):
                I_B3 = data[:, i]
            elif (name == "bistratified3"):
                I_BS3 = data[:, i]
            else:
                I_S = data[:, i]

    # Plot active neuron stats
    fig1, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1, sharex = True, figsize = (9, 9))
    ax1.set_title('CA1 dynamics')
    ax1.plot(t, E_CA1)
    ax1.set_ylabel('Pyramidal')
    ax2.plot(t, I_B1)
    ax2.set_ylabel('Basket')
    ax3.plot(t, I_BS1)
    ax3.set_ylabel('Bistratified')
    ax4.plot(t, I_CA1I)
    ax4.set_ylabel('Interneurons')
    ax5.plot(t, I_CA1P)
    ax5.set_ylabel('hippocampo-septal')
    ax6.plot(t, I_BP)
    ax6.set_ylabel('Backprojection')
    ax7.plot(t, I_S)
    ax7.set_ylabel('Septum')
    ax7.set_xlabel('time (ms)')

    fig2, (ax11, ax22, ax33, ax44, ax55, ax66) = plt.subplots(6, 1, sharex = True, figsize = (9, 9))
    ax11.set_title('CA3 dynamics')
    ax11.plot(t, E_CA3)
    ax11.set_ylabel('Pyramidal')
    ax22.plot(t, I_B3)
    ax22.set_ylabel('Basket')
    ax33.plot(t, I_BS3)
    ax33.set_ylabel('Bistratified')
    ax44.plot(t, I_CA3I)
    ax44.set_ylabel('Interneurons')
    ax55.plot(t, I_CA3P)
    ax55.set_ylabel('hippocampo-septal')
    ax66.plot(t, I_S)
    ax66.set_ylabel('Septum')
    ax66.set_xlabel('time (ms)')
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
