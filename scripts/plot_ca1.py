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
    I_B = []
    I_BS = []
    I_BP = []

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "pyramidal"):
                E_CA1 = data[:, i]
            elif (name == "hippocamposeptal"):
                I_CA1P = data[:, i]
            elif (name == "interneurons"):
                I_CA1I = data[:, i]
            elif (name == "basket"):
                I_B = data[:, i]
            elif (name == "bistratified"):
                I_BS = data[:, i]
            elif (name == "backprojection"):
                I_BP = data[:, i]
            else:
                I_S = data[:, i]

    '''
    E_CA1_ref = []
    I_CA1P_ref = []
    I_CA1I_ref = []
    I_S_ref = []
    I_B_ref = []
    I_BS_ref = []
    I_BP_ref = []

    for i in range(len(row)):
        if (row[i].find("_ref") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "pyramidal"):
                E_CA1_ref = data[:, i]
            elif (name == "hippocamposeptal"):
                I_CA1P_ref = data[:, i]
            elif (name == "interneurons"):
                I_CA1I_ref = data[:, i]
            elif (name == "basket"):
                I_B_ref = data[:, i]
            elif (name == "bistratified"):
                I_BS_ref = data[:, i]
            else:
                I_S_ref = data[:, i]

    E_CA1_inac = []
    I_CA1P_inac = []
    I_CA1I_inac = []
    I_S_inac = []
    I_B_inac = []
    I_BS_inac = []
    I_BP_inac = []

    for i in range(len(row)):
        if (row[i].find("_inactive") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "pyramidal"):
                E_CA1_inac = data[:, i]
            elif (name == "hippocamposeptal"):
                I_CA1P_inac = data[:, i]
            elif (name == "interneurons"):
                I_CA1I_inac = data[:, i]
            elif (name == "basket"):
                I_B_inac = data[:, i]
            elif (name == "bistratified"):
                I_BS_inac = data[:, i]
            else:
                I_S_inac = data[:, i]
    '''

    # Plot active neuron stats
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1, sharex = True, figsize = (9, 9))
    ax1.set_title('CA1 dynamics')
    ax1.plot(t, E_CA1)
    ax1.set_ylabel('Pyramidal')
    ax2.plot(t, I_B)
    ax2.set_ylabel('Basket')
    ax3.plot(t, I_BS)
    ax3.set_ylabel('Bistratified')
    ax4.plot(t, I_CA1I)
    ax4.set_ylabel('Interneurons')
    ax5.plot(t, I_CA1P)
    ax5.set_ylabel('Hippocampo-septal')
    ax6.plot(t, I_BP)
    ax6.set_ylabel('Backprojection')
    ax7.plot(t, I_S)
    ax7.set_ylabel('Septum')
    ax7.set_xlabel('time (ms)')
    plt.show()

    '''
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, sharex = True, figsize = (9, 9))
    ax1.set_title('CA1 dynamics')
    ax1.plot(t, E_CA1, color='red')
    ax1.plot(t, E_CA1_ref, color='blue')
    ax1.plot(t, E_CA1_inac, color='green')
    ax1.set_ylabel('Pyramidal')
    ax2.plot(t, I_B, color='red')
    ax2.plot(t, I_B_ref, color='blue')
    ax2.plot(t, I_B_inac, color='green')
    ax2.set_ylabel('Basket')
    ax3.plot(t, I_BS, color='red')
    ax3.plot(t, I_BS_ref, color='blue')
    ax3.plot(t, I_BS_inac, color='green')
    ax3.set_ylabel('Bistratified')
    ax4.plot(t, I_CA1I, color='red')
    ax4.plot(t, I_CA1I_ref, color='blue')
    ax4.plot(t, I_CA1I_inac, color='green')
    ax4.set_ylabel('Interneurons')
    ax5.plot(t, I_CA1P, color='red')
    ax5.plot(t, I_CA1P_ref, color='blue')
    ax5.plot(t, I_CA1P_inac, color='green')
    ax5.set_ylabel('hippocampo-septal')
    ax6.plot(t, I_S, color='red')
    ax6.plot(t, I_S_ref, color='blue')
    ax6.plot(t, I_S_inac, color='green')
    ax6.set_ylabel('Septum')
    ax6.set_xlabel('time (ms)')
    plt.show()
    '''

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
