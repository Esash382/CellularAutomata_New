from numpy import genfromtxt
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

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
    CA3 = []
    PS = []
    EC = []

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
            elif (name == "ca3"):
                CA3 = data[:, i]
            elif (name == "ec"):
                EC = data[:, i]
            elif (name == "ps"):
                PS = data[:, i]
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
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9) = plt.subplots(9, 1, sharex = True, figsize = (9, 9))
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
    ax6.plot(t, I_S)
    ax6.set_ylabel('Septum')
    ax7.plot(t, CA3)
    ax7.set_ylabel('CA3')
    ax8.plot(t, EC)
    ax8.set_ylabel('EC')
    ax9.plot(t, PS)
    ax9.set_ylabel('Ext Septum')
    ax9.set_xlabel('time (ms)')
    plt.show()

    # FFT
    # Number of sample points
    N = len(t)
    # sample spacing
    T = 1.0 / len(t)
    yf = fft(E_CA1)
    xf = fftfreq(N, T)[:N//2]
    #plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
    plt.semilogy(xf[1:N//2], 2.0/N * np.abs(yf[1:N//2]), '-b')
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.grid()
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

with open('results/ca_bin_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')

    dataR = genfromtxt('results/ca_bin_stats.csv', delimiter='\t')
    data = dataR.T

    data = np.delete(data, 0, axis=0)
    data = np.delete(data, (len(data)-1), axis=0)

    bins = []
    E = []
    I_CA1P = []
    I_CA1I = []
    I_S = []
    I_B = []
    I_BS = []
    CA3 = []
    PS = []
    EC = []

    for row in reader:
        if (row[0] == "pyramidal"):
            E = row[1:-1]
            E = [int(i) for i in E]
        elif (row[0] == "basket"):
            I_B = row[1:-1]
            I_B = [int(i) for i in I_B]
        elif (row[0] == "bistratified"):
            I_BS = row[1:-1]
            I_BS = [int(i) for i in I_BS]
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
        elif (row[0] == "ec"):
            EC = row[1:-1]
            EC = [int(i) for i in EC]
        elif (row[0] == "bins"):
            bins = row[1:-1]
            bins = [int(i) for i in bins]

    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, sharex = True, figsize = (9, 9))
    ax1.set_title('CA1 dynamics')
    ax1.bar(bins, E)
    ax1.plot(bins, E)
    ax1.set_ylabel('Pyramidal')

    ax2.bar(bins, I_B)
    ax2.plot(bins, I_B)
    ax2.set_ylabel('Basket')

    ax3.bar(bins, I_BS)
    ax3.plot(bins, I_BS)
    ax3.set_ylabel('Bistratified')

    ax4.bar(bins, I_CA1P)
    ax4.plot(bins, I_CA1P)
    ax4.set_ylabel('Hippocampo-septal')

    ax5.bar(bins, I_CA1I)
    ax5.plot(bins, I_CA1I)
    ax5.set_ylabel('Interneurons')

    ax6.bar(bins, I_S)
    ax6.plot(bins, I_S)
    ax6.set_ylabel('Septum')
    ax6.set_xlabel('time (ms)')

    plt.tight_layout()
    plt.show()

    # FFT
    xt = np.linspace(0.0, 1.0/(2.0*len(bins)), int(len(bins)/2))
    yt = scipy.fftpack.fft(E)
    plt.figure(figsize=(8, 5))
    plt.semilogx(xt[1:], 2.0/len(bins) * np.abs(yt[0:int(len(bins)/2)])[1:])
    plt.title('FFT plot: Pyramidal cell population')
    plt.xlabel('time')
    plt.ylabel('frequency')
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
