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
from scipy.fftpack import fft, fftfreq
import os
import sys

t = []

def my_sin(x, freq, amp, phase):
    return amp * np.sin(((2 * np.pi * freq * x / 1000) + phase))

curPath = os.getcwd()
if ("scripts" in curPath):
    os.chdir("../")

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
    HIPP = []
    S = []
    B = []
    BS = []
    CA3 = []
    PS = []
    EC = []

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "ex"):
                E = data[:, i]
            elif (name == "in"):
                I = data[:, i]
            elif (name == "hipp"):
                HIPP = data[:, i]
            elif (name == "septum"):
                S = data[:, i]
            elif (name == "bas"):
                B = data[:, i]
            elif (name == "bis"):
                BS = data[:, i]
            elif (name == "ec"):
                EC = data[:, i]
            elif (name == "ps"):
                PS = data[:, i]
            else:
                CA3 = data[:, i]

    # Plot active neuron stats
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9) = plt.subplots(9, 1, figsize=(10, 10), sharex=True)
    ax1.title.set_text('Cellular automata simulation of CA1 circuit')
    if (len(E) > 0):
        ax1.plot(t, E)
        ax1.set_ylim(0, 0.3)
        ax1.set_ylabel('E_CA1(t)')
        guess_phase= 0
        guess_amplitude = 0.25
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        e_guess = my_sin(t, *p0)
        ax1.plot(t, e_guess, color='green')

    if (len(B) > 0):
        ax2.plot(t, B)
        ax2.set_ylim(0, 0.3)
        ax2.set_ylabel('I_B(t)')
        guess_phase = 0
        guess_amplitude = 0.2
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        ip_guess = my_sin(t, *p0)
        ax2.plot(t, ip_guess, color='green')

    if (len(BS) > 0):
        ax3.plot(t, BS)
        ax3.set_ylim(0, 0.3)
        ax3.set_ylabel('I_BS(t)')
        guess_phase = 0
        guess_amplitude = 0.25
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        i_guess = my_sin(t, *p0)
        ax3.plot(t, i_guess, color='green')

    if (len(I) > 0):
        ax4.plot(t, I)
        ax4.set_ylim(0, 0.3)
        ax4.set_ylabel('I_CA1I(t)')
        guess_phase = -2.5
        guess_amplitude = 0.25
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax4.plot(t, s_guess, color='green')

    if (len(HIPP) > 0):
        ax5.plot(t, HIPP)
        ax5.set_ylim(0, 0.3)
        ax5.set_ylabel('I_CA1P(t)')
        guess_phase = -0.2
        guess_amplitude = 0.25
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax5.plot(t, s_guess, color='green')
        
    if (len(S) > 0):
        ax6.plot(t, S)
        ax6.set_ylim(0, 0.25)
        ax6.set_ylabel('I_S(t)')
        guess_phase = 1.3
        guess_amplitude = 0.2
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax6.plot(t, s_guess, color='green')

    if (len(EC) > 0):
        ax7.plot(t, EC)
        ax7.set_ylabel('EC')

    if (len(CA3) > 0):
        ax8.plot(t, CA3)
        ax8.set_ylabel('CA3')

    if (len(PS) > 0):
        ax9.plot(t, PS)
        ax9.set_ylabel('PS')
        ax9.set_xlabel('time, t(ms)')

    plt.tight_layout()
#    plt.savefig('ca1_theta_gamma.png', dpi=500)
    plt.show()

'''
    plt.figure(figsize=(8, 3))
    dataR = genfromtxt('results/ex.csv', delimiter='\t')
    dataPRand = genfromtxt('results/ca_p_rand_stats.csv', delimiter='\t')
    dataT = dataR.T
    dataS = np.delete(dataT, 0, axis=0)
    dataS[ dataS ==-1 ] = np.nan
    dataS2 = dataS.copy()

    totalNeurons = len(dataR[0])-1
    totalSubsetNeurons = len(dataPRand[0])-1
    tNeurons = 0
    fNeurons = 0

    dots = None
    stars = None
    for (i, k) in zip(dataS, dataS2):
        for j in range(len(i)):
            if (i[j] != np.nan and i[j] not in dataPRand):
                i[j] = np.nan
            else:
                k[j] = np.nan
        tCount = np.count_nonzero(~np.isnan(i))
        fCount = np.count_nonzero(~np.isnan(k))
        if (tCount > 0):
            tNeurons = tNeurons + 1
        if (fCount > 0):
            fNeurons = fNeurons + 1
        dots = plt.scatter(t, i, marker=".", s=150, color='b')
        stars = plt.scatter(t, k, marker="*", s=20, color='r')

    print("Truely recalled neurons = ", (tNeurons / totalSubsetNeurons) * 100, "%")
    print("Falsely recalled neurons = ", (fNeurons / (totalNeurons - totalSubsetNeurons)) * 100, "%")

    plt.title('Spike raster plot')
    plt.xlabel('time (ms)')
    plt.ylabel('Excitatory population: neuron number')
    plt.legend((dots, stars), ('Truely recalled neurons', 'Falsely recalled neurons'))

    plt.figure(figsize=(8, 6))
    N = len(t)
    T = 1.0 / len(t)
    xf = fftfreq(N, T)[:N//2]

    yf = fft(E)
    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Pyramidal cells')

#    yf = fft(B)
#    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Basket cells')

#    yf = fft(BS)
#    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Bistratified cells')

    plt.legend()
    plt.grid()
#    plt.savefig('ca1_theta_gamma_fft.png', dpi=500)
    plt.show()
'''
    # plt.figure()
    # N = len(t)
    # T = 1.0 / len(t)
    # xf = fftfreq(N, T)[:N//2]
    # yf = fft(E)
    # plt.plot(xf, 1.0 / 10 * np.abs(yf[0:N//2]))
    # plt.show()

    
'''
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

E_error = rmse(e_guess, E)
IP_error = rmse(ip_guess, HIPP)
I_error = rmse(i_guess, I)
S_error = rmse(s_guess, S)

print("RMSE of excitatory population = ", E_error)
print("RMSE of hippocampo-septal population = ", IP_error)
print("RMSE of inhibitory population = ", I_error)
print("RMSE of septal population = ", S_error)

fig = plt.figure(figsize=(8, 6))
fig.suptitle('Phase space diagrams')
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax1.plot3D(ip_guess, i_guess, e_guess)
ax1.set_xlabel('I_CA1P(t)', fontsize = 5.0)
ax1.set_ylabel('I_CA1I(t)', fontsize = 5.0)
ax1.set_zlabel('E_CA1(t)', fontsize = 5.0)
ax1.tick_params(axis='x', labelsize= 5.0)
ax1.tick_params(axis='y', labelsize= 5.0)
ax1.tick_params(axis='z', labelsize= 5.0)
ax1.view_init(-150, 60)

ax2 = fig.add_subplot(2, 2, 2, projection='3d')
ax2.plot3D(ip_guess, s_guess, e_guess)
ax2.set_xlabel('I_CA1P(t)', fontsize = 5.0)
ax2.set_ylabel('I_S(t)', fontsize = 5.0)
ax2.set_zlabel('E_CA1(t)', fontsize = 5.0)
ax2.tick_params(axis='x', labelsize= 5.0)
ax2.tick_params(axis='y', labelsize= 5.0)
ax2.tick_params(axis='z', labelsize= 5.0)
# ax2.view_init(-140, 60)

ax3 = fig.add_subplot(2, 2, 3, projection='3d')
ax3.plot3D(i_guess, s_guess, ip_guess)
ax3.set_xlabel('I_CA1I(t)', fontsize = 5.0)
ax3.set_ylabel('I_S(t)', fontsize = 5.0)
ax3.set_zlabel('I_CA1P(t)', fontsize = 5.0)
ax3.tick_params(axis='x', labelsize= 5.0)
ax3.tick_params(axis='y', labelsize= 5.0)
ax3.tick_params(axis='z', labelsize= 5.0)
ax3.view_init(-140, 60)

ax4 = fig.add_subplot(2, 2, 4, projection='3d')
ax4.plot3D(i_guess, s_guess, e_guess)
ax4.set_xlabel('I_CA1I(t)', fontsize = 5.0)
ax4.set_ylabel('I_S(t)', fontsize = 5.0)
ax4.set_zlabel('E_CA1(t)', fontsize = 5.0)
ax4.tick_params(axis='x', labelsize= 5.0)
ax4.tick_params(axis='y', labelsize= 5.0)
ax4.tick_params(axis='z', labelsize= 5.0)
ax4.view_init(-140, 60)

plt.subplots_adjust(hspace=0.12)
# plt.savefig('phase_space_diagrams.png', dpi=500)
plt.show()

'''
'''

    plt.figure()
    plt.title('E-I phase plot')
    plt.plot(E[0:200], S[0:200])
    plt.xlabel('E(t)')
    plt.ylabel('I(t)')

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 6), sharex=True)
ex_data = genfromtxt('results/ex.csv', delimiter='\t')
ex_data_50th = ex_data[:, 10]
ax1.scatter(t, ex_data_50th)

hipp_data = genfromtxt('results/hipp.csv', delimiter='\t')
hipp_data_50th = hipp_data[:, 10]
ax2.scatter(t, hipp_data_50th)

in_data = genfromtxt('results/in.csv', delimiter='\t')
in_data_50th = in_data[:, 10]
ax3.scatter(t, in_data_50th)

sept_data = genfromtxt('results/septum.csv', delimiter='\t')
sept_data_50th = sept_data[:, 10]
ax4.scatter(t, sept_data_50th)

plt.figure()
plt.title('E-I phase plot of the 50th neuron')
plt.plot(ex_data_50th, in_data_50th)
plt.xlabel('E(t) 50th neuron')
plt.ylabel('I(t) 50th neuron')

plt.tight_layout()
plt.show()


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
ax2.set_ylabel('External pseudo population: neuron number')
ax2.set_xlabel('time (ms)')
#plt.savefig('/home/ashraya/Desktop/1.png', dpi=250)

with open('results/ca_bin_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')

    dataR = genfromtxt('results/ca_bin_stats.csv', delimiter='\t')
    data = dataR.T

    data = np.delete(data, 0, axis=0)
    data = np.delete(data, (len(data)-1), axis=0)

    bins = []
    E = []
    EXT = []

    for row in reader:
        if (row[0] == "ex"):
            E = row[1:-1]
            E = [int(i) for i in E]
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

    ax2.bar(bins, EXT)
    ax2.plot(bins, EXT)
    ax2.set_xlabel('time, t(ms)')
    ax2.set_ylabel('External input')
    ax2.set_title('External input')

    plt.tight_layout()
    plt.show()
'''
