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
import matplotlib.cm as cm
from scipy.fftpack import fft, fftfreq
import os
import sys
from scipy.signal import butter, lfilter
from scipy.signal import freqz
from basic_units import radians, degrees, cos

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
        guess_freq = 6.4

        p0=[guess_freq, guess_amplitude, guess_phase]

        e_guess = my_sin(t, *p0)
        ax1.plot(t, e_guess, color='green')
        ax1.grid()

    if (len(B) > 0):
        ax2.plot(t, B)
        ax2.set_ylim(0, 0.3)
        ax2.set_ylabel('I_B(t)')
        guess_phase = -2.5
        guess_amplitude = 0.2
        guess_freq = 6.4

        p0=[guess_freq, guess_amplitude, guess_phase]

        ip_guess = my_sin(t, *p0)
        ax2.plot(t, ip_guess, color='green')
        ax2.grid()

    if (len(BS) > 0):
        ax3.plot(t, BS)
        ax3.set_ylim(0, 0.3)
        ax3.set_ylabel('I_BS(t)')
        guess_phase = 0
        guess_amplitude = 0.25
        guess_freq = 6.4

        p0=[guess_freq, guess_amplitude, guess_phase]

        i_guess = my_sin(t, *p0)
        ax3.plot(t, i_guess, color='green')
        ax3.grid()

    if (len(I) > 0):
        ax4.plot(t, I)
        ax4.set_ylim(0, 0.3)
        ax4.set_ylabel('I_CA1I(t)')
        guess_phase = -2.5
        guess_amplitude = 0.25
        guess_freq = 6.4
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax4.plot(t, s_guess, color='green')
        ax4.grid()

    if (len(HIPP) > 0):
        ax5.plot(t, HIPP)
        ax5.set_ylim(0, 0.3)
        ax5.set_ylabel('I_CA1P(t)')
        guess_phase = -0.2
        guess_amplitude = 0.25
        guess_freq = 6.4
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax5.plot(t, s_guess, color='green')
        ax5.grid()
        
    if (len(S) > 0):
        ax6.plot(t, S)
        ax6.set_ylim(0, 0.25)
        ax6.set_ylabel('I_S(t)')
        guess_phase = 1.3
        guess_amplitude = 0.2
        guess_freq = 6.4
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax6.plot(t, s_guess, color='green')
        ax6.grid()

    if (len(EC) > 0):
        ax7.plot(t, EC)
        ax7.set_ylabel('EC')
        ax7.grid()

    if (len(CA3) > 0):
        ax8.plot(t, CA3)
        ax8.set_ylabel('CA3')
        ax8.grid()

    if (len(PS) > 0):
        ax9.plot(t, PS)
        ax9.set_ylabel('PS')
        ax9.set_xlabel('time, t(ms)')
        ax9.grid()

    plt.tight_layout()
#    plt.savefig('figs/recall_fin/nonintersecting patterns/1_ca1_nonint_activity.png', dpi=500)
#    plt.show()


#    dt = 0.01
#    Fs = 1 / dt
#    plt.figure()
#    plt.grid()
#    spec, freq, _ = plt.phase_spectrum(E, color ='green', Fs=Fs, Fc=0)

    plt.figure(figsize=(8, 3))
    dataR = genfromtxt('results/ex.csv', delimiter='\t')
    dataPRand = genfromtxt('results/ca_p_rand_stats.csv', delimiter='\t')
    dataT = dataR.T
    dataS = np.delete(dataT, 0, axis=0)
    dataS[dataS == -1] = np.nan
    dataS1 = dataS.copy()
    dataS2 = dataS.copy()
    dataS3 = dataS.copy()
    dataS4 = dataS.copy()
    dataS5 = dataS.copy()
    dataS6 = dataS.copy()

    colors = cm.rainbow(np.linspace(0, 1, 6))

    for i in range(len(dataS)):
        for j in range(len(dataS[i])):
            if (np.isnan(dataS[i][j])):
                continue
            for p in range(len(dataPRand)):
                if (dataS[i][j] not in dataPRand[p]): # have only recalled neurons
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
                else:
                    dataS6[i][j] = np.nan   # have only spurious active neurons

        plt.scatter(t, dataS1[i], marker="o", s=5, color=colors[0])
        plt.scatter(t, dataS2[i], marker="o", s=5, color=colors[1])
        plt.scatter(t, dataS3[i], marker="o", s=5, color=colors[2])
        plt.scatter(t, dataS4[i], marker="o", s=5, color=colors[3])
        plt.scatter(t, dataS5[i], marker="o", s=5, color=colors[4])
        plt.scatter(t, dataS6[i], marker="o", s=50, color='red')

    plt.title('Spike raster plot')
    plt.xlabel('time (ms)')
    plt.ylabel('Excitatory population: neuron number')
    plt.ylim(0, 100)
#    plt.savefig('figs/recall_fin/nonintersecting patterns/2_ca1_nonint_raster.png', dpi=500)
    # plt.legend((dots, stars), ('Truely recalled neurons', 'Falsely recalled neurons'))

plt.figure(figsize=(8, 6))
N = len(t)
T = 1.0 / len(t)
xf = fftfreq(N, T)[:N//2]

yf = fft(E)
plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Pyramidal cells')

#yf = fft(B)
#plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Basket cells')

#yf = fft(BS)
#plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Bistratified cells')

plt.xlabel("Frequency (Hz)")
plt.legend()
plt.grid()
#plt.savefig('figs/recall_fin/nonintersecting patterns/3_ca1_nonint_fft.png', dpi=500)

# bar plot for % recall of patterns
plt.figure(figsize=(8, 6))
data = genfromtxt('results/recall_percentage.csv', dtype=int, delimiter='\t')
last_data = []
for d in data:
    if (d[0] == 1000):
        last_data.append(d)

last_data = np.array(last_data)
pattern_index = last_data[:, 1]
recall_percent = last_data[:, 2]
plt.xlabel('pattern index')
plt.ylabel('% recall')
plt.xticks(pattern_index, pattern_index)
plt.bar(pattern_index, recall_percent)
#plt.savefig('figs/recall_fin/nonintersecting patterns/4_ca1_nonint_bar.png', dpi=500)

# bar plot for recall correlation
plt.figure(figsize=(8, 6))
data = genfromtxt('results/recall_correlation.csv', dtype=float, delimiter='\t')
last_data = []
for d in data:
    last_data.append(d)

last_data = np.array(last_data)
pattern_index = last_data[:, 0]
recall_percent = last_data[:, 1]
plt.xlabel('pattern index')
plt.ylabel('recall correlation')
plt.xticks(pattern_index, pattern_index)
plt.bar(pattern_index, recall_percent)
#plt.savefig('figs/recall_fin/nonintersecting patterns/4_ca1_nonint_recall_corr.png', dpi=500)


plt.show()
