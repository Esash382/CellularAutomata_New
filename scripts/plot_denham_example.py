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
    HIPP = []
    S = []
    CA3 = []
    PS = []

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
            elif (name == "ps"):
                PS = data[:, i]
            else:
                CA3 = data[:, i]

    # Plot active neuron stats
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, figsize=(8, 6), sharex=True)
    ax1.title.set_text('Cellular automata simulation of septal pacemaker circuit')
    if (len(E) > 0):
        ax1.plot(t, E)
        ax1.set_ylabel('E_CA1(t)')
#        ax1.set_ylim(0, 0.25)
        guess_fs = 1000
        guess_phase= 0
        guess_amplitude = 0.25
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase, guess_fs]

        e_guess = my_sin(t, *p0)
#        ax1.plot(t, e_guess, color='green')

    if (len(HIPP) > 0):
        ax2.plot(t, HIPP)
        ax2.set_ylabel('I_CA1P(t)')
        
#        ax2.set_ylim(0, 0.25)
        guess_fs = 1010
        guess_phase = -200
        guess_amplitude = 0.2
        guess_freq = 6.5
            
        p0=[guess_freq, guess_amplitude, guess_phase, guess_fs]

        ip_guess = my_sin(t, *p0)
#        ax2.plot(t, ip_guess, color='green')

    if (len(I) > 0):
        ax3.plot(t, I)
#        ax3.set_ylim(0, 0.25)
        ax3.set_ylabel('I_CA1I(t)')
        
        guess_fs = 1035
        guess_phase = -1800
        guess_amplitude = 0.2
        guess_freq = 6.5
            
        p0=[guess_freq, guess_amplitude, guess_phase, guess_fs]

        i_guess = my_sin(t, *p0)
#        ax3.plot(t, i_guess, color='green')

    if (len(S) > 0):
        ax4.plot(t, S)
#        ax4.set_ylim(0, 0.25)
        ax4.set_ylabel('I_S(t)')
        
        guess_fs = 1100
        guess_phase = 2200
        guess_amplitude = 0.2
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase, guess_fs]

        s_guess = my_sin(t, *p0)
#        ax4.plot(t, s_guess, color='green')

    if (len(CA3) > 0):
        ax5.plot(t, CA3)
        ax5.set_ylabel('CA3')

    if (len(PS) > 0):
        ax6.plot(t, PS)
        ax6.set_xlabel('time, t(ms)')
        ax6.set_ylabel('PS')

    plt.tight_layout()
    plt.rcParams.update({'font.size': 22})
    plt.savefig('denham_results.png', dpi=300)
#    plt.show()


    # FFT
    N = len(t)
    T = 1.0 / len(t)
    xf = fftfreq(N, T)[:N//2]

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 6), sharex=True)

    yf = fft(E)
    ax1.set_ylabel('E_CA1(t)')
    ax1.plot(xf[:10], 1.0 / 10 * np.abs(yf[0:N//100]))

    yf = fft(HIPP)
    ax2.set_ylabel('I_CA1P(t)')
    ax2.plot(xf[:10], 1.0 / 10 * np.abs(yf[0:N//100]))

    yf = fft(I)
    ax3.set_ylabel('I_CA1I(t)')
    ax3.plot(xf[:10], 1.0 / 10 * np.abs(yf[0:N//100]))

    yf = fft(S)
    ax4.set_ylabel('S(t)')
    ax4.set_xlabel('time, t(ms)')
    plt.plot(xf[:10], 1.0 / 10 * np.abs(yf[0:N//100]))

    plt.xlabel('Frequency')
    plt.savefig('denham_results_fft.png', dpi=300)
    plt.show()


    
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
plt.savefig('/home/ashraya/Desktop/1.png', dpi=250)

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
