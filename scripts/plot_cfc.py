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
import seaborn as sns
from scipy.fftpack import fft, fftfreq
import os
import sys
from basic_units import radians, degrees, cos
from scipy import signal
from pylab import *
import scipy.stats

t = []
E = []
I = []
HIPP = []
S = []
B = []
BS = []
CA3 = []
PS = []
EC = []

def my_sin(x, freq, amp, phase):
    return amp * np.sin(((2 * np.pi * freq * x / 1000) + phase))

curPath = os.getcwd()
if ("scripts" in curPath):
    os.chdir("../")

with open('results/exp/ca_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    row = next(reader)

    data = genfromtxt('results/exp/ca_stats.csv', delimiter='\t')
    data = np.delete(data, (0), axis=0)

    t = data[:, 0]
    time_start = t[0]
    time_end = t[-1]

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
        guess_phase= 0.5
        guess_amplitude = 0.25
        guess_freq = 6.3

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
        guess_freq = 6.3

        p0=[guess_freq, guess_amplitude, guess_phase]

        ip_guess = my_sin(t, *p0)
        ax2.plot(t, ip_guess, color='green')
        ax2.grid()

    if (len(BS) > 0):
        ax3.plot(t, BS)
        ax3.set_ylim(0, 0.3)
        ax3.set_ylabel('I_BS(t)')
        guess_phase = 0.5
        guess_amplitude = 0.25
        guess_freq = 6.3

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
        guess_freq = 6.3
            
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
        guess_freq = 6.3
            
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
        guess_freq = 6.5
            
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

dt = 0.001
fNQ = 1 / dt / 2                     # ... and Nyquist frequency. 

# calculate for E

Wn = [5,7];                          # Set the passband [5-7] Hz,
n = 100;                             # ... and filter order,
b = signal.firwin(n, Wn, nyq=fNQ, pass_zero=False, window='hamming');
Vlo_E = signal.filtfilt(b, 1, E);    # ... and apply it to the data.

Wn = [50, 100];                      # Set the passband [50-100] Hz,
n = 100;                             # ... and filter order,
b = signal.firwin(n, Wn, nyq=fNQ, pass_zero=False, window='hamming');
Vhi_E = signal.filtfilt(b, 1, E);    # ... and apply it to the data.

figure(figsize=(14, 4))         # Create a figure with a specific size.
plot(t, E)                    # Plot the original data vs time.
plot(t, Vlo_E)                    # Plot the low-frequency filtered data vs time.
plot(t, Vhi_E)                    # Plot the high-frequency filtered data vs time.
xlabel('Time [s]')
plt.title('low and high frequency bands of E')
legend(['E', 'Vlo', 'Vhi']);  # Add a legend.

phi = angle(signal.hilbert(Vlo_E))     # Compute phase of low-freq signal
amp = abs(signal.hilbert(Vhi_E))       # Compute amplitude of high-freq signal
phi_hi = angle(signal.hilbert(Vhi_E))  # Compute phase of high-freq signal

p_bins = arange(-pi,pi,0.1)          # To compute CFC, define phase bins,
a_mean = zeros(size(p_bins)-1)       # ... variable to hold the amplitude,
p_mean = zeros(size(p_bins)-1)       # ... and variable to hold the phase.
h = max(a_mean)-min(a_mean)

for k in range(size(p_bins)-1):      # For each phase bin,
    pL = p_bins[k]                   #... get lower phase limit,
    pR = p_bins[k+1]                 #... get upper phase limit.
    indices=(phi_hi>=pL) & (phi_hi<pR)     #Find phases falling in this bin,
    a_mean[k] = mean(amp[indices])   #... compute mean amplitude,
    p_mean[k] = mean([pL, pR])       #... save center phase.

# calculate for B

Wn = [5,7];                          # Set the passband [5-7] Hz,
n = 100;                             # ... and filter order,
b = signal.firwin(n, Wn, nyq=fNQ, pass_zero=False, window='hamming');
Vlo_B = signal.filtfilt(b, 1, B);    # ... and apply it to the data.

Wn = [50, 100];                      # Set the passband [50-100] Hz,
n = 100;                             # ... and filter order,
b = signal.firwin(n, Wn, nyq=fNQ, pass_zero=False, window='hamming');
Vhi_B = signal.filtfilt(b, 1, B);    # ... and apply it to the data.

figure(figsize=(14, 4))         # Create a figure with a specific size.
plot(t, B)                    # Plot the original data vs time.
plot(t, Vlo_B)                    # Plot the low-frequency filtered data vs time.
plot(t, Vhi_B)                    # Plot the high-frequency filtered data vs time.
xlabel('Time [s]')
plt.title('low and high frequency bands of B')
legend(['B', 'Vlo', 'Vhi']);  # Add a legend.

phi_b = angle(signal.hilbert(Vlo_B))     # Compute phase of low-freq signal
amp_b = abs(signal.hilbert(Vhi_B))       # Compute amplitude of high-freq signal
phi_hi_b = angle(signal.hilbert(Vhi_B))  # Compute phase of high-freq signal

p_bins_b = arange(-pi,pi,0.1)          # To compute CFC, define phase bins,
a_mean_b = zeros(size(p_bins)-1)       # ... variable to hold the amplitude,
p_mean_b = zeros(size(p_bins)-1)       # ... and variable to hold the phase.
h = max(a_mean_b)-min(a_mean_b)

for k in range(size(p_bins_b)-1):      # For each phase bin,
    pL_b = p_bins_b[k]                   #... get lower phase limit,
    pR_b = p_bins_b[k+1]                 #... get upper phase limit.
    indices_b=(phi_hi_b>=pL_b) & (phi_hi_b<pR_b)     #Find phases falling in this bin,
    a_mean_b[k] = mean(amp_b[indices_b])   #... compute mean amplitude,
    p_mean_b[k] = mean([pL_b, pR_b])       #... save center phase.

# calculate for BS

Wn = [5,7];                          # Set the passband [5-7] Hz,
n = 100;                             # ... and filter order,
b = signal.firwin(n, Wn, nyq=fNQ, pass_zero=False, window='hamming');
Vlo_BS = signal.filtfilt(b, 1, BS);    # ... and apply it to the data.

Wn = [50, 100];                      # Set the passband [50-100] Hz,
n = 100;                             # ... and filter order,
b = signal.firwin(n, Wn, nyq=fNQ, pass_zero=False, window='hamming');
Vhi_BS = signal.filtfilt(b, 1, BS);    # ... and apply it to the data.

figure(figsize=(14, 4))         # Create a figure with a specific size.
plot(t, BS)                    # Plot the original data vs time.
plot(t, Vlo_BS)                    # Plot the low-frequency filtered data vs time.
plot(t, Vhi_BS)                    # Plot the high-frequency filtered data vs time.
xlabel('Time [s]')
plt.title('low and high frequency bands of BS')
legend(['BS', 'Vlo', 'Vhi']);  # Add a legend.

phi_bs = angle(signal.hilbert(Vlo_BS))     # Compute phase of low-freq signal
amp_bs = abs(signal.hilbert(Vhi_BS))       # Compute amplitude of high-freq signal
phi_hi_bs = angle(signal.hilbert(Vhi_BS))  # Compute phase of high-freq signal

p_bins_bs = arange(-pi,pi,0.1)          # To compute CFC, define phase bins,
a_mean_bs = zeros(size(p_bins_bs)-1)       # ... variable to hold the amplitude,
p_mean_bs = zeros(size(p_bins_bs)-1)       # ... and variable to hold the phase.
h = max(a_mean_bs)-min(a_mean_bs)

for k in range(size(p_bins_bs)-1):      # For each phase bin,
    pL_bs = p_bins_bs[k]                   #... get lower phase limit,
    pR_bs = p_bins_bs[k+1]                 #... get upper phase limit.
    indices_bs=(phi_hi_bs>=pL_bs) & (phi_hi_bs<pR_bs)     #Find phases falling in this bin,
    a_mean_bs[k] = mean(amp_bs[indices_bs])   #... compute mean amplitude,
    p_mean_bs[k] = mean([pL_bs, pR_bs])       #... save center phase.

phi = np.rad2deg(phi)
phi_b = np.rad2deg(phi_b)
phi_bs = np.rad2deg(phi_bs)

phi_hi = np.rad2deg(phi_hi)
phi_hi_b = np.rad2deg(phi_hi_b)
phi_hi_bs = np.rad2deg(phi_hi_bs)


# Plot of phase distribution
plt.figure()
plt.title('Plot with 7Hz phase')
plt.plot(t[50:350], phi[50:350], color='red', alpha=0.5, label='E')
plt.plot(t[50:350], phi_b[50:350], color='blue', alpha=0.5, label='B')
plt.plot(t[50:350], phi_bs[50:350], color='green', alpha=0.5, label='BS')
plt.xlabel('time(ms)')
plt.ylabel('phase(deg)')
plt.legend()


#fig, axes = plt.subplots(3, 1, sharex=True)
#fig.suptitle('Phase distribution low frequency')
#axes[0].set_title('Pyramidal population')
#axes[1].set_title('Basket population')
#axes[2].set_title('Bistratified population')
#sns.distplot(phi[0:250], label='E', ax=axes[0], kde=True, bins=50)
#sns.distplot(phi_b[0:250], label='B', ax=axes[1], kde=True, bins=50)
#sns.distplot(phi_bs[0:250], label='BS', ax=axes[2], kde=True, bins=50)
#sns.distplot(phi[0:250], label='E', ax=axes[0], kde=True, hist=False)
#sns.distplot(phi_b[0:250], label='B', ax=axes[1], kde=True, hist=False)
#sns.distplot(phi_bs[0:250], label='BS', ax=axes[2], kde=True, hist=False)
#plt.xticks([-120, -110, -100, -90, -80, -70, -60, -50, 0, 50, 100, 150])


'''
# Plot of CFC
plt.figure()
plot(p_mean, a_mean, color='blue', label='Pyramidal')                 #Plot the phase versus amplitude,
plot(p_mean_b, a_mean_b, color='red', label='Basket')                 #Plot the phase versus amplitude,
plot(p_mean_bs, a_mean_bs, color='green', label='Bistratified')       #Plot the phase versus amplitude,
ylabel('High-frequency amplitude')   #... with axes labeled.
xlabel('Low-frequency phase')
title('CFC');
legend();
grid();

width = np.pi / 4 * np.random.rand(len(E))
colors = cm.rainbow(np.linspace(0, 1, 6))

plt.figure()
plt.title('Plot with 50-60Hz phase')
plt.plot(phi_hi, t, color='red', alpha=0.5, label='E')
plt.plot(phi_hi_b, t, color='blue', alpha=0.5, label='B')
plt.plot(phi_hi_bs, t, color='green', alpha=0.5, label='BS')
plt.xlabel('time(ms)')
plt.ylabel('phase(rad)')
plt.legend()

fig, axes = plt.subplots(3, 1, sharex=True)
fig.suptitle('Phase distribution low frequency')
axes[0].set_title('Pyramidal population')
axes[1].set_title('Basket population')
axes[2].set_title('Bistratified population')
sns.distplot(phi[0:250], label='E', ax=axes[0], kde=True, bins=50)
sns.distplot(phi_b[0:250], label='B', ax=axes[1], kde=True, bins=50)
sns.distplot(phi_bs[0:250], label='BS', ax=axes[2], kde=True, bins=50)
plt.xticks([-120, -110, -100, -90, -80, -70, -60, -50, 0, 50, 100, 150])

fig, axes = plt.subplots(3, 1, sharex=True)
fig.suptitle('Phase distribution high frequency')
axes[0].set_title('Pyramidal population')
axes[1].set_title('Basket population')
axes[2].set_title('Bistratified population')
sns.distplot(phi_hi, label='E', ax=axes[0], kde=True, bins=50)
sns.distplot(phi_hi_b, label='B', ax=axes[1], kde=True, bins=50)
sns.distplot(phi_hi_bs, label='BS', ax=axes[2], kde=True, bins=50)
plt.xticks([-120, -110, -100, -90, -80, -70, -60, -50, 0, 50, 100, 150])

plt.figure()
ax = plt.subplot(projection='polar')
ax.set_title('polar E plot with 7Hz phase')
ax.bar(E, phi, width=width, bottom=0.0, color='red', alpha=0.5)

plt.figure()
ax = plt.subplot(projection='polar')
ax.set_title('polar B plot with 7Hz phase')
ax.bar(B, phi_b, width=width, bottom=0.0, color='green', alpha=0.5)

plt.figure()
ax = plt.subplot(projection='polar')
ax.set_title('polar BS plot with 7Hz phase')
ax.bar(BS, phi_bs, width=width, bottom=0.0, color='blue', alpha=0.5)

plt.figure()
ax = plt.subplot(projection='polar')
ax.set_title('polar E plot with 60Hz phase')
ax.bar(E, phi_hi, width=width, bottom=0.0, color='red', alpha=0.5)

plt.figure()
ax = plt.subplot(projection='polar')
ax.set_title('polar B plot with 60Hz phase')
ax.bar(B, phi_hi_b, width=width, bottom=0.0, color='green', alpha=0.5)

plt.figure()
ax = plt.subplot(projection='polar')
ax.set_title('polar BS plot with 60Hz phase')
ax.bar(BS, phi_hi_bs, width=width, bottom=0.0, color='blue', alpha=0.5)

plt.figure()
p = scipy.stats.norm.pdf(p_mean_b)
ax = sns.histplot(p, bins=100)
ax.set(ylabel='Probability Distribution of B phases', xlabel='Frequency')

plt.figure()
p = scipy.stats.norm.pdf(p_mean_bs)
ax = sns.histplot(p, bins=100)
ax.set(ylabel='Probability Distribution of BS phases', xlabel='Frequency')

plt.figure()
p = scipy.stats.norm.pdf(p_mean)
ax = sns.histplot(p, bins=100)
ax.set(ylabel='Probability Distribution of E phases', xlabel='Frequency')
'''

plt.show()
