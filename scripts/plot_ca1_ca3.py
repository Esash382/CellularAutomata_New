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

def my_sin(x, freq, amp, phase):
    return amp * np.sin(((2 * np.pi * freq * x / 1000) + phase))

with open('results/ca_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    row = next(reader)

    data = genfromtxt('results/ca_stats.csv', delimiter='\t')
    data = np.delete(data, (0), axis=0)

    t = data[:, 0]
    time_start = t[0]
    time_end = t[-1]

    E3 = []
    I3 = []
    HIPP3 = []
    B3 = []
    BS3 = []

    E1 = []
    I1 = []
    HIPP1 = []
    B1 = []
    BS1 = []

    S = []
    DG = []
    PS = []
    EC = []

    for i in range(len(row)):
        if (row[i].find("_active") > 0):
            name = row[i][:row[i].find('_')]
            if (name == "ex3"):
                E3 = data[:, i]
            elif (name == "in3"):
                I3 = data[:, i]
            elif (name == "hipp3"):
                HIPP3 = data[:, i]
            elif (name == "bas3"):
                B3 = data[:, i]
            elif (name == "bis3"):
                BS3 = data[:, i]

            elif (name == "ex1"):
                E1 = data[:, i]
            elif (name == "in1"):
                I1 = data[:, i]
            elif (name == "hipp1"):
                HIPP1 = data[:, i]
            elif (name == "bas1"):
                B1 = data[:, i]
            elif (name == "bis1"):
                BS1 = data[:, i]

            elif (name == "septum"):
                S = data[:, i]
            elif (name == "ec"):
                EC = data[:, i]
            elif (name == "dg"):
                DG = data[:, i]
            elif (name == "ps"):
                PS = data[:, i]

    # Plot active neuron stats
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9) = plt.subplots(9, 1, figsize=(10, 10), sharex=True)
    ax1.title.set_text('Cellular automata simulation of CA3 circuit')
    if (len(E3) > 0):
        ax1.plot(t, E3)
        ax1.set_ylim(0, 0.25)
        ax1.set_ylabel('E_CA3(t)')
        guess_phase= 6.5
        guess_amplitude = 0.2
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        e_guess = my_sin(t, *p0)
        ax1.plot(t, e_guess, color='green')

    if (len(B3) > 0):
        ax2.plot(t, B3)
        ax2.set_ylim(0, 0.4)
        ax2.set_ylabel('I_B(t)')
        guess_phase = 5
        guess_amplitude = 0.2
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        ip_guess = my_sin(t, *p0)
        ax2.plot(t, ip_guess, color='green')

    if (len(BS3) > 0):
        ax3.plot(t, BS3)
        ax3.set_ylim(0, 0.055)
        ax3.set_ylabel('I_BS(t)')
        guess_phase = 6
        guess_amplitude = 0.05
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        i_guess = my_sin(t, *p0)
        ax3.plot(t, i_guess, color='green')

    if (len(I3) > 0):
        ax4.plot(t, I3)
        ax4.set_ylim(0, 0.4)
        ax4.set_ylabel('I_CA3I(t)')
        guess_phase = 4
        guess_amplitude = 0.3
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        i_guess = my_sin(t, *p0)
        ax4.plot(t, i_guess, color='green')
 
    if (len(HIPP3) > 0):
        ax5.plot(t, HIPP3)
        ax5.set_ylim(0, 0.3)
        ax5.set_ylabel('I_CA3P(t)')
        guess_phase = 6
        guess_amplitude = 0.25
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        i_guess = my_sin(t, *p0)
        ax5.plot(t, i_guess, color='green')

    if (len(S) > 0):
        ax6.plot(t, S)
        ax6.set_ylim(0, 0.3)
        ax6.set_ylabel('I_S(t)')
        guess_phase = 8
        guess_amplitude = 0.25
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        i_guess = my_sin(t, *p0)
        ax6.plot(t, i_guess, color='green')
        
    if (len(EC) > 0):
        ax7.plot(t, EC)
        ax7.set_ylabel('EC')

    if (len(DG) > 0):
        ax8.plot(t, DG)
        ax8.set_ylabel('DG')

    if (len(PS) > 0):
        ax9.plot(t, PS)
        ax9.set_ylabel('PS')
        ax9.set_xlabel('time, t(ms)')

    plt.tight_layout()
    # plt.savefig('ca1_ca3_ca3_theta_gamma.png', dpi=500)

    # FFT
    N = len(t)
    T = 1.0 / len(t)
    yf = fft(E3)
    xf = fftfreq(N, T)[:N//2]
    plt.figure()
    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Pyramidal cells')

    yf = fft(B3)
    xf = fftfreq(N, T)[:N//2]
    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Basket cells')

    yf = fft(BS3)
    xf = fftfreq(N, T)[:N//2]
    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Bistratified cells')
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.title('FFT of CA3')
    plt.grid()
    plt.legend()
    # plt.savefig('ca1_ca3_ca3_theta_gamma_fft.png', dpi=500)
    plt.show()

    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(8, 1, figsize=(10, 10), sharex=True)
    ax1.title.set_text('Cellular automata simulation of CA1 circuit')
    if (len(E1) > 0):
        ax1.plot(t, E1)
        ax1.set_ylim(0, 0.25)
        ax1.set_ylabel('E_CA1(t)')
        guess_phase= 7.5
        guess_amplitude = 0.2
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        e_guess = my_sin(t, *p0)
        ax1.plot(t, e_guess, color='green')

    if (len(B1) > 0):
        ax2.plot(t, B1)
        ax2.set_ylim(0, 0.3)
        ax2.set_ylabel('I_B(t)')
        guess_phase = 5
        guess_amplitude = 0.2
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        ip_guess = my_sin(t, *p0)
        ax2.plot(t, ip_guess, color='green')

    if (len(BS1) > 0):
        ax3.plot(t, BS1)
        ax3.set_ylim(0, 0.2)
        ax3.set_ylabel('I_BS(t)')
        guess_phase = 5.5
        guess_amplitude = 0.1
        guess_freq = 6.5

        p0=[guess_freq, guess_amplitude, guess_phase]

        i_guess = my_sin(t, *p0)
        ax3.plot(t, i_guess, color='green')

    if (len(I1) > 0):
        ax4.plot(t, I1)
        ax4.set_ylim(0, 0.3)
        ax4.set_ylabel('I_CA1I(t)')
        guess_phase = 4
        guess_amplitude = 0.25
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax4.plot(t, s_guess, color='green')

    if (len(HIPP1) > 0):
        ax5.plot(t, HIPP1)
        ax5.set_ylim(0, 0.3)
        ax5.set_ylabel('I_CA1P(t)')
        guess_phase = 6
        guess_amplitude = 0.25
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax5.plot(t, s_guess, color='green')
        
    if (len(S) > 0):
        ax6.plot(t, S)
        ax6.set_ylim(0, 0.25)
        ax6.set_ylabel('I_S(t)')
        guess_phase = 2
        guess_amplitude = 0.2
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase]

        s_guess = my_sin(t, *p0)
        ax6.plot(t, s_guess, color='green')

    if (len(EC) > 0):
        ax7.plot(t, EC)
        ax7.set_ylabel('EC')

    if (len(PS) > 0):
        ax8.plot(t, PS)
        ax8.set_ylabel('PS')
        ax8.set_xlabel('time, t(ms)')

    plt.tight_layout()
    # plt.savefig('ca1_ca3_ca1_theta_gamma.png', dpi=500)

    plt.figure(figsize=(8, 6))
    N = len(t)
    T = 1.0 / len(t)
    xf = fftfreq(N, T)[:N//2]

    yf = fft(E1)
    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Pyramidal cells')

    yf = fft(B1)
    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Basket cells')

    yf = fft(BS1)
    plt.plot(xf[:100], 1.0 / 10 * np.abs(yf[0:N//10]), label='Bistratified cells')

    plt.legend()
    plt.grid()
    # plt.savefig('ca1_ca3_ca1_theta_gamma_fft.png', dpi=500)
    plt.show()
