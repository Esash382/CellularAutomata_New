#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 18:45:39 2021

@author: ashraya
"""

from numpy import genfromtxt
import csv
import numpy as np
import matplotlib.pyplot as plt

t = []

guess_freq = 1
guess_amplitude = 0.1
guess_phase = 0
guess_offset = 0
guess_fs = 1

def my_sin(x, freq, amp, phase, fs):
    return amp * np.sin(((2 * np.pi * freq * x / 1000) + phase))

with open('../results/ca_stats.csv') as f:
    reader = csv.reader(f, delimiter='\t')
    row = next(reader)

    data = genfromtxt('../results/ca_stats.csv', delimiter='\t')
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
        ax1.set_ylim(0, 0.3)
        guess_fs = 1000
        guess_phase= 0
        guess_amplitude = 0.2
        guess_freq = 7

        p0=[guess_freq, guess_amplitude, guess_phase, guess_fs]

        e_guess = my_sin(t, *p0)
        ax1.plot(t, e_guess, color='green')

    if (len(HIPP) > 0):
        ax2.plot(t, HIPP)
        ax2.set_ylabel('I_CA1P(t)')
        
        ax2.set_ylim(0, 0.25)
        guess_fs = 1010
        guess_phase = -0.5
        guess_amplitude = 0.2
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase, guess_fs]

        ip_guess = my_sin(t, *p0)
        ax2.plot(t, ip_guess, color='green')

    if (len(I) > 0):
        ax3.plot(t, I)
        ax3.set_ylim(0, 0.3)
        ax3.set_ylabel('I_CA1I(t)')
        
        guess_fs = 1045
        guess_phase = -91
        guess_amplitude = 0.25
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase, guess_fs]

        i_guess = my_sin(t, *p0)
        ax3.plot(t, i_guess, color='green')

    if (len(S) > 0):
        ax4.plot(t, S)
        ax4.set_ylim(0, 0.25)
        ax4.set_ylabel('I_S(t)')
        
        guess_fs = 1100
        guess_phase = -92.5
        guess_amplitude = 0.2
        guess_freq = 7
            
        p0=[guess_freq, guess_amplitude, guess_phase, guess_fs]

        s_guess = my_sin(t, *p0)
        ax4.plot(t, s_guess, color='green')

    if (len(CA3) > 0):
        ax5.plot(t, CA3)
        ax5.set_ylabel('CA3')

    if (len(PS) > 0):
        ax6.plot(t, PS)
        ax6.set_xlabel('time, t(ms)')
        ax6.set_ylabel('PS')

    plt.tight_layout()
    # plt.savefig('curve_fit_ca_denham_results.png', dpi=500)
    plt.show()
    
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

# E_error = rmse(e_guess, E)
# IP_error = rmse(ip_guess, HIPP)
# I_error = rmse(i_guess, I)
# S_error = rmse(s_guess, S)

# print("RMSE of excitatory population = ", E_error)
# print("RMSE of hippocampo-septal population = ", IP_error)
# print("RMSE of inhibitory population = ", I_error)
# print("RMSE of septal population = ", S_error)

# e1_guess = [0 if x < 0 else x for x in e_guess]
# ip1_guess = [0 if x < 0 else x for x in ip_guess]
# i1_guess = [0 if x < 0 else x for x in i_guess]
# s1_guess = [0 if x < 0 else x for x in s_guess]

# fig = plt.figure(figsize=(8, 6))
# fig.suptitle('Phase space diagrams')
# ax1 = fig.add_subplot(2, 2, 1, projection='3d')
# ax1.plot3D(ip1_guess, i1_guess, e1_guess)
# ax1.set_xlabel('I_CA1P(t)', fontsize = 5.0)
# ax1.set_ylabel('I_CA1I(t)', fontsize = 5.0)
# ax1.set_zlabel('E_CA1(t)', fontsize = 5.0)
# ax1.set_xlim(0, 0.5)
# ax1.set_ylim(0, 0.5)
# ax1.set_zlim(0, 0.5)
# ax1.tick_params(axis='x', labelsize= 5.0)
# ax1.tick_params(axis='y', labelsize= 5.0)
# ax1.tick_params(axis='z', labelsize= 5.0)
# ax1.view_init(-150, 50)

# ax2 = fig.add_subplot(2, 2, 2, projection='3d')
# ax2.plot3D(ip1_guess, s1_guess, e1_guess)
# ax2.set_xlabel('I_CA1P(t)', fontsize = 5.0)
# ax2.set_ylabel('I_S(t)', fontsize = 5.0)
# ax2.set_zlabel('E_CA1(t)', fontsize = 5.0)
# ax2.set_xlim(0, 0.5)
# ax2.set_ylim(0, 0.5)
# ax2.set_zlim(0, 0.5)
# ax2.tick_params(axis='x', labelsize= 5.0)
# ax2.tick_params(axis='y', labelsize= 5.0)
# ax2.tick_params(axis='z', labelsize= 5.0)
# ax2.view_init(-150, 50)

# ax3 = fig.add_subplot(2, 2, 3, projection='3d')
# ax3.plot3D(i1_guess, s1_guess, ip1_guess)
# ax3.set_xlabel('I_CA1I(t)', fontsize = 5.0)
# ax3.set_ylabel('I_S(t)', fontsize = 5.0)
# ax3.set_zlabel('I_CA1P(t)', fontsize = 5.0)
# ax3.set_xlim(0, 0.5)
# ax3.set_ylim(0, 0.5)
# ax3.set_zlim(0, 0.5)
# ax3.tick_params(axis='x', labelsize= 5.0)
# ax3.tick_params(axis='y', labelsize= 5.0)
# ax3.tick_params(axis='z', labelsize= 5.0)
# ax3.view_init(-150, 50)

# ax4 = fig.add_subplot(2, 2, 4, projection='3d')
# ax4.plot3D(i1_guess, s1_guess, e1_guess)
# ax4.set_xlabel('I_CA1I(t)', fontsize = 5.0)
# ax4.set_ylabel('I_S(t)', fontsize = 5.0)
# ax4.set_zlabel('E_CA1(t)', fontsize = 5.0)
# ax4.set_xlim(0, 0.5)
# ax4.set_ylim(0, 0.5)
# ax4.set_zlim(0, 0.5)
# ax4.tick_params(axis='x', labelsize= 5.0)
# ax4.tick_params(axis='y', labelsize= 5.0)
# ax4.tick_params(axis='z', labelsize= 5.0)
# ax4.view_init(-150, 50)

# plt.subplots_adjust(hspace=0.12)
# # plt.savefig('phase_space_diagrams.png', dpi=500)
# plt.show()