#--PlotFits.py - Version 1 - 04/02/2016
#--Author: Bradley J Kavanagh
#--Summary: Code for plotting BG-only and signal+BG
#--fits to the ATLAS 750 GeV excess
#--Run ATLASfits.py with saveResults=1 to get the necessary files
#--Change 'include_signal' lower down in the options to
#--plot with or without signal
#--Please report any problems to: bradley.kavanagh@lpthe.jussieu.fr

print "----Plotting fits to ATLAS diphoton data---"

import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
import sys

from scipy.stats import chi2
font = {'family' : 'sans-serif',
        'size'   : 16}
mpl.rc('font', **font)

from ATLASfits_utils import getBestFit, f_ATLAS0, f_ATLAS1, f_wide

#---Functions---
#Calculate Poisson CLs on the data
def calcError1(k, CL):
    return k - 0.5*chi2.ppf((1-CL)/2, 2*k)
def calcError2(k, CL):
    return 0.5*chi2.ppf((1+CL)/2.0, 2*k + 2) - k


#---Options---
#Include signal in the fits (set to 0 for BG-only)
include_signal = 0


#---Load in the ATLAS data and calculate bins
data = np.loadtxt("ATLASdata1.txt")[:,1]
nbins = data.size

#Calculate bin edges
m_edges = np.linspace(150,150+40*(nbins),nbins+1)
m_centres = m_edges[:-1] + 20


#---Main Procedure----
m = np.linspace(200,1850,100)

#Loading best fit values
sig_str = 'BG-only'
if (include_signal): sig_str = 'wide'

bf_A0 = np.loadtxt('Fits_BG=0_'+sig_str+'.txt')
bf_A1 = np.loadtxt('Fits_BG=1_'+sig_str+'.txt')
bf_A0N = np.loadtxt('Fits_BG=3_'+sig_str+'.txt')
bf_A1N = np.loadtxt('Fits_BG=4_'+sig_str+'.txt')

#Calculate BG curves for plotting
y_A0 = f_ATLAS0(m, bf_A0[0], bf_A0[1])
y_A1 = f_ATLAS1(m, bf_A1[0], bf_A1[1], bf_A1[2])
y_A0N = 10**bf_A0N[0]*f_ATLAS0(m, bf_A0N[1], bf_A0N[2])
y_A1N = 10**bf_A1N[0]*f_ATLAS1(m, bf_A1N[1], bf_A1N[2], bf_A1N[3])

#Add signal (if you fancy it)
if (include_signal):
    y_A0 += 40*10**bf_A0[-1]*f_wide(m, bf_A0[-3], bf_A0[-2]/100.0)
    y_A1 += 40*10**bf_A1[-1]*f_wide(m, bf_A1[-3], bf_A1[-2]/100.0)
    y_A0N += 40*10**bf_A0N[-1]*f_wide(m, bf_A0N[-3], bf_A0N[-2]/100.0)
    y_A1N += 40*10**bf_A1N[-1]*f_wide(m, bf_A1N[-3], bf_A1N[-2]/100.0)



pl.figure()

#Plot the fits    
pl.semilogy(m, y_A0, 'r--', linewidth=2)
pl.semilogy(m, y_A1, 'b--', linewidth=2)
pl.semilogy(m, y_A0N, 'r-', linewidth=2, label=r'$k = 0$')
pl.semilogy(m, y_A1N, 'b-', linewidth=2, label=r'$k = 1$')
 
#Plot the ATLAS data
pl.errorbar(m_centres[np.where(data > 1e-5)],data[np.where(data > 1e-5)],
    xerr=20, yerr=[ calcError1(data[np.where(data > 1e-5)], 0.68), 
    calcError2(data[np.where(data > 1e-5)], 0.68)], 
    capsize=0, color='k', fmt='o', zorder=3)

pl.plot([], [], 'k-', linewidth=2.0, label='Free norm.')
pl.plot([], [], 'k--', linewidth=2.0, label='Fixed norm.')

pl.legend(fontsize=16.0, frameon=False)

if (include_signal):
    signalstr = 'Signal+Background fit\n'
else:
    signalstr = 'Background-only fit\n'

pl.text(250, 0.2, signalstr )
pl.text(250, 0.18, r'$n_\mathrm{bins} = 40$' )

pl.ylim(1e-1, 1e3)
pl.xlim(150,1750)

outfile = 'ATLASfits_BG-only.pdf'
if (include_signal):
    outfile = 'ATLASfits_BG+Sig.pdf'

pl.xlabel(r'$m_{\gamma \gamma}$ [GeV]', fontsize=18.0)
pl.ylabel('Events / 40 GeV')

pl.savefig(outfile, bbox_inches='tight')
pl.show()