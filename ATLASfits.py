#--ATLASfits.py - Version 1 - 04/02/2016
#--Author: Bradley J Kavanagh
#--Summary: Code for fitting the ATLAS diphoton data
#--and calculating the significance of the 750 GeV excess
#--Note: Requires emcee (http://dan.iel.fm/emcee/current/)
#--Please report any problems to: bradley.kavanagh@lpthe.jussieu.fr



print "----Likelihood fits to ATLAS diphoton data---"

import numpy as np
from ATLASfits_utils import getBestFit

#----Options----
#Print best-fit points to file, for plotting
#Note: This will generate 18 files...
saveResults = 1 

#----Main procedure-----

BG_ID = ['k = 0 (fixed norm.)',
         'k = 1 (fixed norm.)',
         'k = 2 (fixed norm.)',
         'k = 0 (free norm.)',
         'k = 1 (free norm.)',
         'k = 2 (free norm.)',]

SIGNAL_ID = ['Background-only',
             'Signal+BG (NWA)',
             'Signal+BG (Free-width)']

#Loop over possible background parametrisations
for i in range(6):
    print "---------------------"
    print "Background function:" ,BG_ID[i]   

    #Background-only fit
    like_BG, bf_BG = getBestFit(i, 0)
    print "    Background-only"
    print "        lnL:", '{:.2f}'.format(like_BG)
    if (saveResults):
        np.savetxt('Fits_BG=' + str(i) + '_BG-only.txt', bf_BG)
    
    #Narrow width fit
    like_NWA, bf_NWA = getBestFit(i, 1)    
    print "    Signal+BG (NWA)"
    print "        lnL:", '{:.2f}'.format(like_NWA)
    if (saveResults):
        np.savetxt('Fits_BG=' + str(i)+ '_NWA.txt', bf_NWA)

    #Free width fit
    like_wide, bf_wide = getBestFit(i, 2)    
    print "    Signal+BG (Free width)"
    print "        lnL:", '{:.2f}'.format(like_wide)
    if (saveResults):
        np.savetxt('Fits_BG=' + str(i) + '_wide.txt', bf_wide)

    #Calculate significance
    sig_NWA = np.sqrt(2*(like_NWA - like_BG))
    sig_wide = np.sqrt(2*(like_wide - like_BG))

    print " "
    print "    Significance (NWA):", '{:.2f}'.format(sig_NWA)
    print "    Significance (wide):", '{:.2f}'.format(sig_wide)

