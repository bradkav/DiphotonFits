#--AddScatter.py - Version 1 - 04/02/2016
#--Author: Bradley J Kavanagh
#--Summary: Read in 'ATLASdata.txt', add scatter to first
#--ten bins and save as 'ATLASdata1.txt'
#--Please report any problems to: bradley.kavanagh@lpthe.jussieu.fr

import numpy as np

print "----Adding scatter to the ATLAS data to simulate digitisation error---"

#Read in data
xvals = np.loadtxt("ATLASdata.txt")[:,0]
data = np.loadtxt("ATLASdata.txt")[:,1]

#Print initial data
print "ATLAS data (initial):"
print data

#Add in random scatter
for i in range(10):
    #data[i] = np.random.normal(data[i], 0.05*data[i])
    data[i] = round(np.random.uniform(data[i]*0.97, data[i]*1.03))

#Print final data
print "ATLAS data (with digitisation noise):"
print data

#Save to file (DiphotonFits.py reads in 'ATLASdata1.txt' *not* 'ATLASdata.txt')
np.savetxt("ATLASdata1.txt", np.c_[xvals, data])



