# ATLASfits

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45646.svg)](http://dx.doi.org/10.5281/zenodo.45646)

ATLASfits is a collection of python code (and a couple of data files) which allows you to perform fits to the (binned, digitised) ATLAS diphoton invariant mass spectrum, as described at http://cds.cern.ch/record/2114853. This code should allow you to reproduce the results of http://arxiv.org/abs/1601.07330. Please cite that paper if you make use of this code in your own work. 

Please send any comments or corrections to: bradkav@gmail.com.

Note: Requires emcee, available at http://dan.iel.fm/emcee/current/.

To Run: Simply run 'python ATLASfits.py' to calculate the best fit points and likelihoods, followed by 'python PlotFits.py' to generate the plots.

Description of files:

- ATLASfits.py: Calculate the best fit points and likelihoods for background-only, narrow signal and wide signal hypotheses (for 6 different background assumptions). Set 'saveResults=1' near the start of the code to save to file the best fit parameter values for plotting later.
- PlotFits.py: Plot the best fit parameter values as calculated by ATLASfits.py. Set 'include_signal=0' near the start of the code to plot the background-only results. Set 'include_signal=1' for the wide resonance signal + background fits.
- AddScatter.py: Run to add random noise to the first 10 ATLAS data bins (as described in arXiv:1601.07330). AddScatter.py reads in 'ATLASdata.txt' and saves the noisy data to 'ATLASdata1.txt'. Note that ATLASfits.py and PlotFits.py perform calculations using 'ATLASdata1.txt'.
- ATLASfits_utils.py: Module containing likelihood and other auxiliary functions for fitting and plotting.
- ATLASdata.txt: Digitised ATLAS diphoton data.
- ATLASdata1.txt: A copy of ATLASdata.txt which is read in by other files. Running AddScatter.py will add digitisation noise to ATLASdata1.txt
- 1601.07330.pdf: The paper in which the results of these calculations are reported.
