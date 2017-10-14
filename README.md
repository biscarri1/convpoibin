# convpoibin

This repository contains code for implementing convolution methods for calculating the Poisson binomial distribution function. The functions themselves are written in C, and then compiled and ported to R.


## Dependencies
The primary dependency for this code is the Faster Fourier Transform in the West (fftw) which is available at http://www.fftw.org/download.html It is highly recommended to check out the installation tutorial here: http://www.fftw.org/fftw3_doc/Installation-and-Customization.html

For UNIX systems, it is highly recommended to add an "--enable-shared" flag to the ./configure command. 

More in depth directions are in development.
