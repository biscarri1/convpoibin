# convpoibin

This repository contains code for implementing convolution methods for calculating the Poisson binomial distribution function, as described in https://www.sciencedirect.com/science/article/pii/S0167947318300082. The functions themselves are written in C, and then compiled and ported to R.

This code is currently **not** parallelized, although it is very clear how it could be, and significant gains in speed could be realized by paralellizing it. This is currently in development.


## Dependencies
The primary dependency for this code is the Fastest Fourier Transform in the West (fftw) which is available at http://www.fftw.org/download.html It is highly recommended to check out the installation tutorial here: http://www.fftw.org/fftw3_doc/Installation-and-Customization.html




## UNIX INSTALLTION
First you will need to install the fftw package. As an example, this can be done by typing the below into the command line. Note that the **--enable-shared** flag must be included.
```
wget http://fftw.org/fftw-3.3.5.tar.gz
tar -xzf fftw-3.3.5.tar.gz
cd fftw-3.3.5
./configure --enable-shared
make
sudo make install
```
Alternatively, it is possible to simple download the .tar.gz without using wget and continuing from the second step. It doesn't really matter how one gets the .tar.gz file, as long as they have it. fftw should be located in ```/usr/local/lib``` to check, you can try running ```ls /usr/local/lib | grep libfftw``` which should return some files with libfftw in the name if they are indeed in ```/usr/local/lib```

**Again, the "--enable-shared" flag seems to be necessary.**

Then, the convpoibin library can be installed using devtools and install.github in R like so:
```R
install.packages('devtools')
devtools::install_github('biscarri1/convpoibin')
```
