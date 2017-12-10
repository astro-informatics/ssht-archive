FROM ubuntu:16.04
MAINTAINER francois.lanusse@gmail.com

# Update Ubuntu Software repository
RUN apt-get update && \
    apt-get -y upgrade

# Install ssht dependencies
RUN apt-get install -y libfftw3-dev build-essential git cython && \
    apt-get install -y python-matplotlib python-numpy 

# Build and install
RUN git clone https://github.com/EiffL/ssht.git && \
    cd ssht && \
    make && \
    export FFTW=/usr && python setup.py install
