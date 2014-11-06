# 
# Dockerfile for 3BEM.
# 

# Pull base image
FROM ubuntu
MAINTAINER T. Ben Thompson <t.ben.thompson@gmail.com>

# Install necessary apt packages (keep in alphabetic order)
RUN apt-get update && apt-get install -y \
    alien\
    cmake\
    freeglut3\ 
    freeglut3-dev\
    git\
    libglew1.5\
    libglew1.5-dev\
    libglu1-mesa\
    libglu1-mesa-dev\
    libgl1-mesa-glx\
    libgl1-mesa-dev\
    libhdf5-mpi-dev\
    libnuma1\
    opencl-headers\
    petsc-dev\
    python-dev\ 
    python-distribute\
    python-pip\
    python-h5py\
    subversion\
    wget


# Install the build tool -- fabricate.py
RUN pip install fabricate

# Setup ssh keys to allow cloning a private github repo.
RUN mkdir -p /root/.ssh
ADD repo-key /root/.ssh/id_rsa
RUN chmod 700 /root/.ssh/id_rsa
RUN echo "Host github.com\n\tStrictHostKeyChecking no\n" >> /root/.ssh/config

# Clone the repository.
RUN git clone git@github.com:tbenthompson/active.git /home/active

# Download the libraries
WORKDIR /home/active
RUN mkdir lib
WORKDIR lib

# Autocheck
RUN git clone git@github.com:tbenthompson/autocheck.git

# Install OpenCL
RUN mkdir /home/active/opencl
WORKDIR /home/active/opencl
RUN wget http://registrationcenter.intel.com/irc_nas/4181/opencl_runtime_14.2_x64_4.5.0.8.tgz 
RUN tar -xvf opencl_runtime_14.2_x64_4.5.0.8.tgz 
WORKDIR pset_opencl_runtime_14.1_x64_4.5.0.8
WORKDIR rpm
RUN alien *.rpm
RUN dpkg -i *.deb
RUN mkdir -p /etc/OpenCL/vendors
RUN ln -s /opt/intel/opencl-1.2-4.5.0.8/lib64/libOpenCL.so /usr/lib
RUN ln -s /opt/intel/opencl-1.2-4.5.0.8/etc/intel64.icd /etc/OpenCL/vendors/

# C++ Actor-framework
WORKDIR /home/active/lib
RUN git clone https://github.com/actor-framework/actor-framework.git
WORKDIR actor-framework
# # Grab the OpenCL CAF submodule
RUN git submodule init
RUN git submodule update
# # Actually build
RUN ./configure --no-cash
RUN make -j12

# UnitTest++
WORKDIR /home/active/lib
RUN svn checkout http://svn.code.sf.net/p/unittest-cpp/code/UnitTest++ unittest-cpp
WORKDIR unittest-cpp
RUN make all

# Build the code
WORKDIR /home/active
ENV PETSC_DIR /etc/alternatives/petsc
ENV PETSC_ARCH linux-gnu-c-opt
RUN python build.py
RUN ./all_tests

# Some python packages that are useful 
# (h5py -- depends on cython, mpi4py)
RUN git clone https://github.com/h5py/h5py.git /h5py
WORKDIR /h5py
RUN pip install cython mpi4py
RUN python configure.py --mpi
RUN python setup.py build
RUN python setup.py install
