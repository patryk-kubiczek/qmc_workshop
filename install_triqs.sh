# Installing software and libraries required by TRIQS
# Make sure you understand what you are doing! 
sudo apt-get install libboost-all-dev cmake git g++ libgfortran3 gfortran openmpi-bin openmpi-common \
     	openmpi-doc libopenmpi-dev libblas-dev liblapack-dev libfftw3-dev libgmp-dev \
     	hdf5-tools libhdf5-serial-dev python-dev \
# You can also install python modules locally using e.g. pip.
sudo apt-get install python-h5py python-numpy python-scipy python-jinja2 python-virtualenv \
	python-matplotlib python-tornado python-zmq python-mpi4py python-mako

# Downloading and installing TRIQS and applications based on TRIQS
# Change to directory where you want to install triqs and...
mkdir -p triqs && cd triqs
git clone https://github.com/TRIQS/triqs.git src
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../install ../src
make -j4
make test
make install
cd .. # go to triqs
mkdir -p cthyb && cd cthyb
git clone https://github.com/TRIQS/cthyb.git src
mkdir -p build && cd build
cmake -DTRIQS_PATH=../../install ../src
make -j4
make test
make install
cd .. & cd .. # go to triqs
mkdir -p som && cd som
git clone https://github.com/krivenko/som src
mkdir -p build && cd build
cmake -DTRIQS_PATH=../../install ../src
make -j4
make test
make install
cd .. & cd .. # go to triqs

# Alternatively, if you connect to a machine working within the physnet, just load a TRIQS module
# module purge
# module load triqs/1.4


