module purge
module load python/3.5-gcc cmake/3.6.3 gcc eigen petsc slepc trilinos boost intel  intellib intelmpi
module list
source /home/cpbhanda/fenics/2017.1.0/share/dolfin/dolfin.conf

export CXX=mpiicpc 
export CC=mpiicc 
export FC=mpiifort
export F77=mpiifort

export CFLAGS="-O3 -mkl -Ofast -fast-transcendentals -ftz -fma -fp-model fast=2 -qopt-calloc  -fp-speculation=fast -fpic -ip -no-prec-div -xCORE-AVX-I  -O3"
export CXXFLAGS=$CFLAGS
export FCFLAGS=$CFLAGS

export INSTALL_DIR=/home/cpbhanda/fenics/2017.1.0
#export INSTALL_DIR=/project/cacds/apps/fenics/1.4.0-gcc
cmake .. \
-DCMAKE_CXX_COMPILER=/share/apps/intel/impi/5.1.1.109/intel64/bin/mpiicpc \
-DCMAKE_C_COMPILER=/share/apps/intel/impi/5.1.1.109/intel64/bin/mpiicc \
-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
-DDOLFIN_DIR=/home/cpbhanda/fenics/2017.1.0/share/dolfin/cmake \
-DXIAR=/share/apps/intel/compilers_and_libraries_2016/linux/bin/intel64/xiar

make
ls -l
