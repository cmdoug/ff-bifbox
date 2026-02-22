#!/usr/bin/env bash
set -euo pipefail
echo "=== CI: install_deps.sh ==="

# Basic system deps
sudo apt update
sudo apt install gcc g++ gfortran m4 patch git wget cmake libhdf5-dev flex bison autoconf automake autotools-dev

# Clone and build PETSc (cmdoug fork)
git clone https://gitlab.com/cmdoug/petsc
cd petsc
./configure PETSC_ARCH=real --with-debugging=0 --with-scalar-type=real --download-mpich --download-mumps --download-parmetis --download-metis --download-hypre --download-superlu --download-slepc --download-hpddm --download-ptscotch --download-suitesparse --download-scalapack --download-tetgen --download-mmg --download-parmmg --download-gsl=https://mirrors.kernel.org/gnu/gsl/gsl-2.8.tar.gz --download-bison=https://mirrors.kernel.org/gnu/bison/bison-3.8.2.tar.gz --download-slepc-configure-arguments="--download-arpack"
make PETSC_ARCH=real
./configure PETSC_ARCH=complex --with-debugging=0 --with-scalar-type=complex --with-mpich-dir=real --with-mumps-dir=real --with-parmetis-dir=real --with-metis-dir=real --with-superlu-dir=real --download-slepc --download-hpddm --download-htool --with-ptscotch-dir=real --with-suitesparse-dir=real --with-scalapack-dir=real --with-tetgen-dir=real --with-gsl-dir=real --with-bison-dir=real
make PETSC_ARCH=complex
cd -

# Clone and build FreeFEM (cmdoug fork)
git clone https://github.com/cmdoug/FreeFem-sources
cd FreeFem-sources
git checkout develop
autoreconf -i
./configure --enable-download --enable-optim --with-petsc="${HOME}/petsc/real" --with-petsc_complex="${HOME}/petsc/complex" --prefix="${HOME}/FreeFem-sources"
./3rdparty/getall -a
./reconfigure
make -j 4
make install
cd -

# Add FreeFEM to $PATH and export necessary environment variables
export PATH=${PATH}:$HOME/FreeFem-sources/src/mpi:$HOME/FreeFem-sources/src/nw
export FF_INCLUDEPATH="$HOME/FreeFem-sources/idp"
export FF_LOADPATH="$HOME/FreeFem-sources/plugin/mpi;;$HOME/FreeFem-sources/plugin/seq"

echo "=== Done: install_deps.sh ==="
