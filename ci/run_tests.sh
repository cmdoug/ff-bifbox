#!/usr/bin/env bash
set -euo pipefail
echo "=== CI: run_tests.sh ==="

# Prepare results directory
mkdir -p test-results

# Sipp and Lebedev example
export workdir=test-results
export nproc=4
ln -sf examples/sipp_lebedev_2007/eqns_sipp_lebedev_2007.idp eqns.idp
ln -sf examples/sipp_lebedev_2007/settings_sipp_lebedev_2007.idp settings.idp\

#create mesh
FreeFem++-mpi -v 0 examples/sipp_lebedev_2007/cylinder.edp -mo $workdir/cylinder
FreeFem++-mpi -v 0 examples/sipp_lebedev_2007/cavity.edp -mo $workdir/cavity

#base flow
ff-mpirun -np $nproc basecompute.edp -v 0 -dir $workdir -mi cylinder.msh -fo cylinder -1/Re 0.1
ff-mpirun -np $nproc basecompute.edp -v 0 -dir $workdir -mi cavity.msh -fo cavity -1/Re 0.1

echo "=== All CI tests completed successfully ==="
