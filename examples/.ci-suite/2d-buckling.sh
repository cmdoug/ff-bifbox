#!/bin/bash
echo " === starting 2d-buckling example === "
# set symlinks
export workdir=examples/.ci-suite/data
export nproc=4
ln -sf examples/.ci-suite/2d-buckling_eqns.idp eqns.idp
ln -sf examples/.ci-suite/2d-buckling_settings.idp settings.idp
# build meshes
FreeFem++-mpi -v 0 examples/.ci-suite/2d-buckling_mesh.edp -mo $workdir/beam
# continue base state
ff-mpirun -np $nproc basecontinue.edp -v 0 -dir $workdir -mi beam.msh -fo buckling -param P -nu 0.3 -snes_rtol 0 -h0 0.1 -amax 15 -tgv -2 -paramtarget 0.1 -maxcount -1 -dmax 0.1 -kmax 0.25
ff-mpirun -np $nproc basecompute.edp -v 0 -dir $workdir -fi buckling_1.base -fo buckling_0 -P 0 -tgv -2
# compute and continue fold points
cd "$workdir" && declare -a foldguesslist=(buckling_*specialpt.base) && cd -
for guess in "${foldguesslist[@]}"; do
export out=$(echo "$guess" | awk -F'specialpt' '{print $1}')
ff-mpirun -np $nproc foldcompute.edp -v 0 -dir $workdir -fi $guess -fo $out -param P -tgv -2 -snes_atol 1e-10
ff-mpirun -np $nproc foldcontinue.edp -v 0 -dir $workdir -fi "$out".fold -fo $out -param P -param2 nu -tgv -2 -snes_atol 1e-10 -maxcount -1 -param2target 0.49
done
# stability analysis
cd "$workdir" && declare -a baselist=(buckling_*[0-9].base) && cd -
for base in "${baselist[@]}"; do
ff-mpirun -np $nproc modecompute.edp -v 0 -dir $workdir -fi $base -so buckling_sym -eps_target 0.1+0i -eps_nev 3 -eps_gen_hermitian -sym 0
ff-mpirun -np $nproc modecompute.edp -v 0 -dir $workdir -fi $base -so buckling_asym -eps_target 0.1+0i -eps_nev 3 -eps_gen_hermitian -sym 1
done
# pitchfork bifurcations
ff-mpirun -np $nproc basecompute.edp -v 0 -dir $workdir -fi buckling_0.base -fo buckling_A -P 0.05 -tgv -2
ff-mpirun -np $nproc modecompute.edp -v 0 -dir $workdir -fi buckling_A.base -fo buckling_A -eps_target 0.1+0i -eps_nev 3 -eps_gen_hermitian -sym 1
ff-mpirun -np $nproc hopfcompute.edp -v 0 -dir $workdir -fi buckling_A.mode -fo buckling_A -param P -tgv -2 -snes_atol 1e-10 -zero 1
ff-mpirun -np $nproc hopfcontinue.edp -v 0 -dir $workdir -fi buckling_A.hopf -fo buckling_A.hopf -param P -param2 nu -tgv -2 -snes_atol 1e-10 -maxcount -1 -param2target 0.49 -zero 1
ff-mpirun -np $nproc modecompute.edp -v 0 -dir $workdir -fi "$out".fold -fo buckling_B -eps_target 0.1+0i -eps_nev 3 -eps_gen_hermitian -sym 1
ff-mpirun -np $nproc hopfcompute.edp -v 0 -dir $workdir -fi buckling_B.mode -fo buckling_B -param P -tgv -2 -snes_atol 1e-10 -zero 1
ff-mpirun -np $nproc hopfcontinue.edp -v 0 -dir $workdir -fi buckling_B.hopf -fo buckling_B.hopf -param P -param2 nu -tgv -2 -snes_atol 1e-10 -maxcount -1 -param2target 0.49 -zero 1
# deflation
ff-mpirun -np $nproc basedeflate.edp -v 0 -dir $workdir -fi buckling_A.base -fo buckling_branch -ndeflate 3 -tgv -2 -snes_linesearch_linesearch l2 -snes_rtol 0 -defp 1 -P 0.04 -snes_divergence_tolerance 1e15 -snes_max_it 1000
