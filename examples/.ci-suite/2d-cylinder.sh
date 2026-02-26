#!/bin/bash
echo " === starting 2d-cylinder example === "
# set symlinks
export workdir=examples/.ci-suite/data
export nproc=4
ln -sf examples/.ci-suite/2d-cylinder_eqns.idp eqns.idp
ln -sf examples/.ci-suite/2d-cylinder_settings.idp settings.idp
# build meshes
FreeFem++-mpi -v 0 examples/.ci-suite/2d-cylinder_mesh.edp -mo $workdir/cylinder_0
# compute and continue base state
ff-mpirun -np $nproc basecompute.edp -v 0 -dir $workdir -mi cylinder_0.msh -1/Re 1 -fo cylinder_0
ff-mpirun -np $nproc basecontinue.edp -v 0 -dir $workdir -fi cylinder_0.base -fo cylinder -h0 -1 -param 1/Re -paramtarget 0.01 -maxcount -1 -mo cylinder -scount 1
ff-mpirun -np $nproc basecompute.edp -v 0 -dir $workdir -fi cylinder_10.base -1/Re 0.021 -fo cylinder50
# stability analysis
ff-mpirun -np $nproc modecompute.edp -v 0 -dir $workdir -fi cylinder50.base -fo cylinder50 -eps_target 0.1+0.8i -sym 1 -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc hopfcompute.edp -v 0 -dir $workdir -fi cylinder50.mode -fo cylinder -param 1/Re -nf 0
ff-mpirun -np $nproc hopfcompute.edp -v 0 -dir $workdir -fi cylinder.hopf -fo cylinderadapt -mo cylinderhopf -adaptto bo -param 1/Re -thetamax 5 -wnl 1 
ff-mpirun -np $nproc porbcontinue.edp -v 0 -dir $workdir -fi cylinder.hopf -fo cylinderNh1 -Nh 1 -param 1/Re -h0 -1 -scount 1 -maxcount 10
ff-mpirun -np $nproc porbcompute.edp -v 0 -dir $workdir -fi cylinderNh1_10.porb -fo cylinder50Nh2 -Nh 2 -1/Re 0.02 -blocks 2
