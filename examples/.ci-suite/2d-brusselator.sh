#!/bin/bash
echo " === starting 2d-brusselator example === "
# set symlinks
export workdir=examples/.ci-suite/data
export nproc=4
ln -sf examples/.ci-suite/2d-brusselator_eqns.idp eqns.idp
ln -sf examples/.ci-suite/2d-brusselator_settings.idp settings.idp
# build meshes
FreeFem++-mpi -v 0 examples/.ci-suite/2d-brusselator_mesh.edp -mo $workdir/square
# LINEAR STABILITY ANALYSIS OF CONSTANT SOLUTION
ff-mpirun -np $nproc modecompute.edp -v 0 -dir $workdir -mi square.msh -fo oe -A 2 -B 5.45 -Dx 0.008 -Dy 0.004 -1/L^2 0.25 -eps_target 0.3+2.1i -eps_nev 5 -eps_pos_gen_non_hermitian -sym 0,0
ff-mpirun -np $nproc modecompute.edp -v 0 -dir $workdir -mi square.msh -fo ee -A 2 -B 5.45 -Dx 0.008 -Dy 0.004 -1/L^2 0.25 -eps_target 0.3+2.1i -eps_nev 5 -eps_pos_gen_non_hermitian -sym 1,0
# COMPUTE HOPF POINTS
ff-mpirun -np $nproc hopfcompute.edp -v 0 -dir $workdir -fi oe.mode -fo brusselator10 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.edp -v 0 -dir $workdir -fi ee.mode -fo brusselator20 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.edp -v 0 -dir $workdir -fi oe_2.mode -fo brusselator30 -param 1/L^2 -snes_rtol 0
# CONTINUE HOPF POINT
ff-mpirun -np $nproc hopfcontinue.edp -v 0 -dir $workdir -fi brusselator10.hopf -fo brusselator10 -param B -param2 A -snes_rtol 0 -maxcount 10
# CONTINUE PERIODIC SOLUTION BRANCHES FROM HOPF POINTS
ff-mpirun -np $nproc porbcontinue.edp -v 0 -dir $workdir -fi brusselator30.hopf -fo brusselator30 -param 1/L^2 -maxcount 3 -h0 -1 -Nh 6
ff-mpirun -np $nproc porbcontinue.edp -v 0 -dir $workdir -fi brusselator30_3.porb -fo brusselator30 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.25 -dmax 1e10 -mono 0 -contorder 0 -count 3
# DEFLATION OF PERIODIC SOLUTION BRANCH
ff-mpirun -np $nproc porbcompute.edp -v 0 -dir $workdir -fi brusselator30_4.porb -fo brusselator_branch-0 -1/L^2 0.25
ff-mpirun -np $nproc porbdeflate.edp -v 0 -dir $workdir -fi brusselator_branch-0.porb -fi2 brusselator_branch-0.porb -fo brusselator_branch -ndeflate 2 -snes_rtol 0 -noise 0.03 -snes_linesearch_type l2 -snes_linesearch_damping 0.8
# FLOQUET ANALYSIS OF PERIODIC ORBITS 
ff-mpirun -np $nproc floqcompute.edp -v 0 -dir $workdir -fi brusselator_branch-0.porb -so brusselator_branch-0 -eps_target 0.3+2.1i -sym 0,0 -eps_pos_gen_non_hermitian -eps_nev 10 -blocks 1
ff-mpirun -np $nproc floqcompute.edp -v 0 -dir $workdir -fi brusselator_branch-0.porb -so brusselator_branch-0 -eps_target 0.3+2.1i -sym 1,0 -eps_pos_gen_non_hermitian -eps_nev 10 -blocks 2
ff-mpirun -np $nproc floqcompute.edp -v 0 -dir $workdir -fi brusselator_branch-0.porb -so brusselator_branch-0 -eps_target 0.3+2.1i -sym 0,1 -eps_pos_gen_non_hermitian -eps_nev 10 -blocks 3
ff-mpirun -np $nproc floqcompute.edp -v 0 -dir $workdir -fi brusselator_branch-0.porb -so brusselator_branch-0 -eps_target 0.3+2.1i -sym 1,1 -eps_pos_gen_non_hermitian -eps_nev 10 -blocks 4
