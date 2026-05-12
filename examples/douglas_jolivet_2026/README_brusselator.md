# Examples: Douglas & Jolivet (202x)
This file shows an example `ff-bifbox` workflow for reproducing the results of the study:
```tex
@article{douglas_jolivet_2026,
  title={ff-bifbox: A scalable, open-source toolbox for bifurcation analysis of nonlinear PDEs},
  author={Douglas, Christopher M. and Jolivet, Pierre R.},
  year={2026},
}
```
The commands below illustrate how to perform a bifurcation analysis of the Brusselator in 3-D using `ff-bifbox`.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/douglas_jolivet_2026/brusselator
export nproc=4
```

# Example 1: 3-D Brusselator

1. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/douglas_jolivet_2026/eqns_douglas_jolivet_2026_brusselator.idp eqns.idp
ln -sf examples/douglas_jolivet_2026/settings_douglas_jolivet_2026_brusselator.idp settings.idp
```

2. Build initial meshes
```sh
FreeFem++-mpi -v 0 examples/douglas_jolivet_2026/mesh_brusselator.md -mo $workdir/cube
```

3. Compute leading modes over the base state at $L=2$.
```sh
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi cube_eighth.mesh -fo sss -A 2 -B 5.45 -Dx 0.008 -Dy 0.004 -1/L^2 0.25 -eps_target 0.3+2.1i -eps_nev 6 -eps_pos_gen_non_hermitian -sym 0,0,0
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi cube_eighth.mesh -fo nss -A 2 -B 5.45 -Dx 0.008 -Dy 0.004 -1/L^2 0.25 -eps_target 0.3+2.1i -eps_nev 6 -eps_pos_gen_non_hermitian -sym 1,0,0
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi cube_eighth.mesh -fo sns -A 2 -B 5.45 -Dx 0.008 -Dy 0.004 -1/L^2 0.25 -eps_target 0.3+2.1i -eps_nev 6 -eps_pos_gen_non_hermitian -sym 0,1,0
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi cube_eighth.mesh -fo nns -A 2 -B 5.45 -Dx 0.008 -Dy 0.004 -1/L^2 0.25 -eps_target 0.3+2.1i -eps_nev 6 -eps_pos_gen_non_hermitian -sym 1,1,0
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi cube_eighth.mesh -fo snn -A 2 -B 5.45 -Dx 0.008 -Dy 0.004 -1/L^2 0.25 -eps_target 0.3+2.1i -eps_nev 6 -eps_pos_gen_non_hermitian -sym 0,1,1
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi cube_eighth.mesh -fo nnn -A 2 -B 5.45 -Dx 0.008 -Dy 0.004 -1/L^2 0.25 -eps_target 0.3+2.1i -eps_nev 6 -eps_pos_gen_non_hermitian -sym 1,1,1
```

4. Compute Hopf points
- 1-D solutions:
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sss.mode -fo brusselator100 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi nss.mode -fo brusselator200 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sss_3.mode -fo brusselator300 -param 1/L^2 -snes_rtol 0
```
- 2-D solutions:
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sns.mode -fo brusselator110 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sss_1.mode -fo brusselator120 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sns_3.mode -fo brusselator130 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi nns.mode -fo brusselator210 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi nss_1.mode -fo brusselator220 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi nns_2.mode -fo brusselator230 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sns_2.mode -fo brusselator310 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sss_5.mode -fo brusselator320 -param 1/L^2 -snes_rtol 0
```
- 3-D solutions:
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi snn.mode -fo brusselator111 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi nnn.mode -fo brusselator211 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sns_1.mode -fo brusselator112 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi nns_1.mode -fo brusselator212 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sss_4.mode -fo brusselator122 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi snn_3.mode -fo brusselator311 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi snn_2.mode -fo brusselator131 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi nss_3.mode -fo brusselator222 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sns_4.mode -fo brusselator312 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi sns_5.mode -fo brusselator132 -param 1/L^2 -snes_rtol 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi nnn_2.mode -fo brusselator231 -param 1/L^2 -snes_rtol 0
```

5. Continue Hopf bifurcations along branches of periodic orbits up to $L=2$
- 1-D solutions:
```sh
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator100.hopf -fo brusselator100 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator100_2.porb -fo brusselator100 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.25 -dmax 1e10 -mono 0 -contorder 0 -count 2
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator200.hopf -fo brusselator200 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator200_2.porb -fo brusselator200 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.1 -dmax 1e10 -mono 0 -contorder 0 -count 2
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator300.hopf -fo brusselator300 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator300_2.porb -fo brusselator300 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.05 -dmax 1e10 -mono 0 -contorder 0 -count 2
```
- 2-D solutions
```sh
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator110.hopf -fo brusselator110 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator110_2.porb -fo brusselator110 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.2 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator120.hopf -fo brusselator120 -param 1/L^2 -maxcount 2 -h0 -2 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator120_2.porb -fo brusselator120 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.05 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator130.hopf -fo brusselator130 -param 1/L^2 -maxcount 2 -h0 -2 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator130_2.porb -fo brusselator130 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.01 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator210.hopf -fo brusselator210 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator210_2.porb -fo brusselator210 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.1 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator220.hopf -fo brusselator220 -param 1/L^2 -maxcount 2 -h0 -2 -Nh 6 -amax 90 -kmax 1e10 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator220_2.porb -fo brusselator220 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.05 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator230.hopf -fo brusselator230 -param 1/L^2 -maxcount 2 -h0 -3 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator230_2.porb -fo brusselator230 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.01 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator310.hopf -fo brusselator310 -param 1/L^2 -maxcount 2 -h0 -2 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator310_2.porb -fo brusselator310 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.05 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator320.hopf -fo brusselator320 -param 1/L^2 -maxcount 2 -h0 -1 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator320_2.porb -fo brusselator320 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.0025 -dmax 1e10 -mono 0 -contorder 0 -count 2
```
- 3-D solutions
```sh
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator111.hopf -fo brusselator111 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator111_2.porb -fo brusselator111 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.1 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator211.hopf -fo brusselator211 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator211_2.porb -fo brusselator211 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.1 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator112.hopf -fo brusselator112 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator112_2.porb -fo brusselator112 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.05 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator212.hopf -fo brusselator212 -param 1/L^2 -maxcount 2 -h0 -5 -Nh 6 -amax 90
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator212_2.porb -fo brusselator212 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.05 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator122.hopf -fo brusselator122 -param 1/L^2 -maxcount 1 -h0 -5 -Nh 6
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator122_1.porb -fo brusselator122 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.02 -dmax 1e10 -mono 0 -contorder 0 -count 1

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator311.hopf -fo brusselator311 -param 1/L^2 -maxcount 2 -h0 -2 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator311_2.porb -fo brusselator311 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.02 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator131.hopf -fo brusselator131 -param 1/L^2 -maxcount 2 -h0 -2 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator131_2.porb -fo brusselator131 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.01 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator222.hopf -fo brusselator222 -param 1/L^2 -maxcount 2 -h0 -1 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator222_2.porb -fo brusselator222 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.005 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator312.hopf -fo brusselator312 -param 1/L^2 -maxcount 2 -h0 -1 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator312_2.porb -fo brusselator312 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.00125 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator132.hopf -fo brusselator132 -param 1/L^2 -maxcount 2 -h0 -1 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator132_2.porb -fo brusselator132 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.0025 -dmax 1e10 -mono 0 -contorder 0 -count 2

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator231.hopf -fo brusselator231 -param 1/L^2 -maxcount 2 -h0 -1 -Nh 6 -kmax 1e10 -amax 90 -mono 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator231_2.porb -fo brusselator231 -param 1/L^2 -paramtarget 0.25 -maxcount -1 -h0 -0.0025 -dmax 1e10 -mono 0 -contorder 0 -count 2
```

6. (OPTIONAL) Floquet analysis calculations along the 1-D solution branches
```sh
cd "$workdir" && declare -a porblist=(brusselator100_*[0-9].porb) && cd -
for porbfile in "${porblist[@]}"; do
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100 -eps_target 0.3+2.1i -sym 0,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100 -eps_target 0.3+2.1i -sym 1,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100 -eps_target 0.3+2.1i -sym 0,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100 -eps_target 0.3+2.1i -sym 1,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100 -eps_target 0.3+2.1i -sym 0,1,1 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100 -eps_target 0.3+2.1i -sym 1,1,1 -eps_pos_gen_non_hermitian -eps_nev 10
done

cd "$workdir" && declare -a porblist=(brusselator200_*[0-9].porb) && cd -
for porbfile in "${porblist[@]}"; do
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200 -eps_target 0.3+2.1i -sym 0,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200 -eps_target 0.3+2.1i -sym 1,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200 -eps_target 0.3+2.1i -sym 0,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200 -eps_target 0.3+2.1i -sym 1,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200 -eps_target 0.3+2.1i -sym 0,1,1 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200 -eps_target 0.3+2.1i -sym 1,1,1 -eps_pos_gen_non_hermitian -eps_nev 10
done

cd "$workdir" && declare -a porblist=(brusselator300_*[0-9].porb) && cd -
for porbfile in "${porblist[@]}"; do
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300 -eps_target 0.3+2.1i -sym 0,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300 -eps_target 0.3+2.1i -sym 1,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300 -eps_target 0.3+2.1i -sym 0,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300 -eps_target 0.3+2.1i -sym 1,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300 -eps_target 0.3+2.1i -sym 0,1,1 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300 -eps_target 0.3+2.1i -sym 1,1,1 -eps_pos_gen_non_hermitian -eps_nev 10
done
```

6. Compute other 1-D branches using deflation and continue them
```sh
ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi brusselator100_10.porb -fo brusselator100_nss -1/L^2 0.65
ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi brusselator100_nss.porb -fo brusselator100_nss -eps_target 0.01+1.81i -sym 1,0,0 -eps_pos_gen_non_hermitian -eps_nev 5
ff-mpirun examples/douglas_jolivet_2026/extend_brusselator.md -v 0 -dir $workdir -fi brusselator100_nss.floq -fo brusselator100_nss -amp 0.1
ff-mpirun -np $nproc porbdeflate.md -v 0 -dir $workdir -mi cube_Qx.mesh -fi brusselator100_nssguess.porb -fi2 brusselator100_nssstart.porb -fo brusselator100_nss -ndeflate 1 -snes_rtol 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator100_nss-1.porb -fo brusselator100_nss -param 1/L^2 -h0 1 -paramtarget 0.25 -maxcount -1

ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi brusselator100_12.porb -fo brusselator100_nss2 -1/L^2 0.28
ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi brusselator100_nss2.porb -fo brusselator100_nss2 -eps_target 0.01+1.84i -sym 1,0,0 -eps_pos_gen_non_hermitian -eps_nev 5
ff-mpirun examples/douglas_jolivet_2026/extend_brusselator.md -v 0 -dir $workdir -fi brusselator100_nss2.floq -fo brusselator100_nss2 -amp 0.1
ff-mpirun -np $nproc porbdeflate.md -v 0 -dir $workdir -mi cube_Qx.mesh -fi brusselator100_nss2guess.porb -fi2 brusselator100_nss2start.porb -fo brusselator100_nss2 -ndeflate 1 -snes_rtol 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator100_nss2-1.porb -fo brusselator100_nss2 -param 1/L^2 -h0 -1 -paramtarget 0.25 -maxcount -1

ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi brusselator200_5.porb -fo brusselator200_sss -1/L^2 0.57
ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi brusselator200_sss.porb -fo brusselator200_sss -eps_target 0+1.98i -sym 0,0,0 -eps_pos_gen_non_hermitian -eps_nev 5
ff-mpirun examples/douglas_jolivet_2026/extend_brusselator.md -v 0 -dir $workdir -fi brusselator200_sss_2.floq -fo brusselator200_sss -amp 0.1
ff-mpirun -np $nproc porbdeflate.md -v 0 -dir $workdir -mi cube_Qx.mesh -fi brusselator200_sssguess.porb -fi2 brusselator200_sssstart.porb -fo brusselator200_sss -ndeflate 1 -snes_rtol 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator200_sss-1.porb -fo brusselator200_sss -param 1/L^2 -h0 -0.005 -paramtarget 0.25 -maxcount -1 -contorder 0 -mono 0 -dmax 1e10

ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi brusselator300_4.porb -fo brusselator300_nss -1/L^2 0.285
ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi brusselator300_nss.porb -fo brusselator300_nss -eps_target 0+2.01i -sym 1,0,0 -eps_pos_gen_non_hermitian -eps_nev 5
ff-mpirun examples/douglas_jolivet_2026/extend_brusselator.md -v 0 -dir $workdir -fi brusselator300_nss_3.floq -fo brusselator300_nss -amp 0.1
ff-mpirun -np $nproc porbdeflate.md -v 0 -dir $workdir -mi cube_Qx.mesh -fi brusselator300_nssguess.porb -fi2 brusselator300_nssstart.porb -fo brusselator300_nss -ndeflate 1 -snes_rtol 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi brusselator300_nss-1.porb -fo brusselator300_nss -param 1/L^2 -h0 -0.005 -paramtarget 0.25 -maxcount -1 -contorder 0 -mono 0 -dmax 1e10
```

6. (OPTIONAL) Floquet analysis calculations along the 1-D solution branches
```sh

cd "$workdir" && declare -a porblist=(brusselator100_nss_*[0-9].porb) && cd -
for porbfile in "${porblist[@]}"; do
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100_nss -eps_target 0.3+1.9i -sym 0,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100_nss -eps_target 0.3+1.9i -sym 0,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100_nss -eps_target 0.3+1.9i -sym 0,1,1 -eps_pos_gen_non_hermitian -eps_nev 10
done

cd "$workdir" && declare -a porblist=(brusselator100_nss2_*[0-9].porb) && cd -
for porbfile in "${porblist[@]}"; do
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100_nss2 -eps_target 0.3+1.9i -sym 0,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100_nss2 -eps_target 0.3+1.9i -sym 0,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator100_nss2 -eps_target 0.3+1.9i -sym 0,1,1 -eps_pos_gen_non_hermitian -eps_nev 10
done

cd "$workdir" && declare -a porblist=(brusselator200_sss_*[0-9].porb) && cd -
for porbfile in "${porblist[@]}"; do
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200_sss -eps_target 0.3+1.85i -sym 0,0,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200_sss -eps_target 0.3+1.85i -sym 0,1,0 -eps_pos_gen_non_hermitian -eps_nev 10
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator200_sss -eps_target 0.3+1.85i -sym 0,1,1 -eps_pos_gen_non_hermitian -eps_nev 10 
done

cd "$workdir" && declare -a porblist=(brusselator300_nss_*[0-9].porb) && cd -
for porbfile in "${porblist[@]}"; do
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300_nss -eps_target 0.3+2i -sym 0,0,0 -eps_pos_gen_non_hermitian -eps_nev 16
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300_nss -eps_target 0.3+2i -sym 0,1,0 -eps_pos_gen_non_hermitian -eps_nev 16
  ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi $porbfile -so brusselator300_nss -eps_target 0.3+2i -sym 0,1,1 -eps_pos_gen_non_hermitian -eps_nev 16 
done
```