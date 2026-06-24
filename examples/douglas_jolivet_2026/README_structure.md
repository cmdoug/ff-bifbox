# Examples: Douglas & Jolivet (2026)
This file shows an example `ff-bifbox` workflow for reproducing the results of the study:
```bibtex
@article{douglas_jolivet_2026,
author = {Douglas, Christopher M. and Jolivet, Pierre},
title = {{f}f-bifbox: {A} scalable, open-source toolbox for bifurcation analysis of nonlinear {PDE}s},
journal = {Computer Physics Communications},
volume = {326},
pages = {110221},
year = {2026},
Publisher = {Elsevier},
issn = {0010-4655},
doi = {10.1016/j.cpc.2026.110221},
url = {https://doi.org/10.1016/j.cpc.2026.110221},
arxiv = {https://arxiv.org/abs/2509.18429},
}
```
The commands below illustrate how to perform a bifurcation analysis of buckling 3-D structure using `ff-bifbox`.

The implementation is based on a weak formulation developed from a total Lagrangian approach, which has been implemented in the equations file for this example: [eqns_douglas_jolivet_2026_structure.idp](./eqns_douglas_jolivet_2026_structure.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/douglas_jolivet_2026/structure
export nproc=4
```

# Example 1: 3-D buckling plate

1. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/douglas_jolivet_2026/eqns_douglas_jolivet_2026_structure.idp eqns.idp
ln -sf examples/douglas_jolivet_2026/settings_douglas_jolivet_2026_structure.idp settings.idp
```

2. Build meshes
```sh
FreeFem++-mpi -v 0 examples/douglas_jolivet_2026/mesh_structure.md -mo $workdir/structure
```

3. Compute $S$ branch at various $P$ values using predictor--corrector continuation method
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -mi structure_S.mesh -fo structure_S -param P -nu 0.3 -h0 0.1 -tgv -2 -paramtarget 1.5e-7 -maxcount -1
```

4. Compute fold points on the $S$ curve
```sh
cd "$workdir" && declare -a foldguesslist=(structure_S_*specialpt.base) && cd -
for guess in "${foldguesslist[@]}"; do
  out=$(echo "$guess" | awk -F'specialpt' '{print $1}')
  ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $guess -fo $out -param P -tgv -2
done
```

5. (OPTIONAL) Compute the relevant portion of the eigenvalue spectrum of the $S$ branch at each point for each possible symmetry (i.e. $S$, $H_x$, $H_y$, $R_z$).
```sh
cd "$workdir" && declare -a baselist=(structure_S_*[0-9].base) && cd -
for baseflow in "${baselist[@]}"; do
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_SS -eps_target 0.1+0i -eps_nev 3 -eps_gen_hermitian -sym 0,0
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_SHx -eps_target 0.1+0i -eps_nev 3 -eps_gen_hermitian -sym 0,1
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_SHy -eps_target 0.1+0i -eps_nev 3 -eps_gen_hermitian -sym 1,0
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_SRz -eps_target 0.1+0i -eps_nev 3 -eps_gen_hermitian -sym 1,1
done
```

6. Compute eigenmodes near exchanges of stability and iterate to find pitchfork bifurcations associated with symmetry breaking of the $S$ state
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_S_4.base -fo structure_SHx -eps_target 0.1+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 0,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_SHx.mode -fo structure_SHx -param P -tgv -2 -zero 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_S_30.base -fo structure_SHx2 -eps_target 0.1+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 0,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_SHx2.mode -fo structure_SHx2 -param P -tgv -2 -zero 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_S_9.base -fo structure_SHy -eps_target 0.1+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_SHy.mode -fo structure_SHy -param P -tgv -2 -zero 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_S_20.base -fo structure_SHy2 -eps_target 0.1+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_SHy2.mode -fo structure_SHy2 -param P -tgv -2 -zero 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_S_14.base -fo structure_SRz -eps_target 0.1+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_SRz.mode -fo structure_SRz -param P -tgv -2 -zero 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_S_22.base -fo structure_SRz2 -eps_target 0.1+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_SRz2.mode -fo structure_SRz2 -param P -tgv -2 -zero 1
```

7. Extend the computed bifurcations across the domain according to their symmetry and trace them along $P$
```sh
ff-mpirun -np 1 examples/douglas_jolivet_2026/extend_structure.md -dir $workdir -fi structure_SHx.hopf -v 0 -fo structure_SHx
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -mi structure_Hx.mesh -fi structure_SHxstart.base -fi2 structure_SHxbranch.base -fo structure_Hx -param P -h0 -0.0071825 -tgv -2 -snes_rtol 0 -mono 0 -maxcount 35
ff-mpirun -np 1 examples/douglas_jolivet_2026/extend_structure.md -dir $workdir -fi structure_SHy.hopf -v 0 -fo structure_SHy
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -mi structure_Hy.mesh -fi structure_SHystart.base -fi2 structure_SHybranch.base -fo structure_Hy -param P -h0 -0.03125 -tgv -2 -snes_rtol 0 -mono 0 -maxcount 40
ff-mpirun -np 1 examples/douglas_jolivet_2026/extend_structure.md -dir $workdir -fi structure_SRz.hopf -v 0 -fo structure_SRz
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -mi structure_Rz.mesh -fi structure_SRzstart.base -fi2 structure_SRzbranch.base -fo structure_Rz -param P -h0 0.015625 -tgv -2 -snes_rtol 0 -maxcount 41
``` 

8. Compute fold points on the bifurcation curves
```sh
cd "$workdir" && declare -a foldguesslist=(structure_[H,R][x,y,z]_*specialpt.base) && cd -
for guess in "${foldguesslist[@]}"; do
  out=$(echo "$guess" | awk -F'specialpt' '{print $1}')
  ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $guess -fo $out -param P -tgv -2
done
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi structure_Hy_34.base -fo structure_Hy_35 -param P -tgv -2
```

9. (OPTIONAL) Compute the relevant portion of the eigenvalue spectrum of the $R_z$ branch at each point for each broken symmetry.
```sh
cd "$workdir" && declare -a baselist=(structure_Hx_*[0-9].base) && cd -
for baseflow in "${baselist[@]}"; do
    ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_HxHx -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 0,1
    ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_HxA -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 1,1
done
cd "$workdir" && declare -a baselist=(structure_Hy_*[0-9].base) && cd -
for baseflow in "${baselist[@]}"; do
    ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_HyHy -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 1,0
    ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_HyA -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 1,1
done
cd "$workdir" && declare -a baselist=(structure_Rz_*[0-9].base) && cd -
for baseflow in "${baselist[@]}"; do
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_RzRz -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 0,0
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_RzA -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 1,1
done
```

10. Compute eigenmodes near exchanges of stability and iterate to find pitchfork bifurcations associated with symmetry breaking of the $H_x$, $H_y$, and $R_z$ states
- Bifurcations from $H_x$ to $A$ states:
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_Hx_14.base -fo structure_HxA -eps_target 0+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_HxA.mode -fo structure_HxA -param P -tgv -2 -zero 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_Hx_16.base -fo structure_HxA2 -eps_target 0+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_HxA2.mode -fo structure_HxA2 -param P -tgv -2 -zero 1
```
- Bifurcations from $H_y$ to $A$ states:
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_Hy_17.base -fo structure_HyA -eps_target 0+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_HyA.mode -fo structure_HyA -param P -tgv -2 -zero 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_Hy_26.base -fo structure_HyA2 -eps_target 0+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_HyA2.mode -fo structure_HyA2 -param P -tgv -2 -zero 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_Hy_34.base -fo structure_HyA3 -eps_target 0+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_HyA3.mode -fo structure_HyA3 -param P -tgv -2 -zero 1
```
- Bifurcations from $R_z$ to $A$ states:
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi structure_Rz_23.base -fo structure_RzA -eps_target 0+0i -eps_nev 1 -strict 1 -eps_gen_hermitian -sym 1,1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi structure_RzA.mode -fo structure_RzA -param P -tgv -2 -zero 1
```

11. Trace bifurcated solutions along $P$
- From $H_x$ to $A$
```sh
ff-mpirun -np 1 examples/douglas_jolivet_2026/extend_structure.md -dir $workdir -fi structure_HxA.hopf -v 0 -fo structure_HxA -pv 1
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -mi structure_A.mesh -fi structure_HxAstart.base -fi2 structure_HxAbranch.base -fo structure_A -param P -h0 -0.015625 -tgv -2 -snes_rtol 0 -mono 0 -maxcount 29
```
- From $H_y$ to $A$
```sh
ff-mpirun -np 1 examples/douglas_jolivet_2026/extend_structure.md -dir $workdir -fi structure_HyA.hopf -v 0 -fo structure_HyA
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -mi structure_A.mesh -fi structure_HyAstart.base -fi2 structure_HyAbranch.base -fo structure_A2 -param P -h0 0.015625 -tgv -2 -snes_rtol 0 -mono 0 -maxcount 36
ff-mpirun -np 1 examples/douglas_jolivet_2026/extend_structure.md -dir $workdir -fi structure_HyA3.hopf -v 0 -fo structure_HyA3 -pv 1
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -mi structure_A.mesh -fi structure_HyA3start.base -fi2 structure_HyA3branch.base -fo structure_A3 -param P -h0 0.0078125 -tgv -2 -snes_rtol 0 -mono 0 -maxcount 20
```
- From $R_z$ to $A$
```sh
ff-mpirun -np 1 examples/douglas_jolivet_2026/extend_structure.md -dir $workdir -fi structure_RzA.hopf -v 0 -fo structure_RzA
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -mi structure_A.mesh -fi structure_RzAstart.base -fi2 structure_RzAbranch.base -fo structure_A4 -param P -h0 -0.03125 -tgv -2 -snes_rtol 0 -maxcount 23
``` 

12. Compute fold points on the bifurcation curves
```sh
cd "$workdir" && declare -a foldguesslist=(structure_A*specialpt.base) && cd -
for guess in "${foldguesslist[@]}"; do
  out=$(echo "$guess" | awk -F'specialpt' '{print $1}')
  ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $guess -fo $out -param P -tgv -2
done
```

13. (OPTIONAL) Compute the relevant portion of the eigenvalue spectrum 
```sh
cd "$workdir" && declare -a baselist=(structure_A_*[0-9].base) && cd -
for baseflow in "${baselist[@]}"; do
    ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_A -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 1,1
done
cd "$workdir" && declare -a baselist=(structure_A2_*[0-9].base) && cd -
for baseflow in "${baselist[@]}"; do
    ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_A2 -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 1,1
done
cd "$workdir" && declare -a baselist=(structure_A3_*[0-9].base) && cd -
for baseflow in "${baselist[@]}"; do
  ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi $baseflow -so structure_A3 -eps_target 0.1+0i -eps_nev 5 -eps_gen_hermitian -sym 1,1
done
```
