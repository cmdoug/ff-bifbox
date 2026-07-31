# Low-Mach V-flame Example: Douglas & Wang, (2026)
This file shows an example `ff-bifbox` workflow for reproducing the results in the study:
```tex
@article{schulke_etal_2026,

}
```
The commands below illustrate how to perform a bifurcation analysis of a lean premixed V-flame in an axisymmetric annular jet using `ff-bifbox`.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/schulke_etal_2026/3D
export nproc=8
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/schulke_etal_2026/eqns_schulke_etal_2026_3D.idp eqns.idp
ln -sf examples/schulke_etal_2026/settings_schulke_etal_2026_3D.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from .geo files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/schulke_etal_2026 -dir $workdir -mi Vflame.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/schulke_etal_2026/Vflame.md -mo $workdir/Vflame
```

## Perform parallel computations using `ff-bifbox`
### Laminar base flow
1. Compute a non-reacting base state with reference parameters on the initial mesh.
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi Vflame.msh -fo nonreacting_0 -Re 70 -Tr 2.3333333333333333 -Da 0 -phi 0.8 -Dhc 128.23469709865222 -snes_linesearch_type l2 -snes_rtol 0 -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting_0.base -fo nonreacting_1 -Re 350 -mo nonreacting_1 -snes_linesearch_type l2 -snes_rtol 0 -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting_1.base -fo nonreacting_2 -Re 1500 -mo nonreacting_2 -snes_linesearch_type l2 -snes_rtol 0 -pv 1
```

2. Turn on chemistry and ignite the U0 = 2.2 m/s flow at an elevated centrebody temperature and lower combustion enthalpy. Then perform continuation back to reference parameters. Coarse meshes are used for computational efficiency and stabilizing artificial dissipation. 
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting_2.base -fo ignite_0 -Tr 3.3333333333333333 -Da 1090638.826627295 -Dhc 12.823469709865222 -mo ignite_0 -snes_rtol 0 -err 0.05 -pv 1
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi ignite_0.base -fo ignite -param Dhc -mo ignite -h0 50 -err 0.05 -scount 5 -paramtarget 128.23469709865222 -maxcount -1 -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_360.base -fo ignited -Tr 2.3333333333333333 -Dhc 128.23469709865222 -mo ignited -snes_rtol 0 -err 0.05 -snes_linesearch_type l2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignited.base -fo phi_0p8_Re_1500 -pv 1 -snes_rtol 0 -snes_linesearch_type l2 -mo phi_0p8_Re_1500 -err 0.01 -hmin 1e-5 -hmax 0.2
```

3. Save base flows over a range of $Re$. Dissipation from mesh coarsening is used to aid convergence at each step before refining the coarse solutions on the reference mesh.
```sh
export phi=8
for Re in {1600..3500..100}; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi phi_0p"$phi"_Re_"$(($Re-100))".base -fo U0inc -Re "$Re" -snes_linesearch_type l2 -mo U0inc -err 0.05
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U0inc.base -fo U0inc -Re "$Re" -snes_linesearch_type l2 -mo U0inc -err 0.05
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U0inc.base -fo phi_0p"$phi"_Re_"$Re".base -Re "$Re" -pv 1 -snes_rtol 0 -snes_linesearch_type l2 -mo phi_0p"$phi"_Re_"$Re" -err 0.01 -hmin 1e-5 -hmax 0.2
done
for Re in {1500..3500..100}; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi phi_0p"$phi"_Re_"$Re".base -fo phi_0p"$phi"_Re_"$Re".base -Re "$Re" -pv 1 -snes_rtol 0
done
```

### Global linear analysis
4. Compute global eigenspectra
```sh
export phi=8
for m in {0..2}; do
for Re in {1600..3200..100}; do
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi phi_0p"$phi"_Re_"$Re".base -so phi_0p"$phi"_Reswp_m_"$m" -eps_nev 25 -eps_target 0.1+1i -ntarget 5 -targetf 0.1+9i  -sym $m -eps_tol 2.22e-14
done
done
```