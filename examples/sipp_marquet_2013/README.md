# 2D Incompressible Flow Example: Sipp and Marquet. TCFD. (2013)
This file shows an example `ff-bifbox` workflow for reproducing the results of the paper:
```tex
@article{sipp_marquet_2013,
  title={Characterization of noise amplifiers with global singular modes: the case of the leading-edge flat-plate boundary layer},
  volume={27},
  DOI={10.1007/s00162-012-0265-y},
  journal={Theoretical and Computational Fluid Dynamics},
  publisher={Springer},
  author={Sipp, Denis and Marquet, Olivier},
  year={2013},
  pages={617–635}
}
```
The commands below illustrate how to perform a resolvent analysis of a flat-plate boundary layer using `ff-bifbox`.

Note that viscosity is parameterized by $1/Re$ instead of $Re$ in this example in order to make the equation system linear with respect to the control parameter. Though such scalings do improve the performance of predictor-corrector methods and weakly-nonlinear analysis, `ff-bifbox` does not require the system to be linear in the parameters.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/sipp_marquet_2013/data
export nproc=10
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/sipp_marquet_2013/eqns_sipp_marquet_2013.idp eqns.idp
ln -sf examples/sipp_marquet_2013/settings_sipp_marquet_2013.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
```sh
FreeFem++-mpi -v 0 examples/sipp_marquet_2013/flatplate.md -mo $workdir/flatplate_0
```

## Perform parallel computations using `ff-bifbox`
### Base flow
1. Compute base state on the created mesh at $Re=100$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi flatplate_0.msh -fo bl_0 -1/Re 0.01
```

2. Continue base state along the parameter $1/Re$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi bl_0.base -fo bl -param 1/Re -h0 -7 -scount 3 -maxcount -1 -paramtarget 0.000002 -mo flatplate -pv 1 -anisomax 2
```

3. Compute base state at $Re=6\times10^5$ with final guess from continuation and readapt the mesh to this solution.
```sh
cd $workdir && export lastfile=$(printf '%s\n' bl_*.base | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo bl_target -1/Re 0.000001666666666666666667 -snes_linesearch_type l2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi bl_target.base -fo bl_target_adapt -mo flatplate_target_adapt -pv 1 -anisomax 2 -hmax 0.1  -snes_linesearch_type l2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi bl_target_adapt.base -fo bl_Re6e5 -mo flatplate_Re6e5 -pv 1 -hmax 0.005 -snes_linesearch_type l2
```

### Resolvent analysis
4. Run resolvent analysis frequency sweep for optimal forcing/response at $Re=6\times10^5$. Choose a range of frequencies from $12 < \omega < 120$ (corresponding to $20 < F < 200$ as defined in the paper).
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi bl_Re6e5.base -so bl_Re6e5 -omega 12 -omegaf 120 -nomega 13 -eps_nev 1
```

5. Compute optimal and first suboptimal forcing and response modes at $Re=6\times10^5$ for $\omega = 60$ (i.e. $F = 100$)
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi bl_Re6e5.base -so bl_Re6e5 -fo bl_Re6e5 -omega 60.0 -eps_nev 2 -pv 1
```
