# 2D Incompressible Flow Example: Sipp and Lebedev. JFM. (2007)
This file shows an example `ff-bifbox` workflow for reproducing the results of the paper:
```bibtex
@article{sipp_lebedev_2007,
  title={Global stability of base and mean flows: a general approach and its applications to cylinder and open cavity flows},
  volume={593},
  DOI={10.1017/S0022112007008907},
  journal={Journal of Fluid Mechanics},
  publisher={Cambridge University Press},
  author={Sipp, Denis and Lebedev, Anton},
  year={2007},
  pages={333–358}
}
```
The commands below illustrate how to perform a weakly nonlinear analysis of the 2D incompressible flow around a cylinder and an open cavity using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{1}{Re}\frac{\partial^2u_i}{\partial x_j^2} &= 0 \\
\frac{\partial u_i}{\partial x_i} &= 0
\end{align*}
$$

together with the boundary conditions:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $u_x=1$, $u_y=0$ |
| Wall, $\Gamma_w$ | $u_x=u_y=0$ |
| Slip, $\Gamma_s$ | $\frac{\partial u_x}{\partial y}=u_y=0$ |
| Axis, $\Gamma_a$| $\frac{\partial u_x}{\partial y}=u_y=0$, if symmetric |
| Axis, $\Gamma_a$| $u_x=\frac{\partial u_y}{\partial y}=0$, if asymmetric |
| Outlet, $\Gamma_o$ | $\frac{1}{Re}\frac{\partial u_i}{\partial x}-p\hat{e}_x= 0$ |

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the planar domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_a+\Gamma_s+\Gamma_o$. Solutions $\vec{q}=\left(u_i,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{p}\right)^T$,

$$
\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} = 0.

$$

This weak formulation has been implemented in the equations file for this example: [eqns_sipp_lebedev_2007.idp](./eqns_sipp_lebedev_2007.idp).

Note that, in this example of Sipp and Lebedev, viscosity is parameterized by $1/Re$ instead of $Re$ in order to make the equation system linear with respect to the control parameter. Though such scalings do improve the performance of predictor-corrector methods and weakly-nonlinear analysis, `ff-bifbox` does not require the system to be linear in the parameters.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/sipp_lebedev_2007/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/sipp_lebedev_2007/eqns_sipp_lebedev_2007.idp eqns.idp
ln -sf examples/sipp_lebedev_2007/settings_sipp_lebedev_2007.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from `.geo` files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/sipp_lebedev_2007 -dir $workdir -mi cylinder.geo
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/sipp_lebedev_2007 -dir $workdir -mi cavity.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/sipp_lebedev_2007/cylinder.md -mo $workdir/cylinder
FreeFem++-mpi -v 0 examples/sipp_lebedev_2007/cavity.md -mo $workdir/cavity
```

## Perform parallel computations using `ff-bifbox`
### Zeroth order
1. Compute base states on the created meshes at $Re=10$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi cylinder.msh -fo cylinder -1/Re 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi cavity.msh -fo cavity -1/Re 0.1
```

2. Continue base state along the parameter $1/Re$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi cylinder.base -fo cylinder -param 1/Re -h0 -1 -scount 2 -maxcount 8 -mo cylinder -thetamax 5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi cavity.base -fo cavity -param 1/Re -h0 -1 -scount 4 -maxcount 16 -mo cavity
```

3. Compute base states at $Re=50$ (cylinder) and $Re=4000$ (cavity) with guess from continuation
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_8.base -fo cylinder50 -1/Re 0.021
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cavity_16.base -fo cavity4000 -1/Re 0.00025
```

### First & second order
1. Compute leading direct eigenmode at $Re=50$ (cylinder) and $Re=4000$ (cavity)
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi cylinder50.base -fo cylinder50 -eps_target 0.1+0.8i -sym 1 -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi cavity4000.base -fo cavity4000 -eps_target 0.1+8.0i -sym 0 -eps_pos_gen_non_hermitian
```
NOTE: Here, the `-sym` argument specifies the asymmetric (1) or symmetric (0) reflective symmetry across the boundary `BCaxis`.

2. Compute the critical point and critical base/direct/adjoint solution
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder50.mode -fo cylinder -param 1/Re -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cavity4000.mode -fo cavity -param 1/Re -nf 0
```

3. Adapt the mesh to the critical solution, save `.vtu` files for Paraview
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder.hopf -fo cylinderadapt -mo cylinderhopf -adaptto bda -param 1/Re -thetamax 5 -pv 1 -wnl 1 
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cavity.hopf -fo cavityadapt -mo cavityhopf -adaptto bda -param 1/Re -pv 1 -wnl 1
```
NOTE: the normalizations for the direct and adjoint eigenmodes (and therefore also the weakly-nonlinear corrections) used by `ff-bifbox` are different than the normalizations used by Sipp and Lebedev. This causes the results to differ by a complex scaling factor.


### Harmonic Balance
1. Continue periodic orbit from initial Hopf bifurcations using 2nd-order Harmonic Balance (Caution: memory intensive!)
```sh
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi cylinder.hopf -fo cylinderNh2 -Nh 2 -mo cylinderporb -param 1/Re -thetamax 5 -h0 -1 -scount 5 -maxcount 10
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi cavity.hopf -fo cavityNh2 -Nh 2 -mo cavityporb -param 1/Re -h0 -1 -scount 4 -maxcount 8
```

2. Compute periodic orbits at $Re=50$ (cylinder) and $Re=5000$ (cavity) using 3rd-order Harmonic Balance with block Jacobi solver (Caution: memory intensive!)
```sh
ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi cylinderNh2_10.porb -fo cylinder50Nh3 -Nh 3 -1/Re 0.02 -blocks 3
ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi cavityNh2_8.porb -fo cavity5000Nh3 -Nh 3 -1/Re 0.0002 -blocks 3
```
NOTE: Sipp & Lebedev do not perform harmonic balance analysis. See Fabre et al., Appl. Mech. Rev. 2018 and Meliga, JFM, 2017 for reference results from the cylinder and cavity geometries, respectively.


### Time-domain nonlinear simulation
1. Compute time series over first 10 time units for cavity case after step increase in Reynolds number, adapt mesh to solution, export files to Paraview. 
```sh
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir -fi cavity4000.base -fo cavity -1/Re 0.0002 -tsdt 0.01 -mo cavitytimeseries -scount 5 -maxcount 1000 -pv 1
```