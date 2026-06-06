# 3D Incompressible Wake Flow Example: Moulin et al., (2019)
This file shows an example `ff-bifbox` workflow for reproducing the results of the study:
```bibtex
@article{moulin_etal_2019,
    Author = {Moulin, Johann and Jolivet, Pierre and Marquet, Olivier},
    Title = {{Augmented Lagrangian Preconditioner for Large-Scale Hydrodynamic Stability Analysis}},
    Year = {2019},
    Volume = {351},
    Pages = {718--743},
    DOI = {10.1016/j.cma.2019.03.052},
    Journal = {Computer Methods in Applied Mechanics and Engineering},
    Publisher = {Elsevier},
    Url = {https://github.com/prj-/moulin2019al}
}
```
The commands below illustrate how to run the perform a stability analysis of 3D wake behind a rectangular plate using `ff-bifbox` with the mAL preconditioner implemented in this study. Note that since this study leverages iterative methods, attempts to solve for bifurcation points (fold, hopf, etc.) will not work properly since the matrix is horribly ill-conditioned near such singularities. This method should only be used for base flow calculations, stability analysis, resolvent analysis, and/or time-domain simulations. (The latter of which could be done faster with other preconditioning strategies.)

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
| Inlet, $\Gamma_i$ | $`u_x=1$, $u_y=u_z=0`$ |
| Wall, $\Gamma_w$ | $`u_x=u_y=u_z=0`$ |
| Slip, $\Gamma_s$ | $`u_i\hat{n}_i=\frac{\partial u_i}{\partial x_l}\epsilon_{ijk}\hat{n}_j\hat{n}_l=0`$ |
| Outlet, $\Gamma_o$ | $`\frac{1}{Re}\frac{\partial u_i}{\partial x}-p\hat{e}_x= 0`$ |

The present implementation is based on a grad–div stabilized weak formulation of these equations. Test functions are introduced, and the equations are integrated over the Cartesian domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_s+\Gamma_o$. Solutions $\vec{q}=\left(u_i,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{p}\right)^T$,

$$
\begin{align*} 
&\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\right)_{\Omega} \\
&+ \left(\gamma\frac{\partial \check{u}_i}{\partial x_i}-\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} = 0.
\end{align*}
$$

This weak formulation, together with the mAL preconditioner approach, has been implemented in the equations file for this example: [eqns_moulin_etal_2019.idp](./eqns_moulin_etal_2019.idp).

IMPORTANT NOTE: The ability to solve 3 dimensional problems in `ff-bifbox` is still under development! In particular, 3D mesh adaptation with `mmg3d` may contain bugs. 

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```

2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/moulin_etal_2019/data
export nproc=4
```

3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/moulin_etal_2019/eqns_moulin_etal_2019.idp eqns.idp
ln -sf examples/moulin_etal_2019/settings_moulin_etal_2019.idp settings.idp
```

## Build initial meshes
In 3D, `ff-bifbox` uses `mshmet`+`mmg` for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine. In this example, we use the mesh `FlatPlate3D.mesh` provided by the authors, which has an aspect ratio of L = 2.5. The original mesh file can be accessed via GitHub at `https://github.com/prj-/moulin2019al`.

## Perform parallel computations using `ff-bifbox`

### Steady dynamics
1. Compute a base state on the mesh at $Re=50$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi FlatPlate3D.mesh -fo 3Dwake -1/Re 0.02 -gamma 0.6
```

2. Continue base state along $1/Re$ from $Re=50$ solution.
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi 3Dwake.base -param 1/Re -h0 -5 -kmax 4 -snes_max_it 20 -scount 2 -maxcount 4
```

3. Compute a base state on the mesh at $Re=100$ with guess from continuation
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi 3Dwake_4.base -fo 3Dwake100 -1/Re 0.01
```

### Unsteady dynamics
4. Compute leading eigenvalue at $Re=100$. This is very slow unless massively parallelized.
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi 3Dwake100.base -fo 3Dwake -eps_target 0.1+0.6i -eps_nev 5 -eps_ncv 15 -eps_tol 1e-6 -recycle 5 -shiftPrecon 1 -st_ksp_converged_reason -eps_pos_gen_non_hermitian
```

5. Compute optimal resolvent gain at $Re=50$. This is very slow unless massively parallelized.
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi 3Dwake.base -fo 3Dwake -omega 1 -recycle 5 -shiftPrecon 1 -eps_tol 1e-6
```
