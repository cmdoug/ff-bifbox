# 3D Incompressible Wake Flow Example: Marquet and Larsson, (2015)
This file shows an example `ff-bifbox` workflow for reproducing the results of the study:
```bibtex
@article{marquet_larsson_2015,
title = {Global wake instabilities of low aspect-ratio flat-plates},
journal = {European Journal of Mechanics - B/Fluids},
volume = {49},
pages = {400-412},
year = {2015},
note = {Trends in Hydrodynamic Instability in honour of Patrick Huerre's 65th birthday},
issn = {0997-7546},
doi = {10.1016/j.euromechflu.2014.05.005},
author = {O. Marquet and M. Larsson},
}
```
The commands below illustrate how to run the perform a stability analysis of 3D wake behind a rectangular plate using `ff-bifbox`. Note that unlike the example in `examples/moulin_etal_2019`, this study leverages direct methods, so bifurcation points (fold, hopf, etc.) can be located.

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
| Inlet, $\Gamma_i$ | $u_x=1$, $u_y=u_z=0$ |
| Wall, $\Gamma_w$ | $u_x=u_y=u_z=0$ |
| xz-plane, $\Gamma_y$ | $`\frac{\partial u_x}{\partial y}=u_y=\frac{\partial u_z}{\partial y}=0`$ if symmetric |
| xz-plane, $\Gamma_y$ | $`u_x=\frac{\partial u_y}{\partial y}=u_z=0`$ if asymmetric |
| xy-plane, $\Gamma_z$ | $`\frac{\partial u_x}{\partial z}=\frac{\partial u_y}{\partial z}=u_z=0`$ if symmetric |
| xy-plane, $\Gamma_z$ | $`u_x=u_y=\frac{\partial u_z}{\partial z}=0`$ if asymmetric |
| Slip, $\Gamma_s$ | $`u_i\hat{n}_i=\frac{\partial u_i}{\partial x_l}\epsilon_{ijk}\hat{n}_j\hat{n}_l=0`$ |
| Outlet, $\Gamma_o$ | $`\frac{1}{Re}\frac{\partial u_i}{\partial x}-p\hat{e}_x= 0`$ |

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the Cartesian domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_y+\Gamma_z+\Gamma_s+\Gamma_o$. Solutions $\vec{q}=\left(u_i,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{p}\right)^T$,

$$
\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} = 0.
$$

This weak formulation has been implemented in the equations file for this example: [eqns_marquet_larsson_2015.idp](./eqns_marquet_larsson_2015.idp).

NOTES: 
- This code uses computational coordinates that differ from the physical coordinates by a scaling factor related to the parameter $L$ (see the `Z()` macro in [settings_marquet_larsson_2015.idp](./settings_marquet_larsson_2015.idp)). Note that ParaView files are exported using the physical coordinates.

- The ability to solve 3 dimensional problems in `ff-bifbox` is still under development. In particular, 3D mesh adaptation with `mmg3d` may contain bugs. 

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```

2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/marquet_larsson_2015/data
export nproc=4
```

3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/marquet_larsson_2015/eqns_marquet_larsson_2015.idp eqns.idp
ln -sf examples/marquet_larsson_2015/settings_marquet_larsson_2015.idp settings.idp
```

## Build initial meshes
In 3D, `ff-bifbox` uses `mshmet`+`mmg` for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine. The example code here does not include arguments for mesh adaptation, as mesh adaptation in 3D using `mshmet`+`mmg` is not as robust as in 2D with `adaptmesh`.
#### Build initial mesh directly from `.geo` files using Gmsh
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/marquet_larsson_2015 -dir $workdir -mi plate.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).

## Perform parallel computations using `ff-bifbox`

### $y$-antisymmetric, $z$-symmetric mode
1. Compute a base state on the mesh at $Re=60$, $L=6$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi plate.mesh -fo wakeL6Re20 -1/Re 0.05 -L 6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi wakeL6Re20.base -fo wakeL6Re60 -1/Re 0.01666666666666666666 -L 6
```

2. Compute the leading eigenmode at $Re=60$, $L=6$ that is anti-symmetric along $y$ and symmetric along $z$.
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi wakeL6Re60.base -fo wakeL6Re60yAzS -eps_target 0.1+0.6i -sym 1,0 -eps_pos_gen_non_hermitian
```

3. Compute the critical point and critical base/direct/adjoint solution
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi wakeL6Re60yAzS.mode -fo wakeL6Re60yAzS -param 1/Re -nf 0
```

4. Continue the neutral Hopf curve in the $(1/Re, L)$-plane
```sh
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi wakeL6Re60yAzS.hopf -fo wakeL6Re60yAzS -param L -param2 1/Re -h0 -4 -scount 4 -maxcount 12
```

### $y$-symmetric, $z$-antisymmetric mode
5. Compute a base state on the mesh at $Re=100$, $L=3$ from guess at $Re=60$, $L=6$
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi wakeL6Re60.base -fo wakeL3Re100 -1/Re 0.01 -L 3
```

6. Compute the leading eigenmode at $Re=100$, $L=3$ that is symmetric along $y$ and anti-symmetric along $z$.
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi wakeL3Re100.base -fo wakeL3Re100ySzA -eps_target 0.1+0.3i -sym 0,1 -eps_pos_gen_non_hermitian
```

7. Compute the critical point and critical base/direct/adjoint solution
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi wakeL3Re100ySzA.mode -fo wakeL3Re100ySzA -param 1/Re -nf 0
```

8. Continue the neutral Hopf curve in the $(1/Re, L)$-plane
```sh
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi wakeL3Re100ySzA.hopf -fo wakeL3Re100ySzA -param L -param2 1/Re -h0 -4 -scount 4 -maxcount 12
```

### $y$-symmetric, $z$-antisymmetric mode
9. Compute a base state on the mesh at $Re\sim105$, $L=1.5$ from guess at $Re=100$, $L=3$
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi wakeL3Re100.base -fo wakeL1p5Re105 -1/Re 0.0095 -L 1.5
```

10. Compute the leading stationary eigenmode at $Re\sim105$, $L=1.5$ that is antisymmetric along $y$ and symmetric along $z$.
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi wakeL1p5Re105.base -fo wakeL1p5Re105yAzS -eps_target 0.1+0.0i -sym 1,0 -eps_pos_gen_non_hermitian
```

11. Compute the critical point and critical base/direct/adjoint solution
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi wakeL1p5Re105yAzS.mode -fo wakeL1p5Re105yAzS -param 1/Re -nf 0 -zero 1
```

12. Continue the neutral Hopf curve in the $(1/Re, L)$-plane
```sh
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi wakeL1p5Re105yAzS.hopf -fo wakeL1p5Re105yAzS -param L -param2 1/Re -h0 4 -scount 4 -maxcount 12 -zero 1
```

### double-Hopf point for $y$/$z$ symmetry/antisymmetry switch
13. Compute a base state on the mesh at $Re=100$, $L=2.5$ from guess at $Re=100$, $L=3$
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi wakeL3Re100.base -fo wakeL2p5Re100 -1/Re 0.01 -L 2.5
```

14. Compute the leading stationary eigenmodes.
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi wakeL2p5Re100.base -fo wakeL2p5Re100yAzS -eps_target 0.1+0.5i -sym 1,0 -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi wakeL2p5Re100.base -fo wakeL2p5Re100ySzA -eps_target 0.1+0.25i -sym 0,1 -eps_pos_gen_non_hermitian
```

15. Compute the critical point and critical base/direct/adjoint solution
```sh
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi wakeL2p5Re100yAzS.mode -fi2 wakeL2p5Re100ySzA.mode -fo wakeL2p5Re100 -param 1/Re -param2 L
```

### Hopf-pitchfork point for $y$/$z$ symmetry/antisymmetry switch
16. Compute a base state on the mesh at $Re=105$, $L=2$ from guess at $Re=100$, $L=3$
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi wakeL2pRe100.base -fo wakeL3Re95 -1/Re 0.0095 -L 2
```

17. Compute the leading stationary eigenmodes.
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi wakeL2p5Re100.base -fo wakeL2p5Re100yAzS -eps_target 0.1+0.5i -sym 1,0 -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi wakeL2p5Re100.base -fo wakeL2p5Re100ySzA -eps_target 0.1+0.25i -sym 0,1 -eps_pos_gen_non_hermitian
```

18. Compute the critical point and critical base/direct/adjoint solution
```sh
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi wakeL2p5Re100yAzS.mode -fi2 wakeL2p5Re100ySzA.mode -fo wakeL2p5Re100 -param 1/Re -param2 L
```
