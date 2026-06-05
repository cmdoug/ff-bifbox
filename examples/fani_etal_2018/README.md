# 2D Compressible Flow Example: Fani et al. PoF. (2018)
This file shows an example `ff-bifbox` workflow for reproducing the results of the paper:
```bibtex
@article{fani_etal_2018,
    author = {Fani, A. and Citro, V. and Giannetti, F. and Auteri, F.},
    title = "{Computation of the bluff-body sound generation by a self-consistent mean flow formulation}",
    journal = {Physics of Fluids},
    volume = {30},
    number = {3},
    year = {2018},
    doi = {10.1063/1.4997536},
}
```
The commands below illustrate how to analyze a 2D compressible flow around a cylinder using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial\rho}{\partial t} + \rho \frac{\partial u_i}{\partial x_i} + u_i\frac{\partial\rho}{\partial x_i} - \tilde{\beta}\left(\rho-1\right)&= 0\\
\rho\frac{\partial u_i}{\partial t} + \rho u_j\frac{\partial u_i}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{1}{Re}\frac{\partial\tau_{ij}}{\partial x_j} - \tilde{\beta}\left(u_i-\hat{e}_x\right) &= 0 \\
\rho\frac{\partial T}{\partial t} + \rho u_i\frac{\partial T}{\partial x_i} +\left(\gamma-1\right)\rho T\frac{\partial u_i}{\partial x_i} - \frac{\gamma\left(\gamma-1\right)M^2}{Re}\tau_{ij}\frac{\partial u_i}{\partial x_j} \\- \frac{\gamma}{Pr Re}\frac{\partial^2T}{\partial x_i^2} - \tilde{\beta}\left(T-1\right) &= 0
\end{align*}
$$

where:
- $`\tau_{ij}=\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}-\frac{2}{3}\delta_{ij}\frac{\partial u_k}{\partial x_k}`$
- $`\rho = \frac{1 + \gamma M^2p}{T}`$
- $`\tilde{\beta}\left(x,y\right)=\begin{cases}
0,& \text{if } x_3 \leq x \leq x_4 \text{ and } \|y\| \leq y_2 \\\\
\left|1-\frac{1}{M}\right|f\left(x_3,x\right),& \text{if } x < x_3 \text{ and } \|y\| \leq y_2 \\\\
\left|1+\frac{1}{M}\right|f\left(x,x_4\right), & \text{if } x > x_4 \text{ and } \|y\| \leq y_2 \\\\
\tilde{\beta}\left(x,y_2\right)+\left|\frac{1}{M}\right|f\left(y,y_4\right),& \text{if } \|y\| > y_2
\end{cases}`$
- $`f\left(a,b\right)=\frac{2\alpha}{l_s^2}\left(a-b\right)`$.

The boundary conditions are:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $u_x=T=1$, $u_y=0$ |
| Wall, $\Gamma_w$| $u_x=u_y=\frac{\partial T}{\partial x_i}\hat{n}_i=0$ |
| Axis, $\Gamma_a$| $\frac{\partial u_x}{\partial y}=u_y=\frac{\partial T}{\partial y}=0$|
| Lateral, $\Gamma_l$ | $\frac{\partial u_x}{\partial y}=u_y=\frac{\partial T}{\partial y}=0$ |
| Outlet, $\Gamma_o$ | $\frac{1}{Re}\tau_{ix}-p\hat{e}_x = \frac{\partial T}{\partial x}=0$ |

NOTE: the original paper also imposes a $\rho=1$ condition along $\Gamma_i$. However, this is not needed, and, in fact, it overconstrains the system due to the algebraic ideal gas equation of state. (There are really only two dynamically independent thermodynamic variables.) The approach taken here corrects this.

The present implementation is based on a weak formulation of these equations. Density is eliminated using the equation of state, test functions are introduced, and the equations are integrated over the planar domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_a+\Gamma_l+\Gamma_o$. Solutions $\vec{q}=\left(u_i,T,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{T},\check{p}\right)^T$,

$$
\begin{align*} 
&\left(\check{u}_i,\frac{1 + \gamma M^2p}{T}\left[\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right] + \tilde{\beta}\left(u_i-\hat{e}_x\right)\right)_{\Omega} - \left(\frac{\partial \check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\tau_{ij}\right)_{\Omega} \\
&+ \left(\check{p},\gamma M^2\left[\frac{\partial p}{\partial t} + u_i\frac{\partial p}{\partial x_i}\right] - \frac{1 + \gamma M^2p}{T}\left[\frac{\partial T}{\partial t} + u_i\frac{\partial T}{\partial x_i} - T\frac{\partial u_i}{\partial x_i}\right] + \tilde{\beta}\left(1 + \gamma M^2p-T\right)\right)_{\Omega}\\
&+\left(\check{T},\frac{1 + \gamma M^2p}{T}\left[\frac{\partial T}{\partial t} + u_i\frac{\partial T}{\partial x_i} + \left(\gamma-1\right)T\frac{\partial u_i}{\partial x_i}\right] - \frac{\gamma\left(\gamma-1\right)M^2}{Re}\tau_{ij}\frac{\partial u_i}{\partial x_j}+\tilde{\beta}\left(T-1\right)\right)_{\Omega} \\
&+\left(\frac{\partial\check{T}}{\partial x_i},\frac{\gamma}{Pr Re}\frac{\partial T}{\partial x_i}\right)_{\Omega} = 0.
\end{align*}
$$

This weak formulation has been implemented in the equations file for this example: [eqns_fani_etal_2018.idp](./eqns_fani_etal_2018.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/fani_etal_2018/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/fani_etal_2018/eqns_fani_etal_2018.idp eqns.idp
ln -sf examples/fani_etal_2018/settings_fani_etal_2018.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from `.geo` files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/fani_etal_2018 -dir $workdir -mi cylinder.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/fani_etal_2018/cylinder.md -mo $workdir/cylinder
```

## Perform parallel computations using `ff-bifbox`
### Zeroth order
1. Compute base state on the created mesh at $Re=10$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi cylinder.msh -fo cylinder -1/Re 0.1 -1/Pr 1.38888888889 -Ma^2 0.04 -gamma 1.4
```

2. Continue base state along the parameter $1/Re$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi cylinder.base -fo cylinder -param 1/Re -h0 -1 -scount 2 -maxcount 14 -mo cylinder -thetamax 5
```

3. Compute base states at $Re\sim50$ and $Re=150$ with guesses from continuation
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_8.base -fo cylinder50 -1/Re 0.021
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_14.base -fo cylinder150 -1/Re 0.0066666666667
```

4. Adapt mesh to the $Re=150$ solution with a maximum triangle size restriction
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder150.base -fo cylinder150 -mo cylinder150 -thetamax 5 -hmax 5 -pv 1
```

### First order
1. Compute leading direct eigenmode at $Re\sim50$ and $Re=150$
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi cylinder50.base -fo cylinder50 -eps_target 0.1+0.7i -sym 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi cylinder150.base -fo cylinder150 -eps_target 0.2+0.8i -sym 1 -pv 1
```
NOTE: Here, the `-sym` argument specifies the asymmetric (1) or symmetric (0) reflective symmetry across the boundary `BCaxis`.

2. Compute the critical point and critical base/direct/adjoint solution
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder50.mode -fo cylinder -param 1/Re -nf 0
```

3. Adapt the mesh to the critical solution, save `.vtu` files for Paraview
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder.hopf -fo cylinder -mo cylinderhopf -adaptto bda -param 1/Re -thetamax 5 -pv 1
```

4. Continue the neutral Hopf curve in the $(1/Re,Ma^2)$-plane with adaptive remeshing
```sh
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi cylinder.hopf -fo cylinder -mo cylinderhopf -adaptto bda -thetamax 5 -param Ma^2 -param2 1/Re -h0 -1 -scount 3 -maxcount 12
```
NOTE: the signs and normalizations of the normal form coefficients used in `hopfcompute.md` are different than those of the Stuart-Landau coefficients in [Sipp and Lebedev JFM (2007)](../sipp_lebedev_2007/).

5. Continue the branch of periodic solutions emanating from the Hopf point along $1/Re$ using harmonic balance.
```sh
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi cylinder.hopf -fo cylinder -mo cylinderporb -adaptto 01 -thetamax 5 -param 1/Re -h0 -1 -scount 5 -maxcount -1 -paramtarget 0.00666667
```
NOTE: the formulation in `ff-bifbox` is more fully self-consistent, and does not neglect the unsteady nonlinear interactions as in the original paper. 