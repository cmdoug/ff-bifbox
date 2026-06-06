# 2D Incompressible Flow Example: Garnaud, PhD Thesis, (2012)
This file shows an example `ff-bifbox` workflow for reproducing the results in Chapter 6.4 of the thesis:
```bibtex
@phdthesis{garnaud_2012,
  TITLE = {{Modes, transient dynamics and forced response of circular jets}},
  AUTHOR = {Garnaud, Xavier},
  URL = {https://theses.hal.science/tel-00740133},
  SCHOOL = {{Ecole Polytechnique}},
  YEAR = {2012},
}
```
The commands below illustrate how to perform a resolvent analysis of an incompressible laminar axisymmetric jet using `ff-bifbox`.

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
| Inlet, $\Gamma_i$ | $u_x=\tanh\left(5-5r\right)$, $u_r=0$ |
| Wall, $\Gamma_w$ | $u_x=u_r=0$ |
| Axis, $\Gamma_a$| $\frac{\partial u_x}{\partial r}=u_r=0$ |
| Open, $\Gamma_o$ | $\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\hat{n}_j-p\hat{n}_i = 0$ |

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_a+\Gamma_o$. Solutions $\vec{q}=\left(u_i,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{p}\right)^T$,

$$
\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} = 0.
$$

This weak formulation has been implemented in the equations file for this example: [eqns_garnaud_2012.idp](./eqns_garnaud_2012.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/garnaud_2012/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/garnaud_2012/eqns_garnaud_2012.idp eqns.idp
ln -sf examples/garnaud_2012/settings_garnaud_2012.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from `.geo` files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/garnaud_2012 -dir $workdir -mi jet.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/garnaud_2012/jet.md -mo $workdir/jet
```

## Perform parallel computations using `ff-bifbox`
### Laminar base flow
1. Compute base states on the created mesh at $Re=10$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi jet.msh -fo jet -1/Re 0.1
```

2. Continue base state along the parameter $1/Re$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi jet.base -fo jet -param 1/Re -h0 -1 -scount 3 -maxcount 15 -mo jet
```

3. Compute base state at $Re=1000$ with guess from continuation
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi jet_15.base -fo jet1000 -1/Re 0.001
```

4. Adapt mesh to $Re=1000$ solution and recompute base state, save `.vtu` file for Paraview
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi jet1000.base -fo jet1000adapt -mo jet1000adapt -pv 1
```

### Resolvent analysis
1. Compute gains across $0.1\leq{}St\leq{}1.0$ frequency range $(\pi/10<\omega<\pi)$ in 10 increments
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi jet1000adapt.base -so jet1000adapt -omega 0.314159 -omegaf 3.14159 -nomega 10
```
2. Compute forcing/response modes at $St=[0.1, 0.48, 0.64, 0.95]$, save `.vtu` files for Paraview
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi jet1000adapt.base -fo jet1000_St0p1 -omega 0.314159 -pv 1
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi jet1000adapt.base -fo jet1000_St0p48 -omega 1.50796 -pv 1
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi jet1000adapt.base -fo jet1000_St0p64 -omega 2.01062 -pv 1
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi jet1000adapt.base -fo jet1000_St0p95 -omega 2.98451 -pv 1
```
