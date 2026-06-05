# 2D Incompressible Swirling Flow Example: Meliga et al, JFM, (2012)
This file shows an example `ff-bifbox` workflow for reproducing the results in the study:
```bibtex
@article{meliga_etal_2012,
  title={A weakly nonlinear mechanism for mode selection in swirling jets},
  volume={699},
  DOI={10.1017/jfm.2012.93},
  journal={Journal of Fluid Mechanics},
  author={Meliga, Philippe and Gallaire, François and Chomaz, Jean-Marc},
  year={2012},
  pages={216–262}}
```
The commands below illustrate how to perform a bifurcation analysis of an incompressible swirling flow using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{\partial}{\partial x_j}\left(\frac{1}{\tilde{Re}}\frac{\partial u_i}{\partial x_j}\right) &= 0 \\
\frac{\partial u_i}{\partial x_i} &= 0
\end{align*}
$$

where:
- $\tilde{Re}(r,z)=\begin{cases}
Re & \text{if } r \leq r_{\max} \text{ and } x \leq x_{\max}\\
Re+\left(Re_s-Re\right)\zeta\left(z, z_{\max}\right) & \text{if } r \leq r_{\max} \text{ and } x > x_{\max}\\
\tilde{Re}\left(r_{\max}, x\right)+\left(Re_s-\tilde{Re}\left(r_{\max},x\right)\right)\zeta\left(r, r_{\max}\right) & \text{if } r > r_{\max}
\end{cases}$
- $\zeta\left(a,b\right)=\frac{1}{2}=\frac{1}{2}\tanh\left\lbrace\tau\tan\left(-\frac{\pi}{2}+\pi\frac{\|a-b\|}{l}\right)\right\rbrace$.

The boundary conditions are:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $u_x=1$, $u_r=0$, $u_{\theta}=S\Psi(r)$ |
| Axis, $\Gamma_a$| $\frac{\partial u_x}{\partial r}=u_r=u_{\theta}=0$, if $m=0$ |
| Axis, $\Gamma_a$| $u_x=\frac{\partial u_r}{\partial r}=\frac{\partial u_{\theta}}{\partial r}=0$, if $\|m\|=1$ |
| Axis, $\Gamma_a$| $u_x=u_r=u_{\theta}=0$, if $\|m\|>1$ |
| Lateral, $\Gamma_l$ | $\frac{\partial u_x}{\partial r}=u_r=u_{\theta}=0$ |
| Open, $\Gamma_o$ | $\frac{1}{Re_s}\frac{\partial u_i}{\partial x}-p\hat{e}_x = 0$ |

where $\Psi(r)=\begin{cases}r(2-r^2) & \text{if }r \leq 1\\\frac{1}{r} & \text{if } r> 1\end{cases}$.

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_a+\Gamma_l+\Gamma_o$. Solutions $\vec{q}=\left(u_i,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{p}\right)^T$,

$$
\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\frac{1}{\tilde{Re}}\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} = 0.
$$

This weak formulation has been implemented in the equations file for this example: [eqns_meliga_etal_2012.idp](./eqns_meliga_etal_2012.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/meliga_etal_2012/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/meliga_etal_2012/eqns_meliga_etal_2012.idp eqns.idp
ln -sf examples/meliga_etal_2012/settings_meliga_etal_2012.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from `.geo` files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/meliga_etal_2012 -dir $workdir -mi vortex.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/meliga_etal_2012/vortex.md -mo $workdir/vortex
```

## Perform parallel computations using `ff-bifbox`

### Steady axisymmetric dynamics
1. Compute base states on the created mesh at $Re=200$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi vortex.msh -fo vortex -1/Re 0.005 -S 0
```

2. Continue base state along the parameter $S$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi vortex.base -fo vortex -param S -h0 1 -scount 4 -maxcount 40 -mo vortexadapt
```

### Unsteady 3-D dynamics
1. Compute the $m=-1$ Hopf point along $S=1$
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi vortex.msh -fo vortexS1m1 -1/Re 0.0061 -S 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fo vortexS1m1 -fi vortexS1m1.base -sym -1 -eps_target 0.1+1.1i -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fo vortexS1m1 -fi vortexS1m1.mode -nf 0 -param 1/Re
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fo vortexS1m1 -fi vortexS1m1.mode -adaptto bda -mo vortexS1m1 -param 1/Re
```

2. Compute base state near the double Hopf point
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi vortex.msh -fo vortexDH -1/Re 0.0139 -S 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi vortexDH.base -fo vortexDH -S 1.44
```

3. Compute near-critical modes
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fo vortexm1 -fi vortexDH.base -sym -1 -eps_target 0+1i -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fo vortexm2 -fi vortexDH.base -sym -2 -eps_target 0+2i -eps_pos_gen_non_hermitian
```

4. Compute Hopf-Hopf point assuming non-resonant interaction
```sh
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fo vortexDH -fi vortexm2.mode -fi2 vortexm1.mode -param S -param2 1/Re -nf 0
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fo vortexDHadapt -fi vortexDH.hoho -param S -param2 1/Re -adaptto bda -mo vortexm1m2adapt
```

5. Compute Hopf-Hopf point assuming $2:1$ resonant interaction
```sh
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fo vortexDHadapt21res -fi vortexDHadapt.hoho -param S -param2 1/Re -res1x 2
```
NOTE: the resonant coefficients $\gamma_{12}$ and $\gamma_{22}$ differ from the original paper by unit-magnitude phase offset due to a different phase reference condition used by `ff-bifbox`.