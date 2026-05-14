# "Turbulent" Swirling Jet Example: Chevalier etal, TCFD, (2024)
This file shows an example `ff-bifbox` workflow for reproducing the results of the study:
```bibtex
@article{chevalier_etal_2024,
  title={Resolvent analysis of a swirling turbulent jet},
  volume={38},
  DOI={10.1007/s00162-024-00704-2},
  journal={Theoretical and Computational Fluid Dynamics},
  publisher={Springer},
  author={Chevalier, Quentin and Douglas, Christopher M. and Lesshafft, Lutz},
  year={2024},
  pages={641-663}
}
```
The commands below illustrate how to perform a mean flow resolvent analysis of an incompressible swirling jet with modeled turbulence (based on the Spalart–Allmaras turbulence model [SA-noft2](https://tmbwg.github.io/turbmodels/spalart.html#sanoft2)) using `ff-bifbox`. The original codes used for the paper can be found on Quentin Chevalier's repository at [github.com/hawkspar](https://github.com/hawkspar).

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial \tilde{\nu}}{\partial t} + u_j\frac{\partial \tilde{\nu}}{\partial x_j} - c_{b1}\tilde{S}\tilde{\nu} + c_{w1}f_w\frac{\tilde{\nu}^2}{d^2} - \frac{1}{\sigma}\left\lbrace\frac{\partial }{\partial x_i}\left[\left(\frac{1}{Re} + \tilde{\nu}\right)\frac{\partial \tilde{\nu}}{\partial x_i}\right] + c_{b2}\frac{\partial \tilde{\nu}}{\partial x_i}\frac{\partial \tilde{\nu}}{\partial x_i}\right\rbrace &= 0 \\
\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{\partial}{\partial x_j}\left[\left(\frac{1}{Re} + \nu_t\right)\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)\right] &= 0 \\
\frac{\partial u_i}{\partial x_i} &= 0 
\end{align*}
$$

where:
- $`\nu_t=\tilde{\nu}f_{v1}`$
- $`f_{v1}=\frac{\tilde{\nu}^3}{\tilde{\nu}^3+\left(\frac{c_{v1}}{Re}\right)^3}`$
- $`\tilde{S}=\sqrt{\frac{1}{2}\left(\frac{\partial u_i}{\partial x_j}-\frac{\partial u_j}{\partial x_i}\right)\left(\frac{\partial u_i}{\partial x_j}-\frac{\partial u_j}{\partial x_i}\right)}+\frac{\tilde{\nu}}{\kappa^2 d^2}f_{v2}`$
- $`f_{v2}=1-\frac{\tilde{\nu}}{\frac{1}{Re}+\tilde{\nu}f_{v1}}`$
- $`f_w=g\left[\frac{1+c_{w3}^6}{g^6+c_{w3}^6}\right]^{1/6}`$
- $`g=\left[1+c_{w2}\left(\min\left(\frac{\tilde{\nu}}{\tilde{S}\kappa^2d^2},10\right)^5-1\right)\right]\min\left(\frac{\tilde{\nu}}{\tilde{S}\kappa^2d^2},10\right)`$

together with the boundary conditions:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $`u_x=\tanh\left[6\left(1-r^2\right)\right]`$, $u_r=0$, $u_{\theta}=Sr\tanh\left[6\left(1-r^2\right)\right]$, $\tilde{\nu}=10^{-6}$ |
| Wall, $\Gamma_w$ | $u_x=u_r=u_{\theta}=\tilde{\nu}=0$ |
| Co-flow, $\Gamma_c$ | $u_x=2\alpha\left(\frac{r-1-\epsilon}{19-\epsilon}\right)\left(1-\frac{r-1-\epsilon}{2\left(19-\epsilon\right)}\right)$, $u_r=u_{\theta}=0$, $\tilde{\nu}=10^{-6}$ |
| Axis, $\Gamma_a$| $`\begin{cases}\frac{\partial u_x}{\partial r}=u_r=u_{\theta}=\frac{\partial \tilde{\nu}}{\partial r}=0, & \text{if } m=0 \\\\ u_x=\frac{\partial u_r}{\partial r}=\frac{\partial u_{\theta}}{\partial r}=\tilde{\nu}=0, & \text{if } \|m\|=1 \\\\ u_x=u_r=u_{\theta}=\tilde{\nu}=0, & \text{if } \|m\|>1\end{cases}`$ |
| Open, $\Gamma_o$ | $\left(\frac{1}{Re}+\nu_t\right)\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)\hat{n}_j-p\hat{n}_i = \frac{\partial\tilde{\nu}}{\partial x_i}\hat{n}_i=0$ |

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_c+\Gamma_a+\Gamma_o$. Solutions $\vec{q}=\left(u_i,\tilde{\nu},p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{\tilde{\nu}},\check{p}\right)^T$,

$$
\begin{align*} 
&\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\left(\frac{1}{Re}+\frac{\tilde{\nu}^4}{\tilde{\nu}^3+\left(\frac{c_{v1}}{Re}\right)^3}\right)\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)\right)_{\Omega} \\
&\left(\check{\tilde{\nu}},\frac{\partial \tilde{\nu}}{\partial t} + u_i\frac{\partial \tilde{\nu}}{\partial x_i} - c_{b1}\tilde{S}\tilde{\nu} + c_{w1}f_w\frac{\tilde{\nu}^2}{d^2} - \frac{c_{b2}}{\sigma}\frac{\partial \tilde{\nu}}{\partial x_i}\frac{\partial \tilde{\nu}}{\partial x_i}\right)_{\Omega} + \left(\frac{\partial \check{\tilde{\nu}}}{\partial x_i},\frac{1}{\sigma}\left(\frac{1}{Re}+\tilde{\nu}\right)\frac{\partial \tilde{\nu}}{\partial x_i}\right)_{\Omega} \\
&- \left(\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} = 0.
\end{align*}
$$

Similarly, for the frozen viscosity flow (used for linear analysis), the same formulation is used, but with $\tilde{\nu}$ removed as an unknown (and the $\check{\tilde{\nu}}$ terms eliminated).

These weak formulations have been implemented in the equations files for this example: [eqns_chevalier_etal_2024_baseflow.idp](./eqns_chevalier_etal_2024_baseflow.idp) and [eqns_chevalier_etal_2024_perturbations.idp](./eqns_chevalier_etal_2024_perturbations.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/chevalier_etal_2024/data
export nproc=4
```
3. Create symbolic links for solver settings.
```sh
ln -sf examples/chevalier_etal_2024/settings_chevalier_etal_2024.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### Build initial mesh directly from `.geo` files using Gmsh
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/chevalier_etal_2024 -dir $workdir -mi nozzle_lg.geo
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/chevalier_etal_2024 -dir $workdir -mi nozzle_sm.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).

## Perform parallel computations using `ff-bifbox`
### Steady axisymmetric dynamics of the mean flow with modeled turbulence

0. Select base flow equations to model turbulent eddy viscosity using the Spalart--Allmaras turbulence model. NOTE: There are good reasons to seriously doubt the choice of the SA model (and the Boussinesq hypothesis entirely) for this flow.
```sh
ln -sf examples/chevalier_etal_2024/eqns_chevalier_etal_2024_baseflow.idp eqns.idp
```

1. Compute base state on the large mesh at $Re=10$, $S=0$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi nozzle_lg.msh -fo S0p0Re10lg -1/Re 0.1 -S 0
```

2. Adapt base state to a coarser mesh for continuation
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p0Re10lg.base -fo jet_adapt_0 -mo nozzle_adapt_0 -err 0.1 -thetamax 0.01
```

3. Continue base state along the parameter $1/Re$. NOTE: This problem becomes very poorly scaled at high $Re$, meaning numerical issues may arise and require workarounds.
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi jet_adapt_0.base -fo jet_adapt -param 1/Re -h0 -1 -scount 4 -maxcount 100 -mo nozzle_adapt -err 0.1 -anisomax 3 -thetamax 0.01 -snes_max_it 50
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi jet_adapt_100.base -fo jet_adapt_101 -mi nozzle_lg.msh -1/Re 5.0e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi jet_adapt_101.base -fo jet_adapt_102 -1/Re 3.0e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi jet_adapt_102.base -fo jet_adapt_103 -1/Re 1.0e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi jet_adapt_103.base -fo jet_adapt_104 -1/Re 8.0e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi jet_adapt_104.base -fo jet_adapt_105 -1/Re 6.0e-6
```

4. Compute the $Re=200,000$, $S=0$ case on the reference mesh
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi jet_adapt_105.base -fo S0p0Re200000lg -1/Re 5.0e-6 -snes_linesearch_type l2 -pv 1 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
```

5. Continue base state along the parameter $S$ in increments of $0.1$ using zeroth-order continuation
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p0Re200000lg.base -fo S0p1Re200000lg -pv 1 -S 0.1 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p1Re200000lg.base -fo S0p2Re200000lg -pv 1 -S 0.2 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p2Re200000lg.base -fo S0p3Re200000lg -pv 1 -S 0.3 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p3Re200000lg.base -fo S0p4Re200000lg -pv 1 -S 0.4 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p4Re200000lg.base -fo S0p5Re200000lg -pv 1 -S 0.5 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p5Re200000lg.base -fo S0p6Re200000lg -pv 1 -S 0.6 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p6Re200000lg.base -fo S0p7Re200000lg -pv 1 -S 0.7 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p7Re200000lg.base -fo S0p8Re200000lg -pv 1 -S 0.8 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p8Re200000lg.base -fo S0p9Re200000lg -pv 1 -S 0.9 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi S0p9Re200000lg.base -fo S1p0Re200000lg -pv 1 -S 1.0 -snes_rtol 0 -snes_atol 1e-12 -snes_stol 0
```

### Dynamics of coherent perturbations to turbulent mean state

0. Select perturbation equations to use a frozen eddy viscosity model
```sh
ln -sf examples/chevalier_etal_2024/eqns_chevalier_etal_2024_perturbations.idp eqns.idp
```

1. Compute eigenvalue spectrum of the $Re=200,000$, $S=1$ flow on the smaller mesh for different azimuthal wavenumbers (to confirm that the flow is globally stable).
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi S1p0Re200000lg.base -so S1p0Re200000 -eps_target 0.1+0.25i -ntarget 13 -targetf 0.1+6.25i -sym 0 -eps_nev 20 -eps_pos_gen_non_hermitian -mi nozzle_sm.msh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi S1p0Re200000lg.base -so S1p0Re200000 -eps_target 0.1-6.25i -ntarget 26 -targetf 0.1+6.25i -sym -1 -eps_nev 20 -eps_pos_gen_non_hermitian -mi nozzle_sm.msh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi S1p0Re200000lg.base -so S1p0Re200000 -eps_target 0.1-6.25i -ntarget 26 -targetf 0.1+6.25i -sym -2 -eps_nev 20 -eps_pos_gen_non_hermitian -mi nozzle_sm.msh
```

2. Compute dominant resolvent gain of the $Re=200,000$, $S=1$ flow on the smaller mesh for different azimuthal wavenumbers. NOTE: Here, forcing is NOT restricted to a set portion of the mesh as it is in the paper, which may lead to numerical artifacts at low frequencies due to spurious variations in the outer coflow surrounding the jet.
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi S1p0Re200000lg.base -so S1p0Re200000 -omega 0 -nomega 64 -omegaf 6.3 -sym 0 -eps_nev 1 -strict 1 -mi nozzle_sm.msh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi S1p0Re200000lg.base -so S1p0Re200000 -omega -6.3 -nomega 127 -omegaf 6.3 -sym -1 -eps_nev 1 -strict 1 -mi nozzle_sm.msh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi S1p0Re200000lg.base -so S1p0Re200000 -omega -6.3 -nomega 127 -omegaf 6.3 -sym -2 -eps_nev 1 -strict 1 -mi nozzle_sm.msh
```

3. Compute dominant resolvent forcing and response modes for $Re=200,000$, $S=1$, $St=0.004$, and $m=\pm2$ as in the paper in Figs 8,9,10. NOTE: $\omega=2\pi{}St$, and, due to the SO(2) symmetry, $(S,m,St)=(S,-m,-St)$. Again, forcing is NOT restricted to the jet region, which can lead to significant differences due to spurious variations in the outer coflow surrounding the jet.
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi S1p0Re200000lg.base -mo S1p0Re200000St0p004 -omega -0.02513274 -nomega 2 -omegaf 0.02513274 -sym -2 -eps_nev 1 -strict 1 -mi nozzle_sm.msh -pv 1
```