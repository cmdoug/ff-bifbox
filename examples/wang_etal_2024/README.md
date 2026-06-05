# Low-Mach V-flame Example: Wang et al., JFM, (2024)
This file shows an example `ff-bifbox` workflow for reproducing the results in the study:
```bibtex
@article{wang_etal_2024,
  title = {{Onset of global instability in a premixed annular V-flame}},
  author = {Wang, Chuhan and Douglas, Christopher M. and Guan, Yu and Xu, Chunxiao and Lesshafft, Lutz},
  journal={Journal of Fluid Mechanics},
  publisher={Cambridge University Press},
  volume={998},
  pages={A23},
  year={2024},
  DOI={10.1017/jfm.2024.869},
}
```
The commands below illustrate how to perform a bifurcation analysis of a lean premixed V-flame in an axisymmetric annular jet using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial\rho}{\partial t} + u_i\frac{\partial \rho}{\partial x_i} + \rho\frac{\partial u_i}{\partial x_i} &= 0 \\
\rho\frac{\partial u_i}{\partial t} + \rho u_j\frac{\partial u_i}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{\partial\tau_{ij}}{\partial x_j} &= 0 \\
\rho\frac{\partial Y_{\mathrm{CH}_4}}{\partial t} + \rho u_i\frac{\partial Y_{\mathrm{CH}_4}}{\partial x_i} - \frac{\partial}{\partial x_i}\left(\frac{\mu}{Sc}\frac{\partial Y_{\mathrm{CH}_4}}{\partial x_i}\right) + W_{\mathrm{CH}_4}\mathcal{Q}&= 0 \\
c_p\rho\frac{\partial T}{\partial t} + c_p\rho u_i\frac{\partial T}{\partial x_i} - \frac{\partial}{\partial x_i}\left(\frac{c_p\mu}{Pr}\frac{\partial T}{\partial x_i}\right) + \Delta h_f^0\mathcal{Q} &= 0
\end{align*}
$$

where:
- $`\tau_{ij}=\mu\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}-\frac{2}{3}\delta_{ij}\frac{\partial u_k}{\partial x_k}\right)`$
- $`p_0 = R_s\rho T`$
- $`\mu=\frac{1}{sg}\frac{A_sT^{1/2}}{1+T_s/T}`$
- $`sg\left(r,x\right)=1 \text{ if } r \leq r_{sg} \text{ and } x \leq x_{sg}`$
- $`sg\left(r,x\right)=1 + \left(\alpha-1\right)\zeta\left(x,x_{sg}\right) \text{ if } r \leq r_{sg} \text{ and } x > x_{sg}`$
- $`sg\left(r,x\right)=sg\left(r_{sg}, x\right) + \left[\alpha-sg\left(r_{sg},x\right)\right]\zeta\left(r,r_{sg}\right) \text{ if } r > r_{sg}`$
- $`\zeta\left(a,b\right)=\frac{1}{2}+\frac{1}{2}\tanh\left\lbrace\tan\left(-\frac{\pi}{2}+\pi\frac{|a-b|}{l}\right)\right\rbrace`$
- $`\mathcal{Q}=A_r\left(\rho\frac{Y_{\mathrm{CH}_4}}{W_{\mathrm{CH}_4}}\right)^{n_{\mathrm{CH}_4}}\left(\rho\frac{Y_{\mathrm{O}_2}}{W_{\mathrm{O}_2}}\right)^{n_{\mathrm{O}_2}}\exp\left(-\frac{T_a}{T}\right)`$.

The boundary conditions are:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $`u_x=93750(2r - 0.003)(0.011 - 2r)U_0\text{ m/s}`$, $`u_r=0`$, $`Y_{\mathrm{CH}_4}=0.04256`$, $`T=300\text{ K}`$ |
| Annular channel wall, $\Gamma_w$ | $`u_x=u_r=\frac{\partial Y_{\mathrm{CH}_4}}{\partial r}=\frac{\partial T}{\partial r}=0`$ |
| Dump plane wall, $\Gamma_{dp}$ | $`u_x=u_r= \frac{\partial Y_{\mathrm{CH}_4}}{\partial x}=0$, $T=300\text{ K}`$ |
| Centrebody wall, $\Gamma_{cb}$ | $`u_x=u_r= \frac{\partial Y_{\mathrm{CH}_4}}{\partial x}=0`$, $`T=T_{cb}\text{ K}`$ |
| Axis, $\Gamma_a$| $`frac{\partial u_x}{\partial r}=u_r=\frac{\partial Y_{\mathrm{CH}_4}}{\partial r}=\frac{\partial T}{\partial r}=0`$ |
| Lateral, $\Gamma_l$ | $`\frac{\partial u_x}{\partial r}=u_r=Y_{\mathrm{CH}_4}=0`$, $`T=300\text{ K}`$ |
| Outlet, $\Gamma_o$ | $`\tau_{ix}-p\hat{e}_x = \frac{\partial Y_{\mathrm{CH}_4}}{\partial x} = \frac{\partial T}{\partial x} = 0`$ |

The present implementation is based on a weak formulation of these equations. Density is eliminated as an unknown using the equation of state, test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_{dp}+\Gamma_{cb}+\Gamma_a+\Gamma_l+\Gamma_o$. Solutions $\vec{q}=\left(u_i,Y_{\mathrm{CH}_4},T,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{Y},\check{T},\check{p}\right)^T$,

$$
\begin{align*}
&\left(\check{u}_i,\rho\left[\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right]\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\mu\left[\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}-\frac{2}{3}\delta_{ij}\frac{\partial u_k}{\partial x_k}\right]\right)_{\Omega} \\
&+\left(\check{Y},\rho\left[\frac{\partial Y_{\mathrm{CH}_4}}{\partial t} + u_i\frac{\partial Y_{\mathrm{CH}_4}}{\partial x_i}\right] + W_{\mathrm{CH}_4}\mathcal{Q}\right)_{\Omega} + \left(\frac{\partial \check{Y}}{\partial x_i},\frac{\mu}{Sc}\frac{\partial Y_{\mathrm{CH}_4}}{\partial x_i}\right)_{\Omega} \\
&+\left(\check{T},\rho\left[\frac{\partial T}{\partial t} + u_i\frac{\partial T}{\partial x_i}\right] + \frac{\Delta h_f^0}{c_p}\mathcal{Q}\right)_{\Omega} + \left(\frac{\partial \check{T}}{\partial x_i},\frac{\mu}{Pr}\frac{\partial T}{\partial x_i}\right)_{\Omega} \\
&+\left(\check{p},\frac{1}{T}\frac{\partial T}{\partial t}+\frac{u_i}{T}\frac{\partial T}{\partial x_i}-\frac{\partial u_i}{\partial x_i}\right)_{\Omega} = 0
\end{align*}
$$

This weak formulation has been implemented in the equations file for this example: [eqns_wang_etal_2024.idp](./eqns_wang_etal_2024.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/wang_etal_2024/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/wang_etal_2024/eqns_wang_etal_2024.idp eqns.idp
ln -sf examples/wang_etal_2024/settings_wang_etal_2024.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from `.geo` files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/wang_etal_2024 -dir $workdir -mi Vflame.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/wang_etal_2024/Vflame.md -mo $workdir/Vflame
```

## Perform parallel computations using `ff-bifbox`
Note that, unlike most of the `ff-bifbox` examples, the results in `wang_etal_2024` are reported in dimensional variables. While non-dimensionalization is preferable for better numerical scaling, the SI units used in the study are retained here for the sake of example.
### Laminar base flow
1. Compute a non-reacting base state with reference parameters on the initial mesh.
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir examples/wang_etal_2024/data -mi Vflame.msh -fo nonreacting_0 -U0 0.1 -Tcb 700 -As 1.67212e-6 -Ts 170.672 -Pr 0.7 -Sc 0.7 -p0 101325 -Rs 264.56013215560904 -Cp 1.3 -YCH4 0.04256 -WCH4 0.016 -YO2 0.2128 -WO2 0.032 -nCH4 1.0 -nO2 0.5 -Ar 0 -Ta 10065.425264217412 -Dh0f -804.084 -alpha 0.05 -xsg 0.15 -rsg 0.03 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting_0.base -fo nonreacting_1 -U0 0.5 -mo nonreacting_1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting_1.base -fo nonreacting -U0 2.2 -pv 1 -mo nonreacting
```

2. Turn on chemistry and ignite the $U_0=2.2$ m/s flow at an elevated centrebody temperature and lower combustion enthalpy. Then perform continuation back to reference parameters. Coarse meshes are used for computational efficiency and stabilizing artificial dissipation. 
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting.base -fo ignite_0 -Tcb 1000 -Ar 1.1e7 -Dh0f -100 -mo ignite_0 -snes_rtol 0 -err 0.05
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi ignite_0.base -fo ignite -param Dh0f -h0 -200 -mo ignite -dmax 100 -err 0.1 -scount 5 -paramtarget -804.084 -maxcount -1 -contorder 2
cd $workdir && export lastfile=$(printf '%s\n' ignite_*.base | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo ignited -Tcb 700 -Dh0f -804.084 -mo ignited -snes_rtol 0 -err 0.05 -snes_linesearch_type l2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignited.base -fo U02p2 -pv 1 -snes_rtol 0 -snes_linesearch_type l2 -mo U02p2
```

3. Save base flows in 0.1 m/s increments for $U_0=2.3$ m/s to $U_0=3.8$ m/s. Dissipation from mesh coarsening is used to aid convergence at each step before refining the coarse solutions on the reference mesh.
```sh
for i in {3..9}
do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi "U02p$(($i-1))".base -fo U0inc -U0 2."$i" -snes_linesearch_type l2 -mo U0inc -err 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U0inc.base -fo U02p"$i" -pv 1 -snes_rtol 0 -snes_linesearch_type l2 -mo U02p"$i"
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U02p"$i".base -fo U02p"$i" -pv 1 -snes_rtol 0 -snes_linesearch_type l2 -mo U02p"$i"
done
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U02p9.base -fo U0inc -U0 3.0 -snes_linesearch_type l2 -mo U0inc -err 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U0inc.base -fo U03p0 -pv 1 -snes_rtol 0 -mo U03p0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U03p0.base -fo U03p0 -pv 1 -snes_rtol 0 -mo U03p0
for i in {1..8}
do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi "U03p$(($i-1))".base -fo U0inc -U0 3."$i" -snes_linesearch_type l2 -mo U0inc -err 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U0inc.base -fo U03p"$i" -pv 1 -snes_rtol 0 -snes_linesearch_type l2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U03p"$i".base -fo U03p"$i" -pv 1 -snes_rtol 0 -snes_linesearch_type l2 -mo U03p"$i"
done
```

### Global linear analysis
4. Compute global eigenspectra at $Re=1978$, $Re=2282$, $Re=2586$, and $Re=2891$. Notably, unlike in the paper, the present results do not identify any criticality of the leading flame-tip eigenmode at $Re=2815$. The source of this disagreement is not known. 
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi U02p6.base -so Re1978 -eps_nev 25 -eps_target 25+250i -ntarget 8 -targetf 50+2000i
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi U03p0.base -so Re2282 -eps_nev 25 -eps_target 25+250i -ntarget 8 -targetf 50+2000i
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi U03p4.base -so Re2586 -eps_nev 25 -eps_target 25+250i -ntarget 8 -targetf 50+2000i
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi U03p8.base -so Re2891 -eps_nev 25 -eps_target 25+250i -ntarget 8 -targetf 50+2000i
```
5. Compute leading flame-tip eigenmode at $Re=2282$ 
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi U03p0.base -fo Re2282flametip -eps_target 1+600i -pv 1
```

6. Compute optimal gain curve in velocity 2-norm at $Re=2282$ and $Re=2586$.
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi U03p0.base -so Re2282 -omega 20 -omegaf 2000 -nomega 100
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi U03p4.base -so Re2586 -omega 20 -omegaf 2000 -nomega 100
```

7. Compute optimal response at $Re=2586$ for $St=0.31$ and $St=0.68$.
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi U03p4.base -fo Re2586St0p31 -omega 602 -pv 1
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi U03p4.base -fo Re2586St0p68 -omega 1320 -pv 1
```

### Nonlinear analysis
8. Compute nonlinear dynamics at $Re=1978$ in time domain. Here, the maximum velocity magnitude of the eigenmode and the base flow are rescaled to provide initial velocity perturbation amplitudes that correspond to the 20% and 50% amplitudes used in the paper. Nonetheless, since the phase was not explicitly specified, variations in the initial phase may cause these results to differ qualitatively from those in the paper.
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi U02p6.base -fo Re1978flametip -eps_target 1+600i -pv 1
ff-mpirun -np 1 examples/wang_etal_2024/moderescale.md -v 0 -dir $workdir -fi Re1978flametip.mode -fo Re1978flametip_scaled
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir -bfi U02p6.base -fi Re1978flametip_scaled.mode -fo Re1978perturbation20perc -amp 0.2 -ts_time_step 0.0001 -ts_adapt none -scount 2 -maxcount 6000 -mo Re1978perturbation20perc
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir -bfi U02p6.base -fi Re1978flametip_scaled.mode -fo Re1978perturbation50perc -amp 0.5 -ts_time_step 0.0001 -ts_adapt none -scount 2 -maxcount 6000 -mo Re1978perturbation50perc
```