# Pool-fire puffing as a hydrodynamic global mode: Moreno-Boza et al. (Combustion and Flame, 2018 )
This example shows how to use `ff-bifbox` to compute the steady axisymmetric base flow and the leading global modes of the low-Mach-number round pool-fire problem considered by Moreno-Boza et al. (2018) in: 

```bibtex
@article{moreno-boza_etal_2018,
  title={On the critical conditions for pool-fire puffing},
  author={Moreno-Boza, Daniel and Coenen, Wilfried and Carpio, Jaime and S{\'a}nchez, Antonio L and Williams, Forman A},
  journal={Combustion and Flame},
  volume={192},
  pages={426--438},
  year={2018},
  publisher={Elsevier}
  doi={10.1016/j.combustflame.2018.02.011},
}
```

The model uses a Lagrange-multiplier formulation for the weak enforcement of the vaporizing-fuel boundary condition at the pool surface (and also the isothermal boundary conditions on the surrounding wall when enforced). The transported scalar variables are the weighted coupling function `Zt` ($\tilde{Z}$) and the excess-enthalpy variable `H` ($\xi$). The density and transport coefficients are reconstructed from the thermochemical closures defined in the equation file. This equation system has numerical challenges because the governing equations have a non-smooth dependence on $\tilde{Z}$. As such, it is possible for the solver to stagnate or diverge when successive iterations push the cusp back-and-forth across quadrature points. To overcome such issues, it is recommended to perturb the parameters and/or the mesh and try again.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial \rho}{\partial t} + u_i\frac{\partial \rho}{\partial x_i} + \rho\frac{\partial u_i}{\partial x_i} &= 0 \\
\frac{Ra}{Pr}\rho\left(\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right) + \frac{\partial p}{\partial x_i} - \left(1-\rho\right)\hat{e}_x - \frac{\partial}{\partial x_j}\left[\mu\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)\right] &= 0 \\
Ra\rho\left(\frac{\partial \xi}{\partial t} + u_j\frac{\partial \xi}{\partial x_j}\right) - \frac{\partial}{\partial x_i}\left(\rho D_T\frac{\partial \xi}{\partial x_i}\right) &= 0 \\
Ra\rho\left(\frac{\partial Z}{\partial t} + u_j\frac{\partial Z}{\partial x_j}\right) - \frac{1}{Le}\frac{\partial}{\partial x_i}\left(\rho D_T\frac{\partial \tilde{Z}}{\partial x_i}\right) &= 0
\end{align*}
$$

where:
- $\mu=\rho D_T=T^\sigma$
- $T=1+\left(T_B-1-\frac{q}{S}\right)\xi+\frac{q}{S}\min\left(1,\frac{\tilde{Z}}{\tilde{Z}_S}\right)$
- $\rho = \frac{1}{T\left[1+\max\left(0, \frac{\tilde{Z}-\tilde{Z}_S}{1-\tilde{Z}_S}\right)\left(\frac{W_A}{W_F}-1\right)\right]}$
- $Z=Z_S\min\left(1,\frac{\tilde{Z}}{\tilde{Z}_S}\right) + \max\left(0, \frac{\tilde{Z} - \tilde{Z}_S}{1 - \tilde{Z}_S}\right)\left(1 - Z_S\right)$.

The boundary conditions are:

| Boundary | Constraints |
| :--- | :--- |
| Pool, $\Gamma_p$ | $-T_B^\sigma\frac{\partial\xi}{\partial x} = \alpha Ra \rho u_x$, $u_r=0$, $\xi=1$, $-T_B^\sigma\frac{\partial\tilde{Z}}{\partial x} =Le_F Ra \rho u_x\left(1-\tilde{Z}\right)$ |
| Wall, $\Gamma_w$ | $`\begin{cases}u_x=u_r=\frac{\partial\tilde{Z}}{\partial x}=\frac{\partial\xi}{\partial x}=0, & \text{if adiabatic} \\\\ u_x=u_r=\frac{\partial\tilde{Z}}{\partial x}=0, \xi=\frac{q/S}{q/S+1-T_B}\min\left(1,\frac{\tilde{Z}}{\tilde{Z}_S}\right), & \text{if isothermal}\end{cases}`$ |
| Axis, $\Gamma_a$| $\frac{\partial u_x}{\partial r}=u_r=\frac{\partial \xi}{\partial r}=\frac{\partial \tilde{Z}}{\partial r}=0$ |
| Lateral, $\Gamma_l$ | $\mu\left(\frac{\partial u_i}{\partial r}+\frac{\partial u_r}{\partial x_i}\right)-p\hat{e}_r = \xi = \tilde{Z} = 0$ |
| Outlet, $\Gamma_o$ | $\mu\left(\frac{\partial u_i}{\partial x}+\frac{\partial u_x}{\partial x_i}\right)-p\hat{e}_x = \frac{\partial \xi}{\partial x} = \frac{\partial \tilde{Z}}{\partial x}=0$ |

The present implementation is based on a weak formulation of these equations. Density is eliminated as an unknown algebraically using the equation of state, test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_p+\Gamma_w+\Gamma_a+\Gamma_l+\Gamma_o$. Solutions $\vec{q}=\left(u_i,p,\tilde{Z},\xi,\lambda\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{p},\check{\tilde{Z}},\check{\xi},\check{\lambda}\right)^T$,

$$
\begin{align*}
&\left(\check{u}_i,\frac{Ra}{Pr}\rho\left[\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right]\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},T^\sigma\left[\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right]\right)_{\Omega} \\
&+\left(\check{u}_x,\rho-1\right)_{\Omega} - \left(\check{u}_x,\lambda\right)_{\Gamma_i} + \left(\check{p}n_i,\rho u_i\right)_{\partial\Omega}-\left(\frac{\partial\check{p}}{\partial x_i},\rho u_i\right)_{\Omega} \\
&+\left(\check{\tilde{Z}},Ra\rho\frac{\partial Z}{\partial \tilde{Z}}\left[\frac{\partial \tilde{Z}}{\partial t} + u_i\frac{\partial \tilde{Z}}{\partial x_i}\right]\right)_{\Omega} + \left(\frac{\partial \check{\tilde{Z}}}{\partial x_i},\frac{T^\sigma}{Le}\frac{\partial \tilde{Z}}{\partial x_i}\right)_{\Omega} \\
&- \left(\check{\tilde{Z}},\frac{T_B^\sigma Le_F}{\alpha Le}\frac{\partial\xi}{\partial x}\left(\tilde{Z}-1\right)\right)_{\Gamma_i} + \left(\check{\xi},Ra\rho\left[\frac{\partial \xi}{\partial t} + u_i\frac{\partial \xi}{\partial x_i}\right]\right)_{\Omega} + \left(\frac{\partial \check{\xi}}{\partial x_i},T^\sigma\frac{\partial \xi}{\partial x_i}\right)_{\Omega} \\
&-\left(\check{\lambda},\frac{T_B^\sigma}{\alpha}\frac{\partial\xi}{\partial x} + Ra\rho u_x\right)_{\Gamma_i} + \underbrace{\left(\check{\lambda},T-1\right)_{\Gamma_w} + \left(\check{\xi},\lambda\right)_{\Gamma_w}}_{\text{removed for adiabatic wall case}} = 0.
\end{align*}
$$

This weak formulation has been implemented in the equations file for this example: [eqns_moreno-boza_etal_2018.idp](./eqns_moreno-boza_etal_2018.idp).

The workflow below targets the methanol parameter set and reproduces the base flow and global spectrum near the onset of the axisymmetric puffing instability using `ff-bifbox`.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/moreno-boza_etal_2018/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/moreno-boza_etal_2018/eqns_moreno-boza_etal_2018.idp eqns.idp
ln -sf examples/moreno-boza_etal_2018/settings_moreno-boza_etal_2018.idp settings.idp
```

## Build initial mesh
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### Build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/moreno-boza_etal_2018/puffing_geom.md -mo $workdir/jet
```

## Perform parallel computations for methanol using `ff-bifbox`
0. Select whether to run the adiabatic or isothermal wall BC case.
```sh
export wallBC=adiabatic
```
OR
```sh
export wallBC=isothermal
```

1. Compute a series of base states for gradually increasing $q/S$ (defined in the input files as `qoS`) and $Ra$ (defined as `Ra`). The target value used here is $Ra = 2 \times 10^4$, which is close to the onset of the axisymmetric puffing instability for the methanol parameter set.
```sh
[[ "$wallBC" == "isothermal" ]] && export isotherm=1 || export isotherm=0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi jet.msh -fo methanol_"$wallBC"_ignite_0 -pv 1 -Ra 1 -Pr 0.70 -sigma 0.70 -qoS 0.2 -S 6.47 -TB 1.12 -lv 2.84 -LeF 1.15 -WFoWA 1.1 -isotherm $isotherm
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_0.base -fo methanol_"$wallBC"_ignite_1 -pv 1 -Ra 5 -qoS 0.4 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_1.base -fo methanol_"$wallBC"_ignite_2 -pv 1 -Ra 10 -qoS 0.6 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_2.base -fo methanol_"$wallBC"_ignite_3 -pv 1 -Ra 30 -qoS 0.7 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_3.base -fo methanol_"$wallBC"_ignite_4 -pv 1 -Ra 100 -qoS 1.2 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_4.base -fo methanol_"$wallBC"_ignite_5 -pv 1 -Ra 200 -qoS 2.5 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_5.base -fo methanol_"$wallBC"_ignite_6 -pv 1 -Ra 400 -qoS 4 -snes_linesearch_type l2 -snes_atol 1e-4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_6.base -fo methanol_"$wallBC"_ignite_7 -pv 1 -Ra 500 -qoS 6 -snes_linesearch_type l2 -snes_atol 1e-4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_7.base -fo methanol_"$wallBC"_ignite_8 -pv 1 -Ra 600 -qoS 7 -snes_linesearch_type l2 -snes_atol 1e-4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_8.base -fo methanol_"$wallBC"_ignite_9 -pv 1 -Ra 800 -qoS 7.7 -snes_linesearch_type l2 -snes_atol 1e-4
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignite_9.base -count 9 -fo methanol_"$wallBC"_ignite -pv 1 -paramtarget 20000 -param Ra -snes_linesearch_type l2 -snes_atol 1e-3 -contorder 0 -kmax 10 -dmax 10 -maxcount -1 -h0 100
cd $workdir && export lastfile=$(printf '%s\n' methanol_"$wallBC"_ignite_*.base | sort -t_ -k4,4n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo methanol_"$wallBC"_ignited -pv 1 -Ra 20000 -snes_linesearch_type l2 -snes_atol 1e-3
```

2. Recompute the final base state on successively adapted meshes. Note that sometimes the solver stagnates here and so re-attempts or modifications to "jitter" the mesh and knock the quadrature points away from the cusp may be necessary.
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_ignited.base -fo methanol_"$wallBC"_adapt_0 -pv 1 -mo adapt_0 -snes_linesearch_type l2 -snes_atol 1e-3 -anisomax 1e6 -hmax 1.0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_adapt_0.base -fo methanol_"$wallBC"_adapt_1 -pv 1 -mo methanol_"$wallBC"_adapt_1 -snes_linesearch_type l2 -snes_atol 1e-3 -anisomax 1e6 -hmax 1.0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_adapt_1.base -fo methanol_"$wallBC"_adapt_2 -pv 1 -mo methanol_"$wallBC"_adapt_2 -snes_linesearch_type l2 -snes_atol 1e-3 -anisomax 1e6 -hmax 1.0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_adapt_2.base -fo methanol_"$wallBC"_adapt_3 -pv 1 -mo methanol_"$wallBC"_adapt_3 -snes_linesearch_type l2 -snes_atol 1e-4 -anisomax 1e6 -hmax 1.0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_adapt_3.base -fo methanol_"$wallBC"_adapt_4 -pv 1 -mo methanol_"$wallBC"_adapt_4 -snes_linesearch_type l2 -snes_atol 1e-5 -anisomax 1e6 -hmax 1.0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_adapt_4.base -fo methanol_"$wallBC"_adapt_5 -pv 1 -mo methanol_"$wallBC"_adapt_5 -snes_linesearch_type l2 -snes_atol 1e-6 -anisomax 1e6 -hmax 1.0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_adapt_5.base -fo methanol_"$wallBC"_Ra20000 -pv 1 -mo methanol_"$wallBC"_Re20000 -snes_linesearch_type l2 -anisomax 2 -hmax 1.0
```

3. Compute the leading direct and adjoint eigenmodes of the methanol base flow. Adapt the mesh to these modes, then recompute it and the global eigenvalue spectrum on the adapted mesh. With the eigenvalue convention used by `modecompute.md`, the leading mode has angular frequency approximately $\omega \approx 0.0112$ and $\omega \approx 0.0119$ for the adiabatic and isothermal cases, respectively. In the scaling of Moreno-Boza et al. (2018), these both correspond to $St = \omega/\pi \approx 3.4 \times 10^{-3}$, in agreement with their Fig 2(a).
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_Ra20000.base -fo methanol_"$wallBC"_Ra20000 -eps_target 0.01+0.01i -eps_nev 1 -pv 1 -eps_two_sided 1 -strict 1
ff-mpirun -np $nproc meshcompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_Ra20000.base -fi2 methanol_"$wallBC"_Ra20000.mode -fi3 methanol_"$wallBC"_Ra20000adj.mode -anisomax 2 -hmax 1.0 -mo methanol_"$wallBC"_Re20000_adapt
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi methanol_"$wallBC"_Re20000_adapt.msh -fi methanol_"$wallBC"_Ra20000.base -fo methanol_"$wallBC"_Ra20000_adapt -pv 1 -snes_linesearch_type l2
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_Ra20000_adapt.base -fo methanol_"$wallBC"_Ra20000_adapt -eps_target 0.01+0.01i -eps_nev 1 -pv 1 -eps_two_sided 1 -strict 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi methanol_"$wallBC"_Ra20000_adapt.base -so methanol_"$wallBC"_Ra20000_adapt -eps_target 0.01+0.01i -targetf 0.01+0.1i -ntarget 10 -eps_nev 10 -pv 1
```

## Perform parallel computations for heptane using `ff-bifbox`
1. The above process can be repeated for the case of heptane fuel by choosing the correct parameters as indicated in the paper. However, the heptane case has stronger discontinuities, making the numerical scheme less robust. An ignition routine for heptane is given below. This may be followed by successive mesh adaptation and stability analysis similar to the approach given above.
```sh
[[ "$wallBC" == "isothermal" ]] && export isotherm=1 || export isotherm=0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi jet.msh -fo heptane_"$wallBC"_ignite_0 -pv 1 -Ra 1 -Pr 0.70 -sigma 0.70 -qoS 0.2 -S 15.2 -TB 1.24 -lv 1.14 -LeF 1.8 -WFoWA 3.45 -isotherm $isotherm
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_0.base -fo heptane_"$wallBC"_ignite_1 -pv 1 -Ra 5 -qoS 0.4 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_1.base -fo heptane_"$wallBC"_ignite_2 -pv 1 -Ra 10 -qoS 0.6 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_2.base -fo heptane_"$wallBC"_ignite_3 -pv 1 -Ra 30 -qoS 0.7 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_3.base -fo heptane_"$wallBC"_ignite_4 -pv 1 -Ra 100 -qoS 1 -snes_linesearch_type l2 -snes_atol 1e-6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_4.base -fo heptane_"$wallBC"_ignite_5 -pv 1 -Ra 200 -qoS 1.5 -snes_linesearch_type l2 -snes_atol 1e-4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_5.base -fo heptane_"$wallBC"_ignite_6 -pv 1 -Ra 400 -qoS 2.3 -snes_linesearch_type l2 -snes_atol 1e-4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_6.base -fo heptane_"$wallBC"_ignite_7 -pv 1 -Ra 500 -qoS 3.6 -snes_linesearch_type l2 -snes_atol 3e-4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_7.base -fo heptane_"$wallBC"_ignite_8 -pv 1 -Ra 600 -qoS 5 -snes_linesearch_type l2 -snes_atol 1e-3
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_8.base -fo heptane_"$wallBC"_ignite_9 -pv 1 -Ra 800 -qoS 5.5 -snes_linesearch_type l2 -snes_atol 1e-3
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi heptane_"$wallBC"_ignite_9.base -count 9 -fo heptane_"$wallBC"_ignite -pv 1 -paramtarget 20000 -param Ra -snes_linesearch_type l2 -snes_atol 1e-3 -contorder 0 -kmax 10 -dmax 10 -maxcount -1 -h0 100
cd $workdir && export lastfile=$(printf '%s\n' heptane_"$wallBC"_ignite_*.base | sort -t_ -k4,4n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo heptane_"$wallBC"_ignited -pv 1 -Ra 20000 -snes_linesearch_type l2 -snes_atol 1e-3
```