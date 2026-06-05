# 2D Reacting Compressible Flow Example: Brokof et al, PROCI. (2024)
This file shows an example `ff-bifbox` workflow for reproducing the results in the study:
```bibtex
@article{brokof_etal_2024,
  title = {The role of hydrodynamic shear in the thermoacoustic response of slit flames},
  journal = {Proceedings of the Combustion Institute},
  volume = {40},
  number = {1},
  pages = {105362},
  year = {2024},
  doi = {10.1016/j.proci.2024.105362},
  author = {Brokof, Philipp and Douglas, Christopher M. and Polifke, Wolfgang},
}
```
The commands below illustrate how to analyze a 2D reacting compressible flow through a duct using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial\rho}{\partial t} + \frac{\partial \left(\rho u_i\right)}{\partial x_i} &= 0 \\
\frac{\partial \left(\rho u_i\right)}{\partial t} + \frac{\partial \left(\rho u_iu_j\right)}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{1}{Re}\frac{\partial}{\partial x_j}\epsilon_{ij} &= 0 \\
\frac{\partial \left(\rho Y\right)}{\partial t} + \frac{\partial \left(\rho u_i Y\right)}{\partial x_i} - \frac{1}{Pe Le}\frac{\partial^2 Y}{\partial x_i^2} + \dot{\omega}_Y &= 0 \\
\frac{\partial \left(\rho h_s-p\gamma Ma^2\right)}{\partial t} + \frac{\partial \left(\rho u_i h_s\right)}{\partial x_i} - \gamma Ma^2\frac{\partial p}{\partial x_i}u_i \\
- \frac{\gamma \Delta T}{\gamma-1}\dot{\omega}_Y - \frac{\gamma}{\gamma-1}\frac{1}{Pe}\frac{\partial^2 T}{\partial x_i^2} - \frac{\gamma Ma^2}{Re} \frac{\partial u_i}{\partial x_j}\epsilon_{ij} &= 0
\end{align*}
$$

where:
- $\epsilon_{ij}=\mu\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}-\delta_{ij}\frac{\partial u_k}{\partial x_k}\right)$
- $\rho = \frac{1+\gamma Ma^2 p}{T}$
- $\dot{\omega}_Y=Da\rho Y \exp\left(\frac{-\left(1+\Delta T\right)^2Ze}{\Delta T T}\right)$
- $h_s = \frac{\gamma}{\gamma-1}T$

The boundary conditions are:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $\bar{u}_x=\frac{2}{5}$, $\bar{u}_y=0$, $\bar{Y}=\bar{T}=1$ for base flow|
| Inlet, $\Gamma_i$ | $\bar{\rho}\acute{u}_x+\acute{\rho}\bar{u}_x=\acute{I}_{\text{mass},i}$, $2\bar{\rho}\bar{u}_x\acute{u}_x+\acute{\rho}\bar{u}_x^2=\acute{I}_{x\text{mom},i}$, $\acute{u}_y=\acute{Y}=\acute{T}=0$ for perturbation|
| Wall, $\Gamma_w$ | $u_x=u_y=\frac{\partial Y}{\partial x_i}\hat{n}_i=\frac{\partial T}{\partial x_i}\hat{n}_i=0$ |
| Symmetry, $\Gamma_s$ | $\frac{\partial u_x}{\partial y}=u_y=\frac{\partial Y}{\partial y}=\frac{\partial T}{\partial y}=0$ |
| Outlet, $\Gamma_o$ | $\bar{\epsilon}_{ix}\hat{e}_x-\bar{p}\hat{e}_x = \frac{\partial \bar{Y}}{\partial x} = \frac{\partial \bar{T}}{\partial x} = 0$ for base flow |
| Outlet, $\Gamma_o$ | $\bar{\rho}\acute{u}_x+\acute{\rho}\bar{u}_x=\acute{I}_{\text{mass},o}$, $2\bar{\rho}\bar{u}_x\acute{u}_x+\acute{\rho}\bar{u}_x^2=\acute{I}_{x\text{mom},o}$, $\acute{\epsilon}_{yx}\hat{e}_x = \frac{\partial \acute{Y}}{\partial x} = \frac{\partial \acute{T}}{\partial x}=0$ for perturbation |

where:
- $\acute{I}_{\text{mass},i}=\frac{1}{2}\left(\acute{\rho}-\frac{\bar{\rho}}{\bar{c}}\acute{u}_x\right)\left[\left(\bar{c}+\bar{u}_x\right)R_{in}+\left(\bar{u}_x - \bar{c}\right)\right]$
- $\acute{I}_{x\text{mom},i}=\frac{1}{2}\left(\acute{\rho}-\frac{\bar{\rho}}{\bar{c}}\acute{u}_x\right)\left[\left(\bar{c}+\bar{u}_x\right)^2R_{in}+\left(\bar{u}_x - \bar{c}\right)^2\right]$
- $\acute{I}_{\text{mass},o}=\frac{1}{2}\left(\acute{\rho}+\frac{\bar{\rho}}{\bar{c}}\acute{u}_x\right)\left[\left(\bar{c}+\bar{u}_x\right)+\left(\bar{u}_x - \bar{c}\right)R_{out}\right]$
- $\acute{I}_{x\text{mom},o}=\frac{1}{2}\left(\acute{\rho}+\frac{\bar{\rho}}{\bar{c}}\acute{u}_x\right)\left[\left(\bar{c}+\bar{u}_x\right)^2+\left(\bar{u}_x - \bar{c}\right)^2R_{out}\right]$
- $c=\frac{\sqrt{\bar{T}}}{Ma}$

The present implementation is based on a weak formulation of these equations. Density is eliminated as an unknown using the equation of state, test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_s+\Gamma_o$. Solutions $\vec{q}=\left(u_i,Y,T,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{Y},\check{T},\check{p}\right)^T$,

$$
\begin{align*}
&\left(\check{u}_i,\frac{\partial \left(\rho u_i\right)}{\partial t}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},-\rho u_i u_j+\frac{1}{Re}\epsilon_{ij}\right)_{\Omega} \\
&+\left(\check{Y},\frac{\partial \left(\rho Y\right)}{\partial t} + \dot{\omega}_Y\right)_{\Omega} + \left(\check{Y}\hat{n}_i,\rho u_i Y\right)_{\partial\Omega} + \left(\frac{\partial \check{Y}}{\partial x_i},-\rho u_i Y + \frac{1}{PeLe}\frac{\partial Y}{\partial x_i}\right)_{\Omega} \\
&+\left(\check{T},Ma^2\left[\frac{\partial p}{\partial t} + u_i\frac{\partial p}{\partial x_i}\right] + \left(1+\gamma Ma^2p\right)\frac{\partial u_i }{\partial x_i} - \Delta T\dot{\omega}_Y - \frac{\left(\gamma-1\right) Ma^2}{Re} \frac{\partial u_i}{\partial x_j}\epsilon_{ij}\right)_{\Omega} \\
&+ \left(\frac{\partial \check{T}}{\partial x_i},\frac{1}{Pe}\frac{\partial T}{\partial x_i}\right)_{\Omega} + \left(\check{p}, \frac{\partial \rho}{\partial t}\right)_{\Omega} - \left(\frac{\partial\check{p}}{\partial x_i}, \rho u_i\right)_{\Omega} + \left(\text{Other terms}\right) = 0
\end{align*}
$$

The other terms are:
- For base flow: $\left(\text{Other terms}\right)=\left(\check{u}_i\hat{n}_j,\rho u_i u_j\right)_{\partial\Omega} + \left(\check{p}\hat{n}_i, \rho u_i\right)_{\partial\Omega}$
- For perturbations: $\left(\text{Other terms}\right)=\left(\check{u}_y\hat{n}_x,\rho u_y u_x\right)_{\partial\Omega} + \left(\check{u}_x\hat{n}_x,I_{x\text{mom},i}\right)_{\Gamma_i} + \left(\check{p}\hat{n}_x, I_{\text{mass},i}\right)_{\Gamma_i} + \left(\check{u}_x\hat{n}_x,I_{x\text{mom},o}\right)_{\Gamma_o} + \left(\check{p}\hat{n}_x, I_{\text{mass},o}\right)_{\Gamma_o}$

This weak formulation has been implemented in the equations file for this example: [eqns_brokof_etal_2024.idp](./eqns_brokof_etal_2024.idp).

NOTE: This code uses computational coordinates that differ from the physical coordinates by a scaling factor related to the parameter $L$ (see the `X()` macro in [settings_brokof_etal_2024.idp](./settings_brokof_etal_2024.idp)). Note that ParaView files are exported using the physical coordinates.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/brokof_etal_2024/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/brokof_etal_2024/eqns_brokof_etal_2024.idp eqns.idp
ln -sf examples/brokof_etal_2024/settings_brokof_etal_2024.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from `.geo` files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/brokof_etal_2024 -dir $workdir -mi duct.geo
```
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/brokof_etal_2024/duct.md -mo $workdir/duct
```

## Perform parallel computations using `ff-bifbox`
### Zeroth order
1. Compute an initial base state at $Re=200$ on the created mesh from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi duct.msh -fo ignite_0 -Re 200 -Pe 70 -Ma 0.01 -gamma 1.4 -dT 5.67 -Da 1 -Ze 0 -L 1
```

2. Gradually ignite the base flow via continuation of $Da$ and $Ze$. (slow!)
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_0.base -fo ignite_1 -Ze 1 -mo ignite_1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_1.base -fo ignite_2 -Ze 2 -Da 2 -mo ignite_2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_2.base -fo ignite_3 -Ze 4 -Da 4 -mo ignite_3
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_3.base -fo ignite_4 -Ze 4.2 -Da 10 -mo ignite_4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_4.base -fo ignite_5 -Ze 4.5 -Da 20 -mo ignite_5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_5.base -fo ignite_6 -Ze 5 -Da 30 -mo ignite_6
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_6.base -fo ignite_7 -Ze 6 -Da 50 -mo ignite_7
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_7.base -fo ignite_8 -Ze 7 -Da 70 -mo ignite_8
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_8.base -fo ignite_9 -Ze 8 -Da 80 -mo ignite_9
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_9.base -fo ignite_10 -Ze 10 -Da 100 -mo ignite_10
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi ignite_10.base -fo ignite -param Da -count 10 -h0 10 -maxcount -1 -scount 5 -mo ignite -paramtarget 1700 -contorder 2
cd $workdir && export lastfile=$(printf '%s\n' ignite_*.base | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo ignited -Da 1700 -mo ignited -hmax 0.1
```

3. Compute $Re=200$, $500$, $800$ base flow fields at $L=0.5$, $1$, $5.0$.
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignited.base -fo Re200L1 -Re 200 -mo Re200L1 -hmax 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re200L1.base -fo Re500L1 -Re 500 -mo Re500L1 -hmax 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re500L1.base -fo Re800L1 -Re 800 -mo Re800L1 -hmax 0.1

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re200L1.base -fo Re200L0p5 -L 0.5 -mo Re200L0p5 -hmax 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re200L0p5.base -fo Re500L0p5 -Re 500 -mo Re500L0p5 -hmax 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re500L0p5.base -fo Re800L0p5 -Re 800 -mo Re800L0p5 -hmax 0.1

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re200L1.base -fo Re200L5 -L 5.0 -mo Re200L5 -hmax 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re200L5.base -fo Re500L5 -Re 500 -mo Re500L5 -hmax 0.1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re500L5.base -fo Re800L5 -Re 800 -mo Re800L5 -hmax 0.1
```

### First order
1. Compute the FTFs and forced response fields at $St=1$. Note that, according to the settings file, the `-sym` argument activates the acoustic characteristic BC in this setup (it does not influence the modes' symmetry).
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re800L0p5.base -fo Re800L0p5 -mo Re800L0p5 -hmax 0.05
ff-mpirun -np $nproc respcompute.md -v 0 -dir $workdir -fi Re800L0p5.base -so Re800L0p5 -Rin 0 -Rout -1 -sym 1 -omega 0 -nomega 64 -omegaf 12.6
ff-mpirun -np $nproc respcompute.md -v 0 -dir $workdir -fi Re800L0p5.base -fo Re800L0p5 -Rin 0 -Rout -1 -sym 1 -omega 6.28318530718 -pv 1

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re800L1.base -fo Re800L1 -mo Re800L1 -hmax 0.05
ff-mpirun -np $nproc respcompute.md -v 0 -dir $workdir -fi Re800L1.base -so Re800L1 -Rin 0 -Rout -1 -sym 1 -omega 0 -nomega 64 -omegaf 12.6
ff-mpirun -np $nproc respcompute.md -v 0 -dir $workdir -fi Re800L1.base -fo Re800L1 -Rin 0 -Rout -1 -sym 1 -omega 6.28318530718 -pv 1

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re800L5.base -fo Re800L5 -mo Re800L5 -hmax 0.01
ff-mpirun -np $nproc respcompute.md -v 0 -dir $workdir -fi Re800L5.base -so Re800L5 -Rin 0 -Rout -1 -sym 1 -omega 0 -nomega 127 -omegaf 12.6
ff-mpirun -np $nproc respcompute.md -v 0 -dir $workdir -fi Re800L5.base -fo Re800L5 -Rin 0 -Rout -1 -sym 1 -omega 6.28318530718 -pv 1
```

2. Compute the eigenspectra at $Re=200$, $500$, $800$ for $L=0.5$ for various $R_{\text{out}}$ values.
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re200L0p5.base -fo Re200L0p5 -mo Re200L0p5 -hmax 0.05
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re200L0p5.base -so Re200L0p5 -Rout 0 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re200L0p5.base -so Re200L0p5 -Rout -0.25 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re200L0p5.base -so Re200L0p5 -Rout -0.5 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re200L0p5.base -so Re200L0p5 -Rout -0.75 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re200L0p5.base -so Re200L0p5 -Rout -1 -sym 1 -eps_target 0.5+6i -eps_nev 50

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re500L0p5.base -fo Re500L0p5 -mo Re500L0p5 -hmax 0.05
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re500L0p5.base -so Re500L0p5 -Rout 0 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re500L0p5.base -so Re500L0p5 -Rout -0.25 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re500L0p5.base -so Re500L0p5 -Rout -0.5 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re500L0p5.base -so Re500L0p5 -Rout -0.75 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re500L0p5.base -so Re500L0p5 -Rout -1 -sym 1 -eps_target 0.5+6i -eps_nev 50

ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re800L0p5.base -so Re800L0p5 -Rout 0 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re800L0p5.base -so Re800L0p5 -Rout -0.25 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re800L0p5.base -so Re800L0p5 -Rout -0.5 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re800L0p5.base -so Re800L0p5 -Rout -0.75 -sym 1 -eps_target 0.5+6i -eps_nev 50
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re800L0p5.base -so Re800L0p5 -Rout -1 -sym 1 -eps_target 0.5+6i -eps_nev 50
```

3. Perform resolvent analysis at $Re=200$, $500$, $800$ for $L=0.5$ for tuned $R_{\text{out}}$ values where $\sigma\sim-0.45$.
```sh
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi Re200L0p5.base -so Re200L0p5 -Rout -0.24125 -sym 1 -omega 0.1 -omegaf 12.6 -nomega 127
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi Re500L0p5.base -so Re500L0p5 -Rout -0.70125 -sym 1 -omega 0.1 -omegaf 12.6 -nomega 127
ff-mpirun -np $nproc rslvcompute.md -v 0 -dir $workdir -fi Re800L0p5.base -so Re800L0p5 -Rout -0.92 -sym 1 -omega 0.1 -omegaf 12.6 -nomega 127
```