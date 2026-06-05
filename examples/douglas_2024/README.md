# 2D Incompressible Swirling Flow Example: Douglas, TCFD, (2024)
This file shows an example `ff-bifbox` workflow for reproducing the results in section III.B of the study:
```bibtex
@article{douglas_2024,
  title={A Balanced Outflow Boundary Condition for Swirling Flows},
  journal={Theoretical and Computational Fluid Dynamics},
  publisher={Springer Nature},
  author={Douglas, Christopher M.},
  year={2024},
  DOI={10.1007/s00162-024-00701-5},
}
```
The commands below illustrate how to perform the analysis of the Grabowski--Berger vortex as in the paper using `ff-bifbox`. Note: a reproducer for the portion of the analysis focused on the rotating pipe flow is given in the Supplementary Materials (included here in [example1_suppmat.md](./example1_suppmat.md)). Note that the reference solution is computed using settings and equations from the [meliga_etal_2012](../meliga_etal_2012/) example.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{1}{Re}\frac{\partial^2u_i}{\partial x_j^2} &= 0 \\
\frac{\partial u_i}{\partial x_i} &= 0
\end{align*}
$$

Depending on the chosen BC set, these are complemented by additional governing equations for the boundary variables:
- $`\frac{\partial p_o}{\partial x_i}\hat{t}_i - \frac{u_{\theta}^2}{r} = 0`$, with modified outflow BC
- $`\frac{\partial^2 \phi_o}{\partial x_i^2} - \frac{\partial}{\partial x_i}\left(-\frac{u_{\theta}^2}{r}\hat{e}_r+\frac{u_ru_\theta}{r}\hat{e}_{\theta}\right) = 0`$, with balanced outflow BC

together with the boundary conditions:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $u_x=1$, $u_r=0$, $u_{\theta}=S\Psi(r)$, $p_o=\phi_o=0$ |
| Axis, $\Gamma_a$| $`\begin{cases}\frac{\partial u_x}{\partial r}=u_r=u_{\theta}=\frac{\partial \phi_o}{\partial r}=0, & \text{if } m=0 \\\\ u_x=\frac{\partial u_r}{\partial r}=\frac{\partial u_{\theta}}{\partial r}=\phi_o=0, & \text{if } \|m\|=1 \\\\ u_x=u_r=u_{\theta}=\phi_o=0, & \text{if } \|m\|>1\end{cases}`$ |
| Open, $\Gamma_o$ | $`\begin{cases}\frac{1}{Re}\frac{\partial u_i}{\partial x}-p\hat{e}_x = 0, & \text{if free outflow BC} \\\\\frac{1}{Re}\frac{\partial u_i}{\partial x}-\left(p-p_o\right)\hat{e}_x = 0, &\text{if modified outflow BC} \\\\ \frac{1}{Re}\frac{\partial u_i}{\partial x}-\left(p-\phi_o\right)\hat{e}_x = 0, & \text{if balanced outflow BC} \\\\ \frac{\partial u_i}{\partial t}+c_j\frac{\partial u_i}{\partial x_j} = 0, & \text{if convective BC}\end{cases}`$ |

where $`\Psi(r)=\begin{cases} r(2-r^2), & \text{if }r \leq 1 \\\\ \frac{1}{r}, & \text{if } r > 1 \end{cases}`$.

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_a+\Gamma_o$. Solutions $\vec{q}$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}$,

$$
\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} + \left(\text{Other terms}\right) = 0.
$$

Here:
- $`\vec{q}=\begin{cases}\left(u_i,p\right)^T, & \text{if free outflow BC} \\\\ \left(u_i,p,p_o\right)^T, & \text{if modified outflow BC} \\\\ \left(u_i,p,\phi_o\right)^T, & \text{if balanced outflow BC} \\\\ \left(u_i,p,c\right)^T, & \text{if convective BC}\end{cases}`$

The "other terms" are:
- $`\left(\text{Other terms}\right)=\begin{cases}\left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega}, & \text{if free outflow} \\\\ \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega}+\left(\check{u}_i,p_o\hat{n}_i\right)_{\Gamma_o} + \left(\check{p}_o,\frac{\partial p_o}{\partial x_i}\hat{t}_i-\frac{u_{\theta}^2}{r}\right)_{\Gamma_o}, & \text{if modified outflow} \\\\ \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\check{u}_i,\phi_o\hat{n}_i\right)_{\Gamma_o} + \left(\frac{\partial \check{\phi}_o}{\partial x_i},\frac{\partial \phi_o}{\partial x_i} + \frac{u_{\theta}^2}{r}\hat{e}_r - \frac{u_ru_\theta}{r}\hat{e}_{\theta}\right)_{\Gamma_o}, & \text{if balanced outflow} \\\\ \left(\check{u}_i,\frac{\partial p}{\partial x_i}\right)_{\Omega} + \left(\check{p},c\right)_{\Omega} + \left(\check{c},\bar{p}\right)_{\Omega} + \left(\check{u}_i,\frac{1}{C_i\hat{n}_i Re}\frac{\partial \hat{u}_i}{\partial t}\right)_{\Gamma_o}, & \text{if convective} \end{cases}`$

These weak formulations have been implemented in the equations files for this example: [eqns_douglas_2024_freeout.idp](./eqns_douglas_2024_freeout.idp), [eqns_douglas_2024_modified.idp](./eqns_douglas_2024_modified.idp), [eqns_douglas_2024_balanced.idp](./eqns_douglas_2024_balanced.idp), [eqns_douglas_2024_convect.idp](./eqns_douglas_2024_convect.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/douglas_2024/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/meliga_etal_2012/eqns_meliga_etal_2012.idp eqns.idp
ln -sf examples/meliga_etal_2012/settings_meliga_etal_2012.idp settings.idp
```
OR
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_convect.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_convect.idp settings.idp
```
OR
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_balanced.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_balanced.idp settings.idp
```
OR
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_modified.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_modified.idp settings.idp
```
OR
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_freeout.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_freeout.idp settings.idp
```

## Build initial meshes
Note: this example does not make use of adaptive meshing.
```sh
FreeFem++-mpi -v 0 examples/douglas_2024/grabowski.md -mo $workdir/G
```
## Run the code provided in the Supplementary Materials
```sh
FreeFem++-mpi -v 0 examples/douglas_2024/example1_suppmat.md
```

## Perform parallel computations for Grabowski--Berger vortex flow using `ff-bifbox`
### Steady axisymmetric dynamics
1. Compute reference base states on the largest mesh at $Re=200$, $S=0.85$ to $S=1.3$ from default guess.
```sh
ln -sf examples/meliga_etal_2012/eqns_meliga_etal_2012.idp eqns.idp
ln -sf examples/meliga_etal_2012/settings_meliga_etal_2012.idp settings.idp
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi G3.msh -fo meligaS0p85 -1/Re 0.005 -S 0.85 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi meligaS0p85.base -fo meligaS0p9 -S 0.9 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi meligaS0p9.base -fo meligaS1 -S 1.0 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi meligaS1.base -fo meligaS1p095 -S 1.095 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi meligaS1p095.base -fo meligaS1p3 -S 1.3 -snes_rtol 0
```

2. Compute base states with convective BC on the truncated mesh.
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_convect.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_convect.idp settings.idp
ff-mpirun -np $nproc examples/douglas_2024/basecomputeaug.md -v 0 -dir $workdir -mi G2.msh -fi meligaS0p85.base -fo convectS0p85 -snes_rtol 0 -pv 1
FreeFem++-mpi -v 0 examples/douglas_2024/computebaseerror.md -fi convect -ci S0p85 -fo convectS0p85err -dir $workdir -pv 1
ff-mpirun -np $nproc examples/douglas_2024/basecomputeaug.md -v 0 -dir $workdir -mi G2.msh -fi meligaS0p9.base -fo convectS0p9 -snes_rtol 0
ff-mpirun -np $nproc examples/douglas_2024/basecomputeaug.md -v 0 -dir $workdir -mi G2.msh -fi meligaS1.base -fo convectS1 -snes_rtol 0 -pv 1
FreeFem++-mpi -v 0 examples/douglas_2024/computebaseerror.md -fi convect -ci S1 -fo convectS1err -dir $workdir -pv 1
ff-mpirun -np $nproc examples/douglas_2024/basecomputeaug.md -v 0 -dir $workdir -mi G2.msh -fi meligaS1p095.base -fo convectS1p095 -snes_rtol 0
ff-mpirun -np $nproc examples/douglas_2024/basecomputeaug.md -v 0 -dir $workdir -mi G2.msh -fi meligaS1p3.base -fo convectS1p3 -snes_rtol 0
```

3. Compute base states with balanced outflow BC on the truncated mesh.
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_balanced.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_balanced.idp settings.idp
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi G2.msh -fo balancedS0p85 -1/Re 0.005 -S 0.85 -snes_rtol 0 -pv 1
FreeFem++-mpi -v 0 examples/douglas_2024/computebaseerror.md -fi balanced -ci S0p85 -fo balancedS0p85err -dir $workdir -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi balancedS0p85.base -fo balancedS0p9 -S 0.9 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi balancedS0p9.base -fo balancedS1 -S 1 -snes_rtol 0 -pv 1
FreeFem++-mpi -v 0 examples/douglas_2024/computebaseerror.md -fi balanced -ci S1 -fo balancedS1err -dir $workdir -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi balancedS1.base -fo balancedS1p095 -S 1.095 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi balancedS1p095.base -fo balancedS1p3 -S 1.3 -snes_rtol 0
```

4. Compute base states with free outflow BC on the truncated mesh.
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_freeout.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_freeout.idp settings.idp
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi G2.msh -fo freeoutS0p5 -1/Re 0.005 -S 0.5 -snes_rtol 0
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi freeoutS0p5.base -fo freeout -param S -h0 10 -paramtarget 0.85 -maxcount -1 -scount 25
cd $workdir && export lastfile=$(printf '%s\n' freeout_*.base | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo freeoutS0p85 -S 0.85 -snes_rtol 0 -pv 1
FreeFem++-mpi -v 0 examples/douglas_2024/computebaseerror.md -fi freeout -ci S0p85 -fo freeoutS0p85err -dir $workdir -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi freeoutS0p85.base -fo freeoutS0p9 -S 0.9 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi freeoutS0p9.base -fo freeoutS1 -S 1 -snes_rtol 0 -pv 1
FreeFem++-mpi -v 0 examples/douglas_2024/computebaseerror.md -fi freeout -ci S1 -fo freeoutS1err -dir $workdir -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi freeoutS1.base -fo freeoutS1p095 -S 1.095 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi freeoutS1p095.base -fo freeoutS1p2 -S 1.2 -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi freeoutS1p2.base -fo freeoutS1p3 -S 1.3 -snes_rtol 0
```

### Linear 3-D dynamics
1. Compute reference eigenvalues/eigenmodes on the largest mesh at $Re=200$, $S=1$ and $S=1.3$.
```sh
ln -sf examples/meliga_etal_2012/eqns_meliga_etal_2012.idp eqns.idp
ln -sf examples/meliga_etal_2012/settings_meliga_etal_2012.idp settings.idp
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi meligaS1.base -fo meligaS1m1 -sym -1 -eps_pos_gen_non_hermitian -eps_target 0.1+1.2i -eps_nev 1 -eps_two_sided 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi meligaS1p3.base -fo meligaS1p3m2 -sym -2 -eps_pos_gen_non_hermitian -eps_target 0.1+2.5i -eps_nev 1 -eps_two_sided 1
```

2. Compute eigenmodes with convective BC on the truncated mesh.
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_convect.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_convect.idp settings.idp
ff-mpirun -np $nproc examples/douglas_2024/modecomputeaug.md -v 0 -dir $workdir -mi G2.msh -fi meligaS1.base -fo convectS1m1 -Cx 0.6 -Cr 0.1 -sym -1 -eps_pos_gen_non_hermitian -eps_target 0.1+1.2i -eps_nev 1
FreeFem++-mpi -v 0 examples/douglas_2024/computemodeerror.md -Cx 0.6 -Cr 0.1 -fi convect -ci S1m1 -fo convectS1m1err -dir $workdir -pv 1
ff-mpirun -np $nproc examples/douglas_2024/modecomputeaug.md -v 0 -dir $workdir -mi G2.msh -fi meligaS1p3.base -fo convectS1p3m2 -Cx 0.6 -Cr 0.1 -sym -2 -eps_pos_gen_non_hermitian -eps_target 0.1+2.5i -eps_nev 1
FreeFem++-mpi -v 0 examples/douglas_2024/computemodeerror.md -Cx 0.6 -Cr 0.1 -fi convect -ci S1p3m2 -fo convectS1p3m2err -dir $workdir -pv 1
```

3. Compute eigemodes with balanced BC on the truncated mesh.
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_balanced.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_balanced.idp settings.idp
FreeFem++-mpi -v 0 examples/douglas_2024/basefieldappend.md -fi meligaS1.base -fo meligaappendS1 -dir $workdir
FreeFem++-mpi -v 0 examples/douglas_2024/modefieldappend.md -fi meligaS1m1.mode -fo meligaappendS1m1 -dir $workdir
FreeFem++-mpi -v 0 examples/douglas_2024/modefieldappend.md -fi meligaS1m1adj.mode -fo meligaappendS1m1adj -dir $workdir
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi G2.msh -fi meligaappendS1.base -fo balancedS1m1 -sym -1 -eps_target 0.1+1.2i -eps_nev 1
FreeFem++-mpi -v 0 examples/douglas_2024/computemodeerror.md -ri meligaappend -fi balanced -ci S1m1 -fo balancedS1m1err -dir $workdir -pv 1
FreeFem++-mpi -v 0 examples/douglas_2024/basefieldappend.md -fi meligaS1p3.base -fo meligaappendS1p3 -dir $workdir
FreeFem++-mpi -v 0 examples/douglas_2024/modefieldappend.md -fi meligaS1p3m2.mode -fo meligaappendS1p3m2 -dir $workdir
FreeFem++-mpi -v 0 examples/douglas_2024/modefieldappend.md -fi meligaS1p3m2adj.mode -fo meligaappendS1p3m2adj -dir $workdir
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi G2.msh -fi meligaappendS1p3.base -fo balancedS1p3m2 -sym -2 -eps_target 0.1+2.5i -eps_nev 1
FreeFem++-mpi -v 0 examples/douglas_2024/computemodeerror.md -ri meligaappend -fi balanced -ci S1p3m2 -fo balancedS1p3m2err -dir $workdir -pv 1
```

4. Compute eigenmodes with free outflow BC on the truncated mesh.
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_freeout.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_freeout.idp settings.idp
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi G2.msh -fi meligaS1.base -fo freeoutS1m1 -sym -1 -eps_pos_gen_non_hermitian -eps_target 0.1+1.2i -eps_nev 1
FreeFem++-mpi -v 0 examples/douglas_2024/computemodeerror.md -fi freeout -ci S1m1 -fo freeoutS1m1err -dir $workdir -pv 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi G2.msh -fi meligaS1p3.base -fo freeoutS1p3m2 -sym -2 -eps_pos_gen_non_hermitian -eps_target 0.1+2.5i -eps_nev 1
FreeFem++-mpi -v 0 examples/douglas_2024/computemodeerror.md -fi freeout -ci S1p3m2 -fo freeoutS1p3m2err -dir $workdir -pv 1
```

5. Compute eigemodes with modified BC on the truncated mesh.
```sh
ln -sf examples/douglas_2024/eqns_douglas_2024_modified.idp eqns.idp
ln -sf examples/douglas_2024/settings_douglas_2024_modified.idp settings.idp
FreeFem++-mpi -v 0 examples/douglas_2024/basefieldappend.md -fi meligaS1.base -fo meligaappend2S1 -dir $workdir
FreeFem++-mpi -v 0 examples/douglas_2024/modefieldappend.md -fi meligaS1m1.mode -fo meligaappend2S1m1 -dir $workdir
FreeFem++-mpi -v 0 examples/douglas_2024/modefieldappend.md -fi meligaS1m1adj.mode -fo meligaappend2S1m1adj -dir $workdir
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi G2.msh -fi meligaappend2S1.base -fo modifiedS1m1 -sym -1 -eps_target 0.1+1.2i -eps_nev 1
FreeFem++-mpi -v 0 examples/douglas_2024/computemodeerror.md -ri meligaappend2 -fi modified -ci S1m1 -fo modifiedS1m1err -dir $workdir -pv 1
FreeFem++-mpi -v 0 examples/douglas_2024/basefieldappend.md -fi meligaS1p3.base -fo meligaappend2S1p3 -dir $workdir
FreeFem++-mpi -v 0 examples/douglas_2024/modefieldappend.md -fi meligaS1p3m2.mode -fo meligaappend2S1p3m2 -dir $workdir
FreeFem++-mpi -v 0 examples/douglas_2024/modefieldappend.md -fi meligaS1p3m2adj.mode -fo meligaappend2S1p3m2adj -dir $workdir
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -mi G2.msh -fi meligaappend2S1p3.base -fo modifiedS1p3m2 -sym -2 -eps_target 0.1+2.5i -eps_nev 1
FreeFem++-mpi -v 0 examples/douglas_2024/computemodeerror.md -ri meligaappend2 -fi modified -ci S1p3m2 -fo modifiedS1p3m2err -dir $workdir -pv 1
```