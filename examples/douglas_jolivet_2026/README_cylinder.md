# Examples: Douglas & Jolivet (2026)
This file shows an example `ff-bifbox` workflow for reproducing the results of the study:
```bibtex
@article{douglas_jolivet_2026,
author = {Douglas, Christopher M. and Jolivet, Pierre},
title = {{f}f-bifbox: {A} scalable, open-source toolbox for bifurcation analysis of nonlinear {PDE}s},
journal = {Computer Physics Communications},
volume = {326},
pages = {110221},
year = {2026},
Publisher = {Elsevier},
issn = {0010-4655},
doi = {10.1016/j.cpc.2026.110221},
url = {https://doi.org/10.1016/j.cpc.2026.110221},
arxiv = {https://arxiv.org/abs/2509.18429},
}
```
The commands below illustrate how to analyze a 2-D compressible flow past a cylinder using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\left(1+\gamma M^2 p\right)\left[\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right] + T\frac{\partial p}{\partial x_i} - \frac{T}{Re}\frac{\partial\tau_{ij}}{\partial x_j} + \beta_sd_s^2\left(u_i-\hat{e}_x\right) &= 0 \\
\left(1+\gamma M^2 p\right)\left[\frac{\partial T}{\partial t} + u_i\frac{\partial T}{\partial x_i} + \left(\gamma-1\right)T\frac{\partial u_i}{\partial x_i}\right] - \frac{\gamma\left(\gamma-1\right)M^2}{Re}T\tau_{ij}\frac{\partial u_i}{\partial x_j} \\- \frac{\gamma T}{Pr Re}\frac{\partial^2T}{\partial x_i^2} + \beta_sd_s^2\gamma\left(T-1\right) &= 0 \\
\gamma M^2 T\left(\frac{\partial p}{\partial t} + u_i\frac{\partial p}{\partial x_i}\right) - \left(1+\gamma M^2 p\right)\left[\frac{\partial T}{\partial t} + u_i\frac{\partial T}{\partial x_i} - T\frac{\partial u_i}{\partial x_i}\right]+ \beta_s d_s^2\gamma M^2 p&= 0
\end{align*}
$$

where $\tau_{ij}=\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}-\delta_{ij}\frac{\partial u_k}{\partial x_k}$.

The boundary conditions are:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $u_x=T=1$, $u_y=0$ |
| Wall, $\Gamma_w$| $u_x=u_y=\frac{\partial T}{\partial x_i}\hat{n}_i=0$ |
| Axis, $\Gamma_a$| $\frac{\partial u_x}{\partial y}=u_y=\frac{\partial T}{\partial y}=0$|
| Lateral, $\Gamma_l$ | $\frac{\partial u_x}{\partial y}=u_y=\frac{\partial T}{\partial y}=0$ |
| Outlet, $\Gamma_o$ | $\frac{1}{Re}\tau_{ix}-p\hat{e}_x = \frac{\partial T}{\partial x}=0$ |

The present implementation is based on a weak formulation of these equations. The equations are integrated over the planar domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_w+\Gamma_a+\Gamma_l+\Gamma_o$. Solutions $\vec{q}=\left(u_i,T,p\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{T},\check{p}\right)^T$,

$$
\begin{align*} 
&\left(\check{u}_i,\left(1 + \gamma M^2p\right)\left[\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right] + \beta_sd_s^2\left(u_i-\hat{e}_x\right)\right)_{\Omega} - \left(\frac{\partial \left(\check{u}_i T\right)}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \left(\check{u}_i T\right)}{\partial x_j},\tau_{ij}\right)_{\Omega} \\
&+\left(\check{T},\left(1 + \gamma M^2p\right)\left[\frac{\partial T}{\partial t} + u_i\frac{\partial T}{\partial x_i} + \left(\gamma-1\right)T \frac{\partial u_i}{\partial x_i}\right] - \frac{\gamma\left(\gamma-1\right)M^2}{Re}T\tau_{ij}\frac{\partial u_i}{\partial x_j}+\beta_sd_s^2\gamma\left(T-1\right)\right)_{\Omega} + \left(\frac{\partial\left(\check{T} T\right)}{\partial x_i},\frac{\gamma}{Pr Re}\frac{\partial T}{\partial x_i}\right)_{\Omega} \\
&+ \left(\check{p},\gamma M^2T\left[\frac{\partial p}{\partial t} + u_i\frac{\partial p}{\partial x_i}\right] - \left(1 + \gamma M^2p\right)\left[\frac{\partial T}{\partial t} + u_i\frac{\partial T}{\partial x_i} - T\frac{\partial u_i}{\partial x_i}\right] + \beta_sd_s^2\gamma M^2 p\right)_{\Omega} = 0.
\end{align*}
$$

This weak formulation has been implemented in the equations file for this example: [eqns_douglas_jolivet_2026_cylinder.idp](./eqns_douglas_jolivet_2026_cylinder.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/douglas_jolivet_2026/cylinder
export nproc=4
```

# Example 3: compressible flow past a cylinder

1. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/douglas_jolivet_2026/eqns_douglas_jolivet_2026_cylinder.idp eqns.idp
ln -sf examples/douglas_jolivet_2026/settings_douglas_jolivet_2026_cylinder.idp settings.idp
```

2. Build initial mesh
```sh
FreeFem++-mpi -v 0 examples/douglas_jolivet_2026/mesh_cylinder.md -mo $workdir/cylinder
```

3. Compute base states on the created meshes at $Re=40,50,70,90$ and $Ma=0,0.4,0.6,0.8$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi cylinder.msh -fo cylinder_Ma0p0_Re40 -1/Re 0.025 -1/Pr 1.38888888889 -Ma^2 0 -gamma 1.4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p0_Re40.base -fo cylinder_Ma0p0_Re50 -1/Re 0.02
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p0_Re50.base -fo cylinder_Ma0p0_Re70 -1/Re 0.014285714285714285
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p0_Re70.base -fo cylinder_Ma0p0_Re90 -1/Re 0.0111111111111111111
for Re in 40 50 70 90; do
  ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p0_Re"$Re".base -fo cylinder_Ma0p4_Re"$Re" -Ma^2 0.16
  ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p4_Re"$Re".base -fo cylinder_Ma0p6_Re"$Re" -Ma^2 0.36
  ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p6_Re"$Re".base -fo cylinder_Ma0p8_Re"$Re" -Ma^2 0.64
  for Ma in 0 4 6 8; do
    ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p"$Ma"_Re"$Re".base -fo cylinder_Ma0p"$Ma"_Re"$Re" -mo cylinder_Ma0p"$Ma"_Re"$Re" -thetamax 0.1 -hmin 1.e-5 -hmax 2
  done
done
```

4. Continue base states along $Re$ or $Ma$ with adaptive remeshing
```sh
for Ma in 0 4 6 8; do
  ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p"$Ma"_Re40.base -fo cylinder_Ma0p"$Ma"_Re20 -1/Re 0.05
  ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi cylinder_Ma0p"$Ma"_Re20.base -fo cylinder_Ma0p"$Ma" -param 1/Re -h0 -1 -scount 2 -maxcount -1 -paramtarget 0.01 -mo cylinder_Ma0p"$Ma" -thetamax 0.1 -hmin 1.e-5 -hmax 2
done
for Re in 40 50 70 90; do
  ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi cylinder_Ma0p0_Re"$Re".base -fo cylinder_Re"$Re" -param Ma^2 -h0 0.25 -scount 2 -maxcount -1 -paramtarget 0.81 -mo cylinder_Re"$Re" -thetamax 0.1 -hmin 1.e-5 -hmax 2
done
```

5. Compute Hopf bifurcation at $Ma=0$ and trace it along the neutral curve
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p0_Re50.base -fo cylinder_Ma0p0_Re47p6 -1/Re 0.021
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi cylinder_Ma0p0_Re47p6.base -fo cylinder_Ma0p0_Re47p6 -eps_target 0.1+0.7i -sym 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder_Ma0p0_Re47p6.mode -fo cylinder_Ma0p0 -param 1/Re -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder_Ma0p0.hopf -fo cylinder_Ma0p0 -mo cylinder_Ma0p0 -param 1/Re -thetamax 0.1 -hmin 1e-5 -hmax 2
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi cylinder_Ma0p0.hopf -fo cylinder -mo cylinder -thetamax 0.1 -hmin 1e-5 -hmax 2 -param 1/Re -param2 Ma^2 -h0 0.1 -scount 3 -paramtarget 0.01 -maxcount -1
```

6. Compute Hopf bifurcations at $Ma=0.4,0.6,0.8$ and $Re=50,70,90$
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder_12.hopf -fo cylinder_Ma0p4 -param 1/Re -Ma^2 0.16
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder_15specialpt.hopf -fo cylinder_Ma0p6 -param 1/Re -Ma^2 0.36
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder_21.hopf -fo cylinder_Ma0p8 -param 1/Re -Ma^2 0.64
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder_15specialpt.hopf -fo cylinder_Re50 -param Ma^2 -1/Re 0.02
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder_18.hopf -fo cylinder_Re70 -param Ma^2 -1/Re 0.014285714285714285
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi cylinder_21.hopf -fo cylinder_Re90 -param Ma^2 -1/Re 0.0111111111111111111
```

7. Continue the branches of periodic solutions emanating from the Hopf points along $Re$ and $Ma$.
```sh
for Ma in 0 4 6 8; do
  ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi cylinder_Ma0p"$Ma".hopf -fo cylinder_Ma0p"$Ma" -mo cyl_Ma0p"$Ma" -thetamax 0.1 -hmin 1e-5 -hmax 2 -param 1/Re -h0 1 -scount 4 -maxcount -1 -paramtarget 0.01 -Nh 3 -fieldsplit_0_fieldsplit_0_mat_mumps_icntl_35 1 -fieldsplit_0_fieldsplit_0_mat_mumps_cntl_7 1.0e-8
done
for Re in 50 70 90; do
  ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi cylinder_Re"$Re".hopf -fo cylinder_Re"$Re" -mo cyl_Re"$Re" -thetamax 0.1 -hmin 1e-5 -hmax 2 -param Ma^2 -h0 1 -scount 4 -maxcount -1 -paramtarget 0.01 -Nh 3 -fieldsplit_0_fieldsplit_0_mat_mumps_icntl_35 1 -fieldsplit_0_fieldsplit_0_mat_mumps_cntl_7 1.0e-8
done
```