# Low-Mach Inverted Conical Flame Example: Schulke, Wang & Douglas, (2026)
This file shows an example `ff-bifbox` workflow for reproducing the results in the study:
```bibtex
@article{schulke_etal_2026,
	title = {A harmonic-balance model of subcritical instability in lean premixed inverted conical flames},
	author = {Schulke, Gretchen and Wang, Chuhan and Douglas, Christopher M.},
	year = {2026},
	journal = {Journal of Engineering for Gas Turbines and Power},
  publisher = {The American Society of Mechanical Engineers},
	note = {accepted manuscript}
}
```
The commands below illustrate how to perform a bifurcation analysis of a lean premixed inverted conical flame in an axisymmetric annular jet using `ff-bifbox`. (See also the previous study [wang_etal_2024](../wang_etal_2024).) The resulting system is very stiff and exhibits strong non-normality, which both represent significant numerical challenges. As such, the workflow given here is relatively fragile and may need to be adapted using various

The dimensionless governing equations are given as 

$$
\begin{aligned}
\frac{\partial u_i}{\partial t} + u_i \frac{\partial u_j}{\partial x_j} + T \frac{\partial p}{\partial x_i} - \frac{T}{Re} \frac{\partial\tau_{ij}}{\partial x_j} &= 0, \\
\frac{\partial Y}{\partial t} + u_i \frac{\partial Y}{\partial x_i} + T \mathcal{Q} - \frac{T}{Re Pr Le} \frac{\partial}{\partial x_i} \left( \mu \frac{\partial Y}{\partial x_i} \right) &= 0, \\
\frac{\partial T}{\partial t} + u_i \frac{\partial T}{\partial x_i} - \frac{\Delta h_c \min\left(1, \phi\right)}{\mathrm{AFR}_{m,\mathrm{st}} + \phi} T \mathcal{Q} - \frac{T}{Re Pr} \frac{\partial}{\partial x_i} \left( \mu \frac{\partial T}{\partial x_i} \right) &= 0, \\
\frac{\partial T}{\partial t} + u_i\frac{\partial T}{\partial x_i} - T \frac{\partial u_i}{\partial x_i} &= 0,
\end{aligned}
$$

where:
- $`\tau_{ij} = \mu\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}-\frac{2}{3}\delta_{ij}\frac{\partial u_k}{\partial x_k}\right)`$
- $`\mu = T^{3/2}\left(\frac{1+S}{T+S}\right)`$
- $`\mathcal{Q} = Da\sqrt{\frac{\min\left(1,\phi\right)}{\mathrm{AFR}_{m,\mathrm{st}}+\phi}}\left[\frac{Y+\max\left(0,\phi-1\right)}{T}\right]\sqrt{\frac{Y+\max\left(0,\phi^{-1}-1\right)}{T}}\exp\left(-\frac{T_a}{T}\right)`$
- $`Da = \frac{\mathcal{A}}{Re}\left(\frac{\mathrm{AFR}_{m,\mathrm{st}}+\phi}{\mathrm{AFR}_{v,\mathrm{st}}+\phi}\right)^{3/2}.`$

Note that the ideal gas equation in the low Mach number limit ($`1 = \rho T`$) has been used to eliminate density from the above equations.

The boundary conditions are 

| Label | Velocity | Temperature | Species |
|---|---|---|---|
| $`\Gamma_{\mathrm{in}}`$ | $`u_x = u_{\mathrm{Pl}}(r)`$, $`u_r=u_\theta=0`$ | $`T = 1`$ | $`Y = 1`$ |
| $`\Gamma_{\mathrm{nozz}}`$ | $`u_x = u_r = u_\theta = 0`$ | $`\hat{n}_i\frac{\partial T}{\partial x_i} = 0`$ | $`\hat{n}_i\frac{\partial Y}{\partial x_i} = 0`$ |
| $`\Gamma_{\mathrm{wall}}`$ | $`u_x = u_r = u_\theta = 0`$ | $`T = 1`$ | $`\hat{n}_i\frac{\partial Y}{\partial x_i} = 0`$|
| $`\Gamma_{\mathrm{cb}}`$ | $`u_x = u_r = u_\theta = 0`$ | $`T = \frac{7}{3}`$ | $`\hat{n}_i\frac{\partial Y}{\partial x_i} = 0`$ |
| $`\Gamma_{\mathrm{axis}}`$ | $`\begin{cases}\frac{\partial u_x}{\partial r}=u_r=u_{\theta}=0, & \mathrm{if\ } m=0 \\\\ u_x=\frac{\partial u_r}{\partial r}=\frac{\partial u_{\theta}}{\partial r}=0, & \mathrm{if\ } \|m\|=1 \\\\ u_x=u_r=u_{\theta}=0, & \mathrm{if\ } \|m\|>1\end{cases}`$ | $`\begin{cases}\frac{\partial T}{\partial r}=0, & \mathrm{if\ } m=0 \\\\ T=0, & \mathrm{if\ } \|m\|\geq1\end{cases}`$ | $`\begin{cases}\frac{\partial Y}{\partial r}=0, & \mathrm{if\ } m=0 \\\\ Y=0, & \mathrm{if\ } \|m\|\geq1\end{cases}`$ |
| $`\Gamma_{\mathrm{lat}}`$ | $`\hat{n}_j\frac{\mu}{Re}\frac{\partial u_i}{\partial x_j}=p\hat{n}_i+\frac{u_i}{2T}\min\left(0,\hat{n}_ju_j\right)`$ | $`T = 1`$ | $`Y = 0`$ |
| $`\Gamma_{\mathrm{out}}`$ | $`\hat{n}_j\frac{\mu}{Re}\frac{\partial u_i}{\partial x_j}=p\hat{n}_i+\frac{u_i}{2T}\min\left(0,\hat{n}_ju_j\right)`$ | $`\hat{n}_i\frac{\partial T}{\partial x_i} = 0`$ | $`\hat{n}_i\frac{\partial Y}{\partial x_i} = 0`$ |

where $`u_{\mathrm{Pl}}(r)=2\left(\frac{\left(1-4r^2\right)\log D_{\mathrm{cb}}-\left(1-D^2_{\mathrm{cb}}\right)\log(2r)}{1-D_{\mathrm{cb}}^2+\left(1+D_{cb}^2\right)\log D_{\mathrm{cb}}}\right)`$.

The present implementaton is based on a weak formulation of these equations. Test functions are introduced and the equations are integrated over the domain $\Omega$ with boundary $`\partial\Omega = \Gamma_{\mathrm{in}} + \Gamma_{\mathrm{nozz}} + \Gamma_{\mathrm{wall}} + \Gamma_{\mathrm{cb}} + \Gamma_{\mathrm{axis}} + \Gamma_{\mathrm{lat}} + \Gamma_{\mathrm{out}}`$. Solutions $`\vec{q} = \left(u_i, Y, T, p\right)^T`$ are sought, in the appropriate spaces, such that for all test functions $`\vec{\check{q}} = \left(\check{u}_i, \check{Y}, \check{T}, \check{p}\right)^T`$,

$$
\begin{aligned}
&\left(\check{u}_i,\frac{\partial u_i}{\partial t}+u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial \left(\check{u}_i T\right)}{\partial x_j},\delta_{ij}p - \frac{1}{Re}\tau_{ij}\right)_{\Omega} - \left(\check{u}_i,\hat{n}_jT\frac{\mu}{Re}\left(\frac{\partial u_j}{\partial x_i} - \tfrac{2}{3}\delta_{ij}\frac{\partial u_k}{\partial x_k}\right)+\tfrac{1}{2}u_i\min\left(0,u_j\hat{n}_j\right)\right)_{\Gamma_{\mathrm{lat}}\cup\Gamma_{\mathrm{out}}} \\
&+ \left(\check{Y},\frac{\partial Y}{\partial t}+u_i\frac{\partial Y}{\partial x_i} + T\mathcal{Q}\right)_{\Omega} + \left(\frac{\partial \left(\check{Y}T\right)}{\partial x_i},\frac{\mu}{Pr Re Le}\frac{\partial Y}{\partial x_i}\right)_{\Omega} \\
&+ \left(\check{T},\frac{\partial T}{\partial t}+u_i\frac{\partial T}{\partial x_i} - \frac{h_c\min(1,\phi)}{\mathrm{AFR}_{m,\mathrm{st}}+\phi} T\mathcal{Q}\right)_{\Omega} + \left(\frac{\partial \left(\check{T}T\right)}{\partial x_i},\frac{\mu}{Re Pr}\frac{\partial T}{\partial x_i}\right)_{\Omega} \\
&+ \left(\check{p},\frac{\partial T}{\partial t}+u_i\frac{\partial T}{\partial x_i} - T\frac{\partial u_i}{\partial x_i}\right)_{\Omega} = 0
\end{aligned}
$$

More details on the formulation can be found in the original paper.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```

2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/schulke_etal_2026/data
export nproc=4
```

3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/schulke_etal_2026/eqns_schulke_etal_2026_axi.idp eqns.idp
ln -sf examples/schulke_etal_2026/settings_schulke_etal_2026_axi.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.

#### CASE 1: Gmsh is installed - build initial mesh directly from .geo files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/schulke_etal_2026 -dir $workdir -mi Vflame.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).

#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/schulke_etal_2026/Vflame.md -mo $workdir/Vflame
```

## Perform parallel computations using `ff-bifbox`
### Laminar base flow
1. Compute a non-reacting base state with reference parameters on the initial mesh.
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi Vflame.msh -fo axi_nonreacting -Re 70 -Tr 2.3333333333333333 -Ar 0 -phi 0.8 -Dhc 128.23469709865222 -snes_linesearch_type secant -snes_rtol 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi axi_nonreacting.base -fo axi_nonreacting -Re 350 -mo nonreacting -snes_linesearch_type secant -snes_rtol 0 -localrefinement 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi axi_nonreacting.base -fo axi_nonreacting -Re 1500 -mo nonreacting -snes_linesearch_type secant -snes_rtol 0 -localrefinement 0
```

2. Turn on chemistry and ignite the $Re = 1500$ flow at an elevated centerbody temperature and lower combustion enthalpy. Then perform continuation back to reference parameters. Coarse meshes are used for computational efficiency and stabilizing artificial dissipation.

```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi axi_nonreacting.base -fo axi_ignite_0 -Tr 3.3333333333333333 -Ar 423521926.87072223 -Dhc 12.823469709865222 -mo ignite_0 -snes_rtol 0 -err 0.05 -localrefinement 0
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi axi_ignite_0.base -fo axi_ignite -param Dhc -mo ignite -h0 50 -err 0.1 -scount 5 -paramtarget 128 -maxcount -1 -localrefinement 0
```

3. Save base flows over a range of $Re$. Dissipation from mesh coarsening is used to aid convergence at each step before refining the coarse solutions on the reference mesh.
```sh
cd $workdir && export lastfile=$(printf '%s\n' axi_ignite_*.base | sort -t_ -k3,3n | tail -1) && cd -
for Re in {1500..3200..100}; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo axi_phi_0p80_Re_"$Re" -Re "$Re" -Tr 2.3333333333333333 -Dhc 128.23469709865222 -snes_linesearch_type secant -mo phi_0p80_Re_"$Re" -err 0.04 -localrefinement 0
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi axi_phi_0p80_Re_"$Re".base -fo axi_phi_0p80_Re_"$Re" -mo phi_0p80_Re_"$Re" -err 0.02
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi axi_phi_0p80_Re_"$Re".base -fo axi_phi_0p80_Re_"$Re" -mo phi_0p80_Re_"$Re" -err 0.01 -hmax 0.2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi axi_phi_0p80_Re_"$Re".base -fo axi_phi_0p80_Re_"$Re" -mo phi_0p80_Re_"$Re" -err 0.01 -anisomax 4 -hmax 0.2 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -pv 1
export lastfile=axi_phi_0p80_Re_"$Re".base
done
```

### Global linear analysis
4. Compute global eigenspectra for axisymmetric case
```sh
for Re in {1500..3200..100}; do
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi axi_phi_0p80_Re_"$Re".base -so phi_0p80_Reswp_m_0 -eps_nev 25 -eps_target 0.1+1i -ntarget 5 -targetf 0.1+9i -eps_tol 2.22e-14
done
```

5. Compute global eigenspectra in 3D. First, link to settings and eqns files for 3D case and import axisymmetric files for 3D analysis.
```sh
ln -sf examples/schulke_etal_2026/eqns_schulke_etal_2026_3D.idp eqns.idp
ln -sf examples/schulke_etal_2026/settings_schulke_etal_2026_3D.idp settings.idp

for Re in {1500..3200..100}; do
ff-mpirun examples/schulke_etal_2026/axi_to_3D.md -v 0 -dir $workdir -fi axi_phi_0p80_Re_"$Re".base -fo 3D_phi_0p80_Re_"$Re".base
for m in {1..2}; do
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi 3D_phi_0p80_Re_"$Re".base -so phi_0p80_Reswp_m_"$m" -eps_nev 25 -eps_target 0.1+1i -ntarget 5 -targetf 0.1+9i -sym $m -eps_tol 2.22e-14
done
done
```

6. Compute leading eigenmodes associated with high- and low-frequency instabilities. First, link back to settings and eqns files for axisymmetric case.
```sh
ln -sf examples/schulke_etal_2026/eqns_schulke_etal_2026_axi.idp eqns.idp
ln -sf examples/schulke_etal_2026/settings_schulke_etal_2026_axi.idp settings.idp

ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi axi_phi_0p80_Re_3000.base -fo axi_highfreq_phi_0p80_Re_3000 -eps_target 0.01+1.6i -eps_nev 1 -strict 1 -eps_two_sided 1 -pv 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi axi_phi_0p80_Re_3100.base -fo axi_lowfreq_phi_0p80_Re_3100 -eps_target 0.01+0.8i -eps_nev 1 -strict 1 -eps_two_sided 1 -pv 1
```

7. Compute Hopf bifurcations associated with critical high- and low-frequency modes.
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi axi_highfreq_phi_0p80_Re_3000.mode -fo axi_highfreq_phi_0p80 -param Re -snes_divergence_tolerance 1e100 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi axi_highfreq_phi_0p80.hopf -fo axi_highfreq_phi_0p80 -param Re -snes_divergence_tolerance 1e100 -mo highfreq_phi_0p80 -err 0.01 -hmax 0.2 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi axi_highfreq_phi_0p80.hopf -fo axi_highfreq_phi_0p80 -param Re -snes_divergence_tolerance 1e100 -mo highfreq_phi_0p80 -err 0.01 -anisomax 4 -hmax 0.2 -nf 1 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -pv 1

ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80_Re_3100.mode -fo axi_lowfreq_phi_0p80 -param Re -snes_divergence_tolerance 1e100 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80.hopf -fo axi_lowfreq_phi_0p80 -param Re -snes_divergence_tolerance 1e100 -mo lowfreq_phi_0p80 -err 0.01 -hmax 0.2 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80.hopf -fo axi_lowfreq_phi_0p80 -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo lowfreq_phi_0p80 -err 0.01 -anisomax 4 -hmax 0.2 -nf 1 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -pv 1
```

8. Continue Hopf bifurcation curves along the $(\phi,Re)$ parameter plane.
```sh
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi axi_highfreq_phi_0p80.hopf -fo axi_highfreq_dec -mo highfreq_dec -param Re -param2 phi -nf 1 -err 0.01 -anisomax 4 -hmax 0.2 -scount 5 -maxcount -1 -paramtarget 3200 -param2target 0.65 -h0 -10 -amax 90
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi axi_highfreq_phi_0p80.hopf -fo axi_highfreq_inc -mo highfreq_inc -param Re -param2 phi -nf 1 -err 0.01 -anisomax 4 -hmax 0.2 -scount 5 -maxcount -1 -paramtarget 3200 -param2target 0.99 -h0 10 -amax 90

ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80.hopf -fo axi_lowfreq_dec -mo lowfreq_dec -param Re -param2 phi -nf 1 -err 0.01 -anisomax 4 -hmax 0.2 -scount 5 -maxcount -1 -paramtarget 3200 -param2target 0.65 -h0 -10 -amax 90
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80.hopf -fo axi_lowfreq_inc -mo lowfreq_inc -param Re -param2 phi -nf 1 -err 0.01 -anisomax 4 -hmax 0.2 -scount 5 -maxcount -1 -paramtarget 3200 -param2target 0.99 -h0 10 -amax 90
```

9. Compute double-Hopf point where high- and low-frequency modes exchange primacy
```sh
export highfreq_guess="axi_highfreq_dec_##.hopf"
export lowfreq_guess="axi_lowfreq_dec_##.hopf"

ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi $highfreq_guess -fi2 $lowfreq_guess -fo axi_highlowfreq -param Re -param2 phi -nf 0
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi axi_highlowfreq.hoho -fo axi_highlowfreq -mo highlowfreq -param Re -param2 phi -err 0.01 -hmax 0.2 -nf 0
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi axi_highlowfreq.hoho -fo axi_highlowfreq -mo highlowfreq -param Re -param2 phi -err 0.01 -anisomax 4 -hmax 0.2 -nf 1 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -pv 1
```

10. Continue along periodic orbits for high- and low-frequency bifurcations using harmonic balance
```sh
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi axi_highfreq_phi_0p80.hopf -fo axi_highfreq_phi_0p80_N_2 -param Re -Nh 2 -maxcount 10 -h0 0.01 -blocks 2 -kmax 2 -amax 90 -scount 5 -stricttangent 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi axi_highfreq_phi_0p80_N_2_10.porb -count 10 -fo axi_highfreq_phi_0p80_N_2 -param Re -Nh 2 -mo highfreq_phi_0p80_N_2 -err 0.01 -anisomax 4 -hmax 0.2 -scount 5 -h0 0.25 -blocks 2 -stricttangent 0 -fieldsplit_fieldsplit_mat_mumps_icntl_35 1 -fieldsplit_fieldsplit_mat_mumps_cntl_7 1.0e-8 -paramtarget 3200 -amax 45 -kmax 1

ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80.hopf -fo axi_lowfreq_phi_0p80_N_2 -param Re -Nh 2 -maxcount 10 -h0 0.05 -blocks 2 -kmax 2 -amax 90 -scount 5 -stricttangent 0
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80_N_2_10.porb -count 10 -fo axi_lowfreq_phi_0p80_N_2 -param Re -Nh 2 -mo lowfreq_phi_0p80_N_2 -err 0.01 -anisomax 4 -hmax 0.2 -scount 5 -h0 -0.5 -blocks 2 -stricttangent 0 -fieldsplit_fieldsplit_mat_mumps_icntl_35 1 -fieldsplit_fieldsplit_mat_mumps_cntl_7 1.0e-8 -paramtarget 2000 -amax 45 -kmax 1

cd $workdir && export lastfile=$(printf '%s\n' axi_lowfreq_phi_0p80_N_2_*.porb | sort -t_ -k7,7n | tail -1) && cd -
ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi $lastfile -fo axi_lowfreq_phi_0p80_Re_2000_N_3 -Re 2000 -Nh 3 -fieldsplit_fieldsplit_mat_mumps_icntl_35 1 -fieldsplit_fieldsplit_mat_mumps_cntl_7 1.0e-8 -blocks 3
ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80_Re_2000_N_3.porb -fo axi_lowfreq_phi_0p80_Re_2000_N_3 -param Re -Nh 3 -mo lowfreq_phi_0p80_Re_2000_N_3 -err 0.01 -anisomax 4 -hmax 0.2 -fieldsplit_fieldsplit_mat_mumps_icntl_35 1 -fieldsplit_fieldsplit_mat_mumps_cntl_7 1.0e-8 -blocks 3
```

### Nonlinear analysis
11. Compute time evolution from slightly perturbed states "above" and "below" the saddle.
```sh
ff-mpirun -np 1 examples/schulke_etal_2026/moderescale.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80_Re_2000_N_3.porb -amp 1.01 -fo axi_saddleplus_phi_0p80_Re_2000
ff-mpirun -np 1 examples/schulke_etal_2026/moderescale.md -v 0 -dir $workdir -fi axi_lowfreq_phi_0p80_Re_2000_N_3.porb -amp 0.99 -fo axi_saddleminus_phi_0p80_Re_2000
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir -fi axi_saddleplus_phi_0p80_Re_2000.base -fo axi_saddleplus_phi_0p80_Re_2000 -ts_time_step 0.005 -ts_type dirk -ts_dirk_type s7511sal -ts_adapt_type basic -scount 5 -maxcount 10000 -mo saddleplus_phi_0p80_Re_2000 -snes_atol 2.22e-14 -err 0.01 -anisomax 4 -hmax 0.2 -pv 1 -ts_adapt_dt_max 0.1 -ts_adapt_dt_min 0.005
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir  -fi axi_saddleminus_phi_0p80_Re_2000.base -fo axi_saddleminus_phi_0p80_Re_2000 -ts_time_step 0.005 -ts_type dirk -ts_dirk_type s7511sal -ts_adapt_type basic -scount 5 -maxcount 10000 -mo axi_saddleminus_phi_0p80_Re_2000 -snes_atol 2.22e-14 -err 0.01 -anisomax 4 -hmax 0.2 -pv 1 -ts_adapt_dt_max 0.1 -ts_adapt_dt_min 0.005
```
