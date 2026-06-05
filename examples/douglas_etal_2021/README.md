# Incompressible Swirling Jet Example: Douglas etal, JFM, (2021)
This file shows an example `ff-bifbox` workflow for reproducing the results of the study:
```bibtex
@article{douglas_etal_2021,
  title={Nonlinear dynamics of fully developed swirling jets},
  volume={924},
  DOI={10.1017/jfm.2021.615},
  journal={Journal of Fluid Mechanics},
  publisher={Cambridge University Press},
  author={Douglas, Christopher M. and Emerson, Benjamin L. and Lieuwen, Timothy C.},
  year={2021},
  pages={A14}
}
```
The commands below illustrate how to perform a bifurcation analysis of an incompressible swirling jet using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*} 
\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{1}{Re}\frac{\partial^2u_i}{\partial x_j^2} &= 0 \\
\frac{\partial u_i}{\partial x_i} &= 0 \\
\frac{\partial p_o}{\partial x_i}\hat{t}_i - \frac{u_{\theta}^2}{r}&= 0
\end{align*}
$$

together with the boundary conditions:

| Boundary | Constraints |
| :--- | :--- |
| Inlet, $\Gamma_i$ | $u_x=2-8r^2$, $\frac{\partial u_r}{\partial r}=0$, $u_{\theta}=2Sr$ |
| Pipe, $\Gamma_p$ | $u_x=u_r=0$, $u_{\theta}=S$ |
| Wall, $\Gamma_w$ | $u_x=u_r=u_{\theta}=p_o=0$ |
| Axis, $\Gamma_a$| $`\begin{cases}\frac{\partial u_x}{\partial r}=u_r=u_{\theta}=0, & \text{if } m=0 \\\\ u_x=\frac{\partial u_r}{\partial r}=\frac{\partial u_{\theta}}{\partial r}=0, & \text{if } \|m\|=1 \\\\ u_x=u_r=u_{\theta}=0, & \text{if } \|m\|>1\end{cases}`$ |
| Open, $\Gamma_o$ | $\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\hat{n}_j-\left(p-p_o\right)\hat{n}_i-\frac{1}{2}u_i\min\left(0,u_j\hat{n}_j\right) = 0$ |

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_p+\Gamma_w+\Gamma_a+\Gamma_o$. Solutions $\vec{q}=\left(u_i,p,p_o\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{p},\check{p}_o\right)^T$,

$$
\begin{align*} 
&\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} \\
&+ \left(\check{u}_i,p_o\hat{n}_i-\frac{1}{2}u_i\min\left(0,u_j\hat{n}_j\right)\right)_{\Gamma_o} + \left(\check{p}_o,\frac{\partial p_o}{\partial x_i}\hat{t}_i-\frac{u_{\theta}^2}{r}\right)_{\Gamma_o} = 0.
\end{align*}
$$

This weak formulation has been implemented in the equations file for this example: [eqns_douglas_etal_2021.idp](./eqns_douglas_etal_2021.idp).

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/douglas_etal_2021/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/douglas_etal_2021/eqns_douglas_etal_2021.idp eqns.idp
ln -sf examples/douglas_etal_2021/settings_douglas_etal_2021.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from `.geo` files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/douglas_etal_2021 -dir $workdir -mi swirljet.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/douglas_etal_2021/swirljet.md -mo $workdir/swirljet
```

## Perform parallel computations using `ff-bifbox`
### Steady axisymmetric dynamics
1. Compute base states on the created mesh at $Re=10$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi swirljet.msh -fo swirljet -1/Re 0.1 -S 0
```

2. Continue base state along the parameter $1/Re$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi swirljet.base -fo swirljet -param 1/Re -h0 -50 -scount 2 -maxcount -1 -paramtarget 0.01 -mo swirljet -thetamax 1
```

3. Compute base state at $Re=100$ with guess from $1/Re$ continuation
```sh
cd $workdir && export lastfile=$(printf '%s\n' swirljet_*.base | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo swirljet100 -1/Re 0.01
```

4. Continue base state at $Re=100$ along the parameter $S$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi swirljet100.base -fo swirljet100 -param S -h0 5 -scount 5 -maxcount -1 -mo swirljet100 -thetamax 1 -paramtarget 3
```

5. Compute backward and forward fold bifurcations from steady solution branch on base-adapted mesh
```sh
cd "$workdir" && set -- swirljet100_*specialpt.base && export B="$1" && export F="$2" && cd -
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $B -fo swirljet100_B -param S -mo swirljet100_B -adaptto b -thetamax 1 -nf 0
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $F -fo swirljet100_F -param S -mo swirljet100_F -adaptto b -thetamax 1 -nf 0
```

6. Adapt the mesh to the critical base/direct/adjoint solutions, save `.vtu` files for Paraview
```sh
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi swirljet100_B.fold -fo swirljet100_B -mo swirljet100_B -adaptto bda -param S -pv 1 -thetamax 1
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi swirljet100_F.fold -fo swirljet100_F -mo swirljet100_F -adaptto bda -param S -pv 1 -thetamax 1
```

7. Continue the neutral fold curve in the $(1/Re,S)$-plane with adaptive remeshing
```sh
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi swirljet100_B.fold -fo swirljet -mo swirljetfold -adaptto bda -thetamax 1 -param 1/Re -param2 S -h0 4 -scount 4 -maxcount 32
```

### Bifurcations to unsteady, 3D dynamics
8. Compute base state at $Re=133$, $S=1.8$ with guess from $Re=100$ continuation along $S$
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi swirljet100_10.base -fo swirljet1p8 -1/Re 0.0075 -S 1.8
```

9. Compute leading $|m|=1$ and $|m|=2$ eigenvalues
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi swirljet1p8.base -fo swirljet1p8m1 -eps_target 0.1-0.8i -sym -1 -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi swirljet1p8.base -fo swirljet1p8m2 -eps_target 0.1+0.4i -sym -2 -eps_pos_gen_non_hermitian
```

10. Compute Hopf bifurcation points
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi swirljet1p8m1.mode -fo swirljetm1 -param 1/Re -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi swirljet1p8m2.mode -fo swirljetm2 -param 1/Re -nf 0
```

11. Adapt the mesh to the critical base/direct/adjoint solutions, save `.vtu` files for Paraview
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi swirljetm1.hopf -fo swirljetm1 -param 1/Re -mo swirljetm1 -adaptto bda -pv 1 -thetamax 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi swirljetm2.hopf -fo swirljetm2 -param 1/Re -mo swirljetm2 -adaptto bda -pv 1 -thetamax 1
```

12. Continue the neutral Hopf curves in the $(1/Re,S)$-plane with adaptive remeshing
```sh
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi swirljetm1.hopf -fo swirljetm1 -mo swirljetm1hopf -adaptto bda -thetamax 1 -param 1/Re -param2 S -h0 4 -scount 4 -maxcount 32
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi swirljetm2.hopf -fo swirljetm2 -mo swirljetm2hopf -adaptto bda -thetamax 1 -param 1/Re -param2 S -h0 4 -scount 4 -maxcount 12
```

13. Compute the Hopf-Hopf point where the $|m|=1$ and $|m|=2$ curves cross
```sh
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi swirljetm2.hopf -fi2 swirljetm1.hopf -fo swirljetm2m1 -param 1/Re -param2 S -nf 0
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi swirljetm2m1.hoho -fo swirljetm2m1 -param 1/Re -param2 S -mo swirljetm2m1 -adaptto bda -pv 1 -thetamax 1
```

14. Compute the fold-Hopf point where the $|m|=1$ curve intersects the fold curve
```sh
cd "$workdir" && set -- swirljetm1_*specialpt.hopf && export fohoguess="$2" && cd -
ff-mpirun -np $nproc fohocompute.md -v 0 -dir $workdir -fi $fohoguess -fo swirljetm1 -param S -param2 1/Re -snes_divergence_tolerance 1e10
```

### Periodic 3D dynamics
15. Continue periodic solutions along $S$ from their initial Hopf points using the harmonic balance method with $N_h=2$.
```sh
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi swirljetm1.hopf -fo swirljetm1 -Nh 2 -mo swirljetm1porb -param S -thetamax 1 -h0 0.5 -scount 4 -maxcount -1 -paramtarget 1.9
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi swirljetm2.hopf -fo swirljetm2 -Nh 2 -mo swirljetm2porb -param S -thetamax 1 -h0 -0.5 -scount 4 -maxcount -1 -paramtarget 1.8
```
NOTE: in the actual paper, $N_h=4$ to $6$ was used to accurately resolve the periodic orbits. $N_h=2$ is used here to reduce computational cost.

16. Compute periodic solutions at $S=1.9$ ($|m|=1$) and $S=1.8$ ($|m|=2$) with $N_h=3$ using a block preconditioner.
```sh
cd $workdir && export m1file=$(printf '%s\n' swirljetm1_*.porb | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi $m1file -fo swirljetm1 -Nh 3 -S 1.9 -blocks 3
cd $workdir && export m2file=$(printf '%s\n' swirljetm2_*.porb | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc porbcompute.md -v 0 -dir $workdir -fi $m2file -fo swirljetm2 -Nh 3 -S 1.8 -blocks 3
```

### Bifurcations to aperiodic 3D dynamics
17. Compute Floquet stability of periodic solutions against each other 
```sh
ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi swirljetm1.porb -fo swirljetm1 -Nh 3 -eps_target 0.1+0.3i -sym -2 -S 1.9 -blocks 3 -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc floqcompute.md -v 0 -dir $workdir -fi swirljetm2.porb -fo swirljetm2 -Nh 3 -eps_target 0.02-0.75i -sym -1 -S 1.8 -blocks 3 -eps_pos_gen_non_hermitian
```