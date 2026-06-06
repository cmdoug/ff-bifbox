# Incompressible Swirling Annual Jet Example: Douglas etal, JFM, (2022)
This file shows an example `ff-bifbox` workflow for reproducing the results in the study:
```bibtex
@article{douglas_etal_2022,
  title={Dynamics and bifurcations of laminar annular swirling and non-swirling jets},
  volume={943},
  DOI={10.1017/jfm.2022.453},
  journal={Journal of Fluid Mechanics},
  publisher={Cambridge University Press},
  author={Douglas, Christopher M. and Emerson, Benjamin L. and Lieuwen, Timothy C.},
  year={2022},
  pages={A35}
}
```
The commands below illustrate how to perform a bifurcation analysis of an incompressible swirling annular jet using `ff-bifbox`.

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
| Inlet, $\Gamma_i$ | $`u_x=\frac{2-8r^2-\log\left(2r\right)\left(1-d^2\right)/\log\left(d\right)}{1+d^2+\left(1-d^2\right)/\log\left(d\right)}`$, $\frac{\partial u_r}{\partial r}=0$, $u_{\theta}=2Sr$ |
| Pipe, $\Gamma_p$ | $u_x=u_r=0$, $u_{\theta}=S$ |
| Wall, $\Gamma_w$ | $u_x=u_r=u_{\theta}=p_o=0$ |
| Axis, $\Gamma_a$| $`\begin{cases}\frac{\partial u_x}{\partial r}=u_r=u_{\theta}=0, & \text{if } m=0 \\\\ u_x=\frac{\partial u_r}{\partial r}=\frac{\partial u_{\theta}}{\partial r}=0, & \text{if } \|m\|=1 \\\\ u_x=u_r=u_{\theta}=0, & \text{if } \|m\|>1\end{cases}`$ |
| Open, $\Gamma_o$ | $`\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\hat{n}_j-\left(p-p_o\right)\hat{n}_i-\frac{1}{2}u_i\min\left(0,u_j\hat{n}_j\right) = 0`$ |

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the axisymmetric domain $\Omega$ with boundary $\partial\Omega=\Gamma_i+\Gamma_p+\Gamma_w+\Gamma_a+\Gamma_o$. Solutions $\vec{q}=\left(u_i,p,p_o\right)^T$ are then sought, in the appropriate spaces, such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{p},\check{p}_o\right)^T$,

$$
\begin{align*} 
&\left(\check{u}_i,\frac{\partial u_i}{\partial t} + u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial\check{u}_i}{\partial x_i},p\right)_{\Omega} + \left(\frac{\partial \check{u}_i}{\partial x_j},\frac{1}{Re}\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\check{p},\frac{\partial u_i}{\partial x_i}\right)_{\Omega} \\
&+ \left(\check{u}_i,p_o\hat{n}_i-\frac{1}{2}u_i\min\left(0,u_j\hat{n}_j\right)\right)_{\Gamma_o} + \left(\check{p}_o,\frac{\partial p_o}{\partial x_i}\hat{t}_i-\frac{u_{\theta}^2}{r}\right)_{\Gamma_o} = 0.
\end{align*}
$$

This weak formulation has been implemented in the equations file for this example: [eqns_douglas_etal_2022.idp](./eqns_douglas_etal_2022.idp).

NOTE: This code uses computational coordinates that differ from the physical coordinates by a scaling factor related to the parameter $d$ (see the `Y()` macro in [settings_douglas_etal_2022.idp](./settings_douglas_etal_2022.idp)). Note that ParaView files are exported using the physical coordinates.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/douglas_etal_2022/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/douglas_etal_2022/eqns_douglas_etal_2022.idp eqns.idp
ln -sf examples/douglas_etal_2022/settings_douglas_etal_2022.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### CASE 1: Gmsh is installed - build initial mesh directly from `.geo` files
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/douglas_etal_2022 -dir $workdir -mi annularjet.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).
#### CASE 2: Gmsh is not installed - build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/douglas_etal_2022/annularjet.md -mo $workdir/annularjet
```

## Perform parallel computations using `ff-bifbox`
### Steady axisymmetric dynamics
1. Compute base states on the created mesh at $Re=20$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi annularjet.msh -fo annularjet -1/Re 0.05 -S 0 -d 0.5
```

2. Continue base state along the parameter $1/Re$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi annularjet.base -fo annularjet -param 1/Re -h0 -100 -scount 2 -maxcount -1 -paramtarget 0.002095 -mo annularjet -thetamax 1
```

3. Compute base state at $Re=100$ with guess from $1/Re$ continuation
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi annularjet_6.base -fo annularjet100 -1/Re 0.01
```

4. Continue base state at $Re=100$ along the parameter $S$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi annularjet100.base -fo annularjet100 -param S -h0 20 -scount 5 -maxcount -1 -mo annularjet100 -thetamax 1 -paramtarget 3
```

5. Compute backward and forward fold bifurcations from steady solution branch on base-adapted mesh
```sh
cd "$workdir" && set -- annularjet100_*specialpt.base && export B="$1" && export F="$2" && cd -
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $B -fo annularjet100_B -param S -mo annularjet100_B -adaptto b -thetamax 1 -nf 0
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $F -fo annularjet100_F -param S -mo annularjet100_F -adaptto b -thetamax 1 -nf 0
```

6. Adapt the mesh to the critical base/direct/adjoint solutions, save `.vtu` files for Paraview
```sh
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi annularjet100_B.fold -fo annularjet100_B -mo annularet100_B -adaptto bda -param S -pv 1 -thetamax 1
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi annularjet100_F.fold -fo annularjet100_F -mo annularjet100_F -adaptto bda -param S -pv 1 -thetamax 1
```

7. Continue the neutral fold curve in the $(1/Re,S)$-plane and $(d,S)$-plane with adaptive remeshing
```sh
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi annularjet100_B.fold -fo annularjet_ReS -mo annularjet_ReSfold -adaptto bda -thetamax 1 -param 1/Re -param2 S -h0 4 -scount 4 -maxcount 32
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi annularjet100_B.fold -fo annularjet_dS -mo annularjet_dSfold -adaptto bda -thetamax 1 -param d -param2 S -h0 4 -scount 4 -maxcount 32
```

### Steady 3D dynamics
8. Compute base state at $Re\sim480$ with guess from $1/Re$ continuation
```sh
cd $workdir && export lastfile=$(printf '%s\n' annularjet_*.base | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo annularjet480 -1/Re 0.002095
```

9. Compute leading $|m|=1$ eigenvalue
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi annularjet480.base -fo annularjet480m1 -eps_target 0.1-0i -sym -1 -eps_pos_gen_non_hermitian
```

10. Compute zero-Hopf bifurcation point
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi annularjet480m1.mode -fo annularjetm1 -zero 1 -param 1/Re -nf 0
```

11. Adapt to zero-Hopf point and compute normal form
```sh
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi annularjetm1.hopf -fo annularjetm1 -param 1/Re -mo annularjetm1 -adaptto bda -pv 1 -thetamax 1 -zero 1
```

12. Continue the neutral zero-Hopf curve in the $(1/Re,d)$-plane with adaptive remeshing
```sh
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi annularjetm1.hopf -fo annularjetm1 -mo annularjetm1hopf -adaptto bda -thetamax 1 -param 1/Re -param2 d -h0 20 -scount 4 -maxcount 32 -zero 1
```