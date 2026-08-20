# Confined Incompressible Swirling Jet Example: Douglas & Lesshafft, JFM, (2022)
This file shows an example `ff-bifbox` workflow for reproducing the results in the study:
```bibtex
@article{douglas_lesshafft_2022,
  title={Confinement effects in laminar swirling jets},
  volume={945},
  DOI={10.1017/jfm.2022.589},
  journal={Journal of Fluid Mechanics},
  author={Douglas, Christopher M. and Lesshafft, Lutz},
  year={2022}, 
  pages={A27}
}
```
The commands below illustrate how to perform a bifurcation analysis of a confined incompressible swirling jet using `ff-bifbox`.

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
| Inlet, $\Gamma_i$ | $`u_x=2-8r^2`$, $\frac{\partial u_r}{\partial r}=0$, $u_{\theta}=2Sr$ |
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

This weak formulation has been implemented in the equations file for this example: [eqns_douglas_lesshafft_2022.idp](./eqns_douglas_lesshafft_2022.idp).

NOTE: This code uses computational coordinates that differ from the physical coordinates by scaling factor related to the parameters $L$ and $C$ (see the `X()` and `Y()` macros in [settings_douglas_lesshafft_2022.idp](./settings_douglas_lesshafft_2022.idp)). Note that ParaView files are exported using the physical coordinates.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/douglas_lesshafft_2022/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/douglas_lesshafft_2022/eqns_douglas_lesshafft_2022.idp eqns.idp
ln -sf examples/douglas_lesshafft_2022/settings_douglas_lesshafft_2022.idp settings.idp
```



## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### Build initial mesh directly from `.geo` files using Gmsh
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/douglas_lesshafft_2022 -dir $workdir -mi jet_flush_unconfined.geo
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/douglas_lesshafft_2022 -dir $workdir -mi jet_flush_confined.geo
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/douglas_lesshafft_2022 -dir $workdir -mi jet_inject_unconfined.geo
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/douglas_lesshafft_2022 -dir $workdir -mi jet_inject_confined.geo
```
Note: since no `-mo` argument is specified, the output files (`.msh`) inherit the names of their parents (`.geo`).

## Compute initial states
1. Compute base states on the created mesh at $Re=40$ from default guess
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi jet_flush_unconfined.msh -fo Re40_S0_L0_unconfined -1/Re 0.025 -S 0 -L 0 -C 1e4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi jet_flush_confined.msh -fo Re40_S0_L0_C8 -1/Re 0.025 -S 0 -L 0 -C 8
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi jet_inject_unconfined.msh -fo Re40_S0_L1_unconfined -1/Re 0.025 -S 0 -L 1 -C 1e4
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi jet_inject_confined.msh -fo Re40_S0_L2_C8 -1/Re 0.025 -S 0 -L 2 -C 8
```

2. Compute base states at $Re=100$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re40_S0_L0_unconfined.base -fo Re100_S0_L0_unconfined -mo Re100_S0_L0_unconfined -1/Re 0.01 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0_unconfined.base -fo Re100_S0_L0_unconfined -mo Re100_S0_L0_unconfined -pv 1 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re40_S0_L1_unconfined.base -fo Re100_S0_L1_unconfined -mo Re100_S0_L1_unconfined -1/Re 0.01 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L1_unconfined.base -fo Re100_S0_L1_unconfined -mo Re100_S0_L1_unconfined -pv 1 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L1_unconfined.base -fo Re100_S0_L0p2_unconfined -L 0.2 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0p2_unconfined.base -fo Re100_S0_L0p2_unconfined -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L1_unconfined.base -fo Re100_S0_L2_unconfined -L 2 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L2_unconfined.base -fo Re100_S0_L2_unconfined -thetamax 1 -hmin 1e-5

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re40_S0_L2_C8.base -fo Re100_S0_L2_C8 -mo Re100_S0_L2_C8 -1/Re 0.01 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L2_C8.base -fo Re100_S0_L2_C8 -mo Re100_S0_L2_C8 -pv 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L2_C8.base -fo Re100_S0_L2_C4 -mo Re100_S0_L2_C4 -C 4 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L2_C4.base -fo Re100_S0_L2_C4 -mo Re100_S0_L2_C4 -pv 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L2_C8.base -fo Re100_S0_L2_C16 -mo Re100_S0_L2_C16 -C 16 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L2_C16.base -fo Re100_S0_L2_C16 -mo Re100_S0_L2_C16 -pv 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L2_C16.base -fo Re100_S0_L2_C40 -mo Re100_S0_L2_C40 -C 40 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L2_C40.base -fo Re100_S0_L2_C40 -mo Re100_S0_L2_C40 -pv 1 -hmin 1e-5

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re40_S0_L0_C8.base -fo Re100_S0_L0_C8 -mo Re100_S0_L0_C8 -1/Re 0.01 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0_C8.base -fo Re100_S0_L0_C8 -mo Re100_S0_L0_C8 -pv 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0_C8.base -fo Re100_S0_L0_C4 -mo Re100_S0_L0_C4 -C 4 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0_C4.base -fo Re100_S0_L0_C4 -mo Re100_S0_L0_C4 -pv 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0_C8.base -fo Re100_S0_L0_C16 -mo Re100_S0_L0_C16 -C 16 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0_C16.base -fo Re100_S0_L0_C16 -mo Re100_S0_L0_C16 -pv 1 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0_C16.base -fo Re100_S0_L0_C40 -mo Re100_S0_L0_C40 -C 40 -hmin 1e-5
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S0_L0_C40.base -fo Re100_S0_L0_C40 -mo Re100_S0_L0_C40 -pv 1 -hmin 1e-5
```

## Computations for Figure 2 (radially unconfined case)
1. Continue base states along the parameter $S$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L0_unconfined.base -fo Re100_L0_unconfined -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L0_unconfined -paramtarget 3 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L0p2_unconfined.base -fo Re100_L0p2_unconfined -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L0p2_unconfined -paramtarget 3 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L1_unconfined.base -fo Re100_L1_unconfined -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L1_unconfined -paramtarget 3 -thetamax 1 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L2_unconfined.base -fo Re100_L2_unconfined -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L2_unconfined -paramtarget 3 -thetamax 1 -hmin 1e-5
```
2. Continue alternate branches of $L=0.2$ and $L=1$ solutions
```sh
cd $workdir && export lastfile=$(printf '%s\n' Re100_L1_unconfined_*.base | sort -t_ -k4,4n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo Re100_S3_L0p2_unconfined -L 0.2 -S 3
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S3_L0p2_unconfined.base -fo Re100_S3_L0p2_unconfined -thetamax 1 -hmin 1e-5 -mo Re100_S3_L0p2_unconfined
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S3_L0p2_unconfined.base -fo Re100_L0p2_unconfined_b2 -param S -h0 -10 -scount 5 -maxcount -1 -mo Re100_L0p2_unconfined_b2 -paramtarget 3.01 -thetamax 1 -hmin 1e-5

cd $workdir && export lastfile=$(printf '%s\n' Re100_L0p2_unconfined_*.base | sort -t_ -k4,4n | tail -n 3 | head -n 1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo Re100_S2p5_L1_unconfined -L 0.3 -S 2.5 -snes_linesearch_type secant
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S2p5_L1_unconfined.base -fo Re100_S2p5_L1_unconfined -L 0.5 -snes_linesearch_type secant
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S2p5_L1_unconfined.base -fo Re100_S2p5_L1_unconfined -L 1 -snes_linesearch_type secant
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S2p5_L1_unconfined.base -fo Re100_S2p5_L1_unconfined -thetamax 1 -hmin 1e-5 -mo Re100_S2p5_L1_unconfined
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S2p5_L1_unconfined.base -fo Re100_L1_unconfined_b2 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L1_unconfined_b2 -paramtarget 2.49 -thetamax 1 -hmin 1e-5 -snes_max_it 20
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S2p5_L1_unconfined.base -fo Re100_L1_unconfined_b3 -param S -h0 -10 -scount 5 -maxcount -1 -mo Re100_L1_unconfined_b3 -paramtarget 2.51 -thetamax 1 -hmin 1e-5 -snes_max_it 20
```

3. Compute backward and forward fold bifurcations from $L=0$ and $L=0.2$ solution branches
```sh
cd "$workdir" && set -- Re100_L0_unconfined_*specialpt.base && export B="$1" && export F="$2" && cd -
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $B -fo Re100_L0_unconfined_B -param S -mo Re100_L0_unconfined_B -thetamax 1 -hmin 1e-5 -nf 0
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi Re100_L0_unconfined_B.fold -fo Re100_L0_unconfined_B -param S -mo Re100_L0_unconfined_B -adaptto bda -thetamax 1 -hmin 1e-5 -nf 1
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $F -fo Re100_L0_unconfined_F -param S -mo Re100_L0_unconfined_F -thetamax 1 -hmin 1e-5 -nf 0
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi Re100_L0_unconfined_F.fold -fo Re100_L0_unconfined_F -param S -mo Re100_L0_unconfined_F -adaptto bda -thetamax 1 -hmin 1e-5 -nf 1

cd "$workdir" && set -- Re100_L0p2_unconfined_*specialpt.base && export B="$1" && export F="$2" && cd -
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $B -fo Re100_L0p2_unconfined_B -param S -mo Re100_L0p2_unconfined_B -thetamax 1 -hmin 1e-5 -nf 0
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi Re100_L0p2_unconfined_B.fold -fo Re100_L0p2_unconfined_B -param S -mo Re100_L0p2_unconfined_B -adaptto bda -thetamax 1 -hmin 1e-5 -nf 1 -pv 1
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $F -fo Re100_L0p2_unconfined_F -param S -mo Re100_L0p2_unconfined_F -thetamax 1 -hmin 1e-5 -nf 0
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi Re100_L0p2_unconfined_F.fold -fo Re100_L0p2_unconfined_F -param S -mo Re100_L0p2_unconfined_F -adaptto bda -thetamax 1 -hmin 1e-5 -nf 1 -pv 1
```

4. Continue the neutral fold curve in the $(S, L)$-plane with adaptive remeshing
```sh
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi Re100_L0p2_unconfined_B.fold -fo Re100_unconfined_B -mo Re100_unconfined_B -adaptto bda -thetamax 1 -hmin 1e-5 -nf 1 -param S -param2 L -h0 4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 1e-3 -amax 90
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi Re100_L0p2_unconfined_B.fold -fo Re100_unconfined_B_b2 -mo Re100_unconfined_B_b2 -adaptto bda -thetamax 1 -hmin 1e-5 -nf 1 -param S -param2 L -h0 -4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 1e-3 -amax 90
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi Re100_L0p2_unconfined_F.fold -fo Re100_unconfined_F -mo Re100_unconfined_F -adaptto bda -thetamax 1 -hmin 1e-5 -nf 1 -param S -param2 L -h0 4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 1e-3 -amax 90
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi Re100_L0p2_unconfined_F.fold -fo Re100_unconfined_F_b2 -mo Re100_unconfined_F_b2 -adaptto bda -thetamax 1 -hmin 1e-5 -nf 1 -param S -param2 L -h0 -4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 1e-3 -amax 90
```

## Computations for Figure 3 (weak axial confinement case)
1. Continue base states along the parameter $S$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L2_C4.base -fo Re100_L2_C4 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L2_C4 -paramtarget 3 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L2_C8.base -fo Re100_L2_C8 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L2_C8 -paramtarget 3 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L2_C16.base -fo Re100_L2_C16 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L2_C16 -paramtarget 3 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L2_C40.base -fo Re100_L2_C40 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L2_C40 -paramtarget 3 -hmin 1e-5
```

2. Continue alternate branch of $C=16$ solutions
```sh
cd $workdir && export lastfile=$(printf '%s\n' Re100_L2_C8_*.base | sort -t_ -k4,4n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo Re100_S3_L2_C16 -C 10 -S 3 -snes_linesearch_type secant
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S3_L2_C16.base -fo Re100_S3_L2_C16 -C 12 -snes_linesearch_type secant
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S3_L2_C16.base -fo Re100_S3_L2_C16 -C 16 -snes_linesearch_type secant
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_S3_L2_C16.base -fo Re100_S3_L2_C16 -hmin 1e-5 -mo Re100_S3_L2_C16
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S3_L2_C16.base -fo Re100_L2_C16_b2 -param S -h0 -10 -scount 5 -maxcount -1 -mo Re100_L2_C16_b2 -paramtarget 3.01 -hmin 1e-5
```

3. Compute backward and forward fold bifurcations from $C=8$ solution branch
```sh
cd "$workdir" && set -- Re100_L2_C8_*specialpt.base && export B="$1" && export F="$2" && cd -
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $B -fo Re100_L2_C8_B -param S -mo Re100_L2_C8_B -hmin 1e-5 -nf 0
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi Re100_L2_C8_B.fold -fo Re100_L2_C8_B -param S -mo Re100_L2_C8_B -adaptto bda -hmin 1e-5 -nf 1
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi $F -fo Re100_L2_C8_F -param S -mo Re100_L2_C8_F -hmin 1e-5 -nf 0
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi Re100_L2_C8_F.fold -fo Re100_L2_C8_F -param S -mo Re100_L2_C8_F -adaptto bda -hmin 1e-5 -nf 1
```

4. Continue the neutral fold curve in the $(S, C)$-plane with adaptive remeshing
```sh
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi Re100_L2_C8_B.fold -fo Re100_L2_B -mo Re100_L2_B -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 -4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 1e-3 -amax 90
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi Re100_L2_C8_B.fold -fo Re100_L2_B_b2 -mo Re100_L2_B_b2 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 1e-3 -amax 90
```

5. Compute the codimension-2 cusp bifurcation
```sh
cd "$workdir" && set -- Re100_L2_B_*specialpt.fold && export guess="$1" && cd -
ff-mpirun -np $nproc cuspcompute.md -v 0 -dir $workdir -fi $guess -fo Re100_L2 -param S -param2 C -nf 0
ff-mpirun -np $nproc cuspcompute.md -v 0 -dir $workdir -fi Re100_L2.cusp -fo Re100_L2 -param S -param2 C -mo Re100_L2 -adaptto bda -thetamax 1 -thetamax 1 -hmin 1e-5 -nf 1
```

6. Compute the Hopf bifurcations
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_L2_C4_5.base -fo Re100_L2_C4_hopfguess -S 1.2
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re100_L2_C4_hopfguess.base -fo Re100_L2_C4_hopfguess -sym -2 -eps_target 0.01+0.2i -eps_nev 1 -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_C4_hopfguess.mode -param S -fo Re100_L2_C4 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_C4.hopf -param S -fo Re100_L2_C4 -mo Re100_L2_C4 -adaptto bda -hmin 1e-5 -nf 1

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_L2_C8_5.base -fo Re100_L2_C8_hopfguess -S 1.3
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re100_L2_C8_hopfguess.base -fo Re100_L2_C8_hopfguess -sym -2 -eps_target 0.01+0.02i -eps_nev 1 -eps_pos_gen_non_hermitian
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_C8_hopfguess.mode -param S -fo Re100_L2_C8 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_C8.hopf -param S -fo Re100_L2_C8 -mo Re100_L2_C8 -adaptto bda -hmin 1e-5 -nf 1

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_L2_C16_5.base -fo Re100_L2_C16_hopfguess -S 1.5
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re100_L2_C16_hopfguess.base -fo Re100_L2_C16_hopfguess -sym -2 -eps_target 0.01+0.002i -eps_nev 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_C16_hopfguess.mode -param S -fo Re100_L2_C16 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_C16.hopf -param S -fo Re100_L2_C16 -mo Re100_L2_C16 -adaptto bda -hmin 1e-5 -nf 1

cd $workdir && export lastfile=$(printf '%s\n' Re100_L2_C16_*.base | sort -t_ -k4,4n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo Re100_L2_S3_hopfguess -S 3
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re100_L2_S3_hopfguess.base -fo Re100_L2_S3_hopfguess -sym -2 -eps_target 0.01+0.02i -eps_nev 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_S3_hopfguess.mode -param C -fo Re100_L2_S3 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_S3.hopf -param C -fo Re100_L2_S3 -mo Re100_L2_S3 -adaptto bda -hmin 1e-5 -nf 1

ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi Re100_L2_S3_hopfguess.base -fo Re100_L2_S2p7_hopfguess -S 2.7
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi Re100_L2_S2p7_hopfguess.base -fo Re100_L2_S2p7_hopfguess -sym -2 -eps_target 0.01+0.002i -eps_nev 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_S2p7_hopfguess.mode -param C -fo Re100_L2_S2p7 -nf 0 -snes_divergence_tolerance 1e10
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi Re100_L2_S2p7.hopf -param C -fo Re100_L2_S2p7 -mo Re100_L2_S2p7 -adaptto bda -hmin 1e-5 -nf 1
```

7. Continue along the Hopf neutral curves
```sh
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_C4.hopf -fo Re100_L2_C4 -mo Re100_L2_C4 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 3.99 -amax 90
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_C4.hopf -fo Re100_L2_C4_b2 -mo Re100_L2_C4_b2 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 -4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 4.01 -amax 90

ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_C8.hopf -fo Re100_L2_C8 -mo Re100_L2_C8 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 7.99 -amax 90
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_C8.hopf -fo Re100_L2_C8_b2 -mo Re100_L2_C8_b2 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 -4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 8.01 -amax 90

ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_C16.hopf -fo Re100_L2_C16 -mo Re100_L2_C16 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 15.99 -amax 90
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_C16.hopf -fo Re100_L2_C16_b2 -mo Re100_L2_C16_b2 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 -4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 16.01 -amax 90

ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_S2p7.hopf -fo Re100_L2_S2p7 -mo Re100_L2_S2p7 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 1e-3 -amax 90
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_S2p7.hopf -fo Re100_L2_S2p7_b2 -mo Re100_L2_S2p7_b2 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 -4 -scount 5 -maxcount -1 -paramtarget 3 -param2target 1e-3 -amax 90

ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi Re100_L2_S3.hopf -fo Re100_L2_S3 -mo Re100_L2_S3 -adaptto bda -hmin 1e-5 -nf 1 -param S -param2 C -h0 4 -scount 5 -maxcount -1 -paramtarget 3.01 -param2target 1e-3 -amax 90
```

## Computations for Figure 4 (strong axial confinement case)
1. Continue base states along the parameter $S$ with adaptive remeshing
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L0_C4.base -fo Re100_L0_C4 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L0_C4 -paramtarget 3 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L0_C8.base -fo Re100_L0_C8 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L0_C8 -paramtarget 3 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L0_C16.base -fo Re100_L0_C16 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L0_C16 -paramtarget 3 -hmin 1e-5
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi Re100_S0_L0_C40.base -fo Re100_L0_C40 -param S -h0 10 -scount 5 -maxcount -1 -mo Re100_L0_C40 -paramtarget 3 -hmin 1e-5
```

(Continues similarly to the above)