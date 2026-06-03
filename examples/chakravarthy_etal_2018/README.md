# Buoyant Jets Chakravarthy etal. JFM (2018)
This file shows an example `ff-bifbox` workflow for reproducing the results of the paper:
```tex
@article{chakravarthy_etal_2018, 
  title={Global stability of buoyant jets and plumes}, 
  volume={835}, 
  DOI={10.1017/jfm.2017.764}, 
  journal={Journal of Fluid Mechanics}, 
  author={Chakravarthy, R. V. K. and Lesshafft, L. and Huerre, P.}, 
  year={2018}, 
  pages={654–673}
} 
```
The paper considers a buoyant jet/plume in a cylindrical (axisymmetric) domain. In the original paper, the governing equations are given as:

$$
\begin{align*} 
\rho\frac{\partial u_i}{\partial t} + \rho u_j\frac{\partial u_i}{\partial x_j} + \frac{\partial p}{\partial x_i} - \frac{1}{Re S}\left[\frac{\partial^2u_i}{\partial x_j^2} + \frac{1}{3}\delta_{ij}\frac{\partial}{\partial x_j}\left(\frac{\partial u_k}{\partial x_k}\right)\right] - \frac{Ri}{S-1}\left(1-\rho\right)\hat{e}_x &= 0 \\
\rho\frac{\partial T}{\partial t} + \rho u_i\frac{\partial T}{\partial x_i} - \frac{1}{Pr Re S}\frac{\partial^2T}{\partial x_i^2} &= 0\\
\frac{\partial\rho}{\partial t} + \frac{\partial\left(\rho u_i\right)}{\partial x_i} &= 0\\
\rho\left[1+T\left(S-1\right)\right] - 1 &= 0.
\end{align*}
$$

In the implementation, we must instead make use of the weak form. First, we eliminate the density variable using the equation of state. Then, we introduce test functions and integrate over the domain $\Omega$. We then seek, in the appropriate spaces, solutions $\vec{q}=\left(u_i,T,p\right)^T$ such that for all test functions $\vec{\check{q}}=\left(\check{u}_i,\check{T},\check{p}\right)^T$,

$$
\begin{align*} 
&\left(\check{u}_i,\frac{\partial u_i}{\partial t}\right)_{\Omega} + \left(\check{u}_i,u_j\frac{\partial u_i}{\partial x_j}\right)_{\Omega} - \left(\frac{\partial}{\partial x_i}\left(\check{u}_i\left[1+T\left(S-1\right)\right]\right),p\right)_{\Omega} - \left(\check{u}_x,Ri T\right)_{\Omega} \\
&+ \left(\frac{\partial}{\partial x_j}\left(\check{u}_i\left[1+T\left(S-1\right)\right]\right),\frac{1}{Re S}\left[\frac{\partial u_i}{\partial x_j} + \frac{1}{3}\delta_{ij}\frac{\partial u_k}{\partial x_k}\right]\right)_{\Omega} - \left(\check{u}_i\hat{n}_i,\frac{1+T\left(S-1\right)}{3 Re S}\frac{\partial u_k}{\partial x_k}\right)_{\partial\Omega}\\
&+ \left(\check{T},\frac{\partial T}{\partial t}\right)_{\Omega} + \left(\check{T},u_i\frac{\partial T}{\partial x_i}\right)_{\Omega} + \left(\frac{\partial}{\partial x_i}\left(\check{T}\left[1+T\left(S-1\right)\right]\right),\frac{1}{Pr Re S}\frac{\partial T}{\partial x_i}\right)_{\Omega} \\
&+ \left(\check{p},\frac{\partial T}{\partial t}\right)_{\Omega} - \left(\check{p},\left(\frac{1}{S-1}+T\right)\frac{\partial u_i}{\partial x_i}\right)_{\Omega} + \left(\check{p},\frac{\partial T}{\partial x_i}u_i\right)_{\Omega} = 0.
\end{align*}
$$

It is these equations that have been implemented in this work.

The commands below illustrate how to perform a weakly nonlinear analysis of the 2D incompressible flow around a cylinder and an open cavity using `ff-bifbox`.

## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/BuoyantJets/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/BuoyantJets/eqns_chakravarthy_etal_2018.idp eqns.idp
ln -sf examples/BuoyantJets/settings_chakravarthy_etal_2018.idp settings.idp
```

## Build initial meshes
`ff-bifbox` uses FreeFEM for adaptive meshing during the solution process, but it needs an initial mesh to adaptively refine.
#### build initial mesh using GMSH
```sh
FreeFem++-mpi -v 0 importgmsh.md -gmshdir examples/BuoyantJets -dir $workdir -mi jet.geo
```

#### build initial mesh using BAMG 
```sh
FreeFem++-mpi -v 0 examples/BuoyantJets/chakravarthy.md -mo $workdir/jet_0
```

## Perform parallel computations using `ff-bifbox`
### Zeroth order
1. Compute jet state at $Pr=0.7$, $S=7$, $Re=200$, and $Ri=10^{-4}$.
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi jet_0.msh -fo jet_0 -Re 1 -Pr 0.7 -Ri 1.0e-4 -S 7
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi jet_0.base -fo jet -paramtarget 200 -param Re -scount 2 -mo jet -maxcount -1
cd $workdir && export lastfile=$(printf '%s\n' jet_*.base | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo jet -Re 200 -mo jet -pv 1
```

2. Compute plume state at $Ri=10^3$
```sh
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi jet.base -fo plume -paramtarget 1e3 -param Ri -scount 5 -mo plume -maxcount -1
cd $workdir && export lastfile=$(printf '%s\n' plume_*.base | sort -t_ -k2,2n | tail -1) && cd -
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi $lastfile -fo plume -Ri 1.0e3 -mo plume -pv 1
```