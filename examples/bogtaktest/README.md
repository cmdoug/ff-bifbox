# Identification of Bogdanov--Takens bifurcation point in model system

This example identifies a Bogdanov-Takens bifurcation for a reaction--diffusion equation in a cylindrical domain using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\begin{align*}
\dot{u} - v - D \nabla^2 u &= 0 \\
\dot{v} - \alpha_1 - \alpha_2 v - u^2 - uv - D \nabla^2 v &= 0
\end{align*}
$$

together with homogeneous Dirichlet boundary conditions $u=v=0$.

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the planar domain $\Omega$ with boundary $\partial\Omega$. Solutions $u,v$ are then sought, in the appropriate spaces, such that for all test functions $\check{u},\check{v}$,

$$
\left(\check{u}, \dot{u}-v\right)_{\Omega} +\left(\frac{\partial \check{u}}{\partial x_i},D\frac{\partial u}{\partial x_i}\right)_{\Omega} + \left(\check{v}, \dot{v}-\alpha_1-\alpha_2v-u^2-uv\right)_{\Omega} +\left(\frac{\partial \check{v}}{\partial x_i},D\frac{\partial v}{\partial x_i}\right)_{\Omega} = 0
$$

This weak formulation has been implemented in the equations file for this example: [eqns_BT.idp](./eqns_BT.idp).


## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Define working directory and number of processors:
```sh
export workdir=examples/bogtaktest/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/bogtaktest/eqns_BT.idp eqns.idp
ln -sf examples/bogtaktest/settings_BT.idp settings.idp
````

## Build initial meshes

#### Build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/FK_problem/vessel.md -mo $workdir/vessel
```

## Perform parallel computations using `ff-bifbox`
### Continue base state along the parameters $\alpha_1$ and $\alpha_2$

```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi vessel.msh -D 0.1 -a1 -0.1 -a2 0 -fo BT
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi BT.base -fo BTf -param a1 -h0 1 -scount 2 -maxcount 16 -amax 10 -dmax 0.5 -kmax 1 -contorder 1
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi BT.base -fo BTh -param a2 -h0 1 -scount 2 -maxcount 10 -amax 10 -dmax 0.5 -kmax 1 -contorder 1
```
This step computes the steady-state bifurcation diagram.


### Compute fold bifurcation curve
```sh
cd "$workdir" && declare -a foldguesslist=(BTf_*specialpt.base) && cd -
for guess in "${foldguesslist[@]}"; do
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi "$guess" -fo BTf -param a1 -pv 1
done
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi BTf.fold -fo BTf -param a2 -param2 a1 -h0 -1 -scount 2 -maxcount 10 -amax 90 -dmax 0.5 -kmax 1 -contorder 1
```

### Compute Hopf bifurcation curve
Compute eigenvalues along the branch:
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fo BTh_10 -fi BTh_10.base -eps_target 1.0+1.0i -eps_pos_gen_non_hermitian -eps_nev 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi BTh_10.mode -fo BTh -param a2
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi BTh.hopf -fo BTh -param a1 -param2 a2 -h0 -1 -scount 2 -maxcount 16 -amax 90 -dmax 0.5 -kmax 1 -contorder 1
```

### Identify Bogdanov-Takens point
```sh
cd "$workdir" && set -- BTh_*specialpt.hopf && export guess="$1" && cd -
ff-mpirun -np $nproc botacompute.md -v 0 -dir $workdir -fi "$guess" -fo BT -param a1 -param2 a2
```