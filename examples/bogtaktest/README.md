# Steady-state theory of thermal explosions in a homogeneous mixture

This example reproduces the classical Frank–Kamenetskii bifurcation diagram of Damköhler number $\mathrm{Da}$ versus temperature $u$ in a cylindrical domain using `ff-bifbox`.

In strong form, the governing equations are given as:

$$
\Delta u + \mathrm{Da} \exp(u) = 0
$$

together with homogeneous Dirichlet boundary conditions $u=0$.

The present implementation is based on a weak formulation of these equations. Test functions are introduced, and the equations are integrated over the planar domain $\Omega$ with boundary $\partial\Omega$. Solutions $u$ are then sought, in the appropriate space, such that for all test functions $\check{u}$,

$$
-\left(\frac{\partial \check{u}}{\partial x_i},\frac{\partial u}{\partial x_i}\right)_{\Omega} + \left(\check{u},\mathrm{Da}\exp\left(u\right)\right)_{\Omega} = 0
$$

This weak formulation has been implemented in the equations file for this example: [eqns_FK.idp](./eqns_FK.idp).


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
FreeFem++-mpi -v 0 examples/bogtaktest/vessel.md -mo $workdir/vessel
```

## Perform parallel computations using `ff-bifbox`
### Continue base state along the parameter $Da$ from trivial solution

```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi vessel.msh -D 1 -a1 1 -a2 1 -fo BT
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi BT.base -fo BT -param a2 -h0 -1 -scount 2 -maxcount 20 -amax 10 -dmax 0.5 -kmax 1 -contorder 1
```
This step computes the steady-state bifurcation diagram.


### Compute fold bifurcation 
```sh
cd "$workdir" && declare -a foldguesslist=(BT_*specialpt.base) && cd -
for guess in "${foldguesslist[@]}"; do
ff-mpirun -np $nproc foldcompute.md -v 0 -dir $workdir -fi "$guess" -fo BT_fold -param a2 -pv 1
done
ff-mpirun -np $nproc foldcontinue.md -v 0 -dir $workdir -fi BT_fold.fold -fo BT_fold -param a2 -param2 a1 -h0 -1 -scount 4 -adaptto bda -mo test -pv 1
````

### Stability analysis 
Compute eigenvalues along the branch:
```sh
cd "$workdir" && declare -a baselist=(FK_*[0-9].base) && cd -
for base in "${baselist[@]}"; do
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi "$base" -so FK -eps_target 1.0+0.0i -eps_gen_hermitian -nev 3
done
```