# Steady-state theory of thermal explosions in a homogeneous mixture

This repository reproduces the classical Frank–Kamenetskii bifurcation diagram using `ff-bifbox`.

We compute steady solutions of

$$
\Delta u + \mathrm{Da} \exp(u) = 0,
$$

in a cylindrical domain with Dirichlet boundary conditions, and perform continuation in the Damköhler number $\mathrm{Da}$.
## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Define working directory and number of processors:
```sh
export workdir=examples/FK_problem/data
export nproc=4
```
3. Create symbolic links for governing equations and solver settings.
```sh
ln -sf examples/FK_problem/eqns_FK.idp eqns.idp
ln -sf examples/FK_problem/settings_FK.idp settings.idp
````

## Build initial meshes

#### Build initial mesh using BAMG in FreeFEM
```sh
FreeFem++-mpi -v 0 examples/FK_problem/cylinder_vessel.edp -mo $workdir/cylinder_vessel
```

## Perform parallel computations using `ff-bifbox`
### Zeroth order
1. Compute base states on the created meshes at $Da=0$ from default guess
```sh
ff-mpirun -np $nproc basecompute.edp -v 0 -dir $workdir -mi cylinder_vessel.msh -fo cylinder_vessel -Da 0 -tgv -2
```

2. Continue base state along the parameter $Da$ 

```sh
ff-mpirun -np $nproc basecontinue.edp -v 0 -dir $workdir -fi cylinder_vessel.base -fo cylinder_vessel -param Da -h0 1 -scount 2 -maxcount 25 -tgv -2 -amax 10 -dmax 0.5 -kmax 1 -contorder 1
```
This step computes the steady-state bifurcation diagram.


# Compute fold bifurcation 
```sh
cd "$workdir" && declare -a foldguesslist=(cylinder_vessel_*specialpt.base) && cd -
for guess in "${foldguesslist[@]}"; do
ff-mpirun -np $nproc foldcompute.edp -v 0 -dir $workdir -fi "$guess" -fo cylinder_vessel_fold -param Da -pv 1 -tgv -2
done
````
Example output:
 alpha[Da] = 2.6252,
 beta = 0.91967,
 Da = 2.00227,
 Tmax = 1.386.

# Stability analysis 

Compute eigenvalues along the branch:
```sh
cd "$workdir" && declare -a baselist=(cylinder_vessel_*[0-9].base) && cd -
for base in "${baselist[@]}"; do
ff-mpirun -np $nproc modecompute.edp -v 0 -dir $workdir -fi "$base" -fo "$base" -eps_target 1.0+0.0i -eps_gen_hermitian -tgv -2
done
```