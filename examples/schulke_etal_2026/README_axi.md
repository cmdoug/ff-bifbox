# Low-Mach Inverted Conical Flame Example: Schulke, Wang & Douglas, (2026)
This file shows an example `ff-bifbox` workflow for reproducing the results in the study:
```tex
@article{schulke_etal_2026,

}
```
The commands below illustrate how to perform a bifurcation analysis of a lean premixed V-flame in an axisymmetric annular jet using `ff-bifbox`.


## Setup environment for `ff-bifbox`
1. Navigate to the main `ff-bifbox` directory.
```sh
cd ~/your/path/to/ff-bifbox/
```
2. Export working directory and number of processors for easy reference.
```sh
export workdir=examples/schulke_etal_2026/axi
export nproc=10
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
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -mi Vflame.msh -fo nonreacting_0 -Re 70 -Tr 2.3333333333333333 -Da 0 -phi 0.8 -Dhc 128.23469709865222 -snes_linesearch_type l2 -snes_rtol 0 -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting_0.base -fo nonreacting_1 -Re 350 -mo nonreacting_1 -snes_linesearch_type l2 -snes_rtol 0 -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting_1.base -fo nonreacting_2 -Re 1500 -mo nonreacting_2 -snes_linesearch_type l2 -snes_rtol 0 -pv 1
```

2. Turn on chemistry and ignite the U0 = 2.2 m/s flow at an elevated centrebody temperature and lower combustion enthalpy. Then perform continuation back to reference parameters. Coarse meshes are used for computational efficiency and stabilizing artificial dissipation. 
```sh
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi nonreacting_2.base -fo ignite_0 -Tr 3.3333333333333333 -Da 1090638.826627295 -Dhc 12.823469709865222 -mo ignite_0 -snes_rtol 0 -err 0.05 -pv 1
ff-mpirun -np $nproc basecontinue.md -v 0 -dir $workdir -fi ignite_0.base -fo ignite -param Dhc -mo ignite -h0 50 -err 0.1 -scount 5 -paramtarget 128.23469709865222 -maxcount -1 -pv 1
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignite_390.base -fo ignited -Tr 2.3333333333333333 -Dhc 128.23469709865222 -mo ignited -snes_rtol 0 -err 0.1 -snes_linesearch_type l2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi ignited.base -fo phi_0p8_Re_1500 -pv 1 -snes_rtol 0 -snes_linesearch_type l2 -mo phi_0p8_Re_1500 -err 0.01 -hmin 1e-5 -hmax 0.1 -snes_stol 0 -snes_atol 2.22e-14
```

3. Save base flows over a range of $Re$. Dissipation from mesh coarsening is used to aid convergence at each step before refining the coarse solutions on the reference mesh.
```sh
export phi=8
for Re in {1600..3200..100}; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi phi_0p"$phi"_Re_"$(($Re-100))".base -fo U0inc -Re "$Re" -snes_linesearch_type l2 -mo U0inc -err 0.05
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U0inc.base -fo U0inc -Re "$Re" -snes_linesearch_type l2 -mo U0inc -err 0.03
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi U0inc.base -fo phi_0p"$phi"_Re_"$Re".base -Re "$Re" -pv 1 -snes_rtol 0 -snes_linesearch_type l2 -mo phi_0p"$phi"_Re_"$Re" -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -snes_stol 0 -snes_atol 2.22e-14
done

export phi=8
for Re in {1600..3200..100}; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi phi_0p"$phi"_Re_"$Re".base -fo phi_0p"$phi"_Re_"$Re" -Re "$Re" -pv 1 -mi phi_0p"$phi"_Re_"$Re".msh -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14
done
```

### Global linear analysis
4. Compute global eigenspectra
```sh
export phi=8
for Re in {1600..3200..100}; do
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi phi_0p"$phi"_Re_"$Re".base -so phi_0p"$phi"_Re1600swp3200 -eps_nev 25 -eps_target 0.1+1i -ntarget 5 -targetf 0.1+9i -eps_tol 2.22e-14
done
```
5. Compute leading critical eigenmodes
```sh
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi phi_0p8_Re_2900.base -fo flametip -eps_target 0.01+1.7i -pv 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi flametip.hopf -fo flametip -eps_target 0.01+1.7i -pv 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.mode -fo flametip -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.hopf -fo flametip -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo flametip -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -nf 0

ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi phi_0p8_Re_3100.base -fo plume -eps_target 0.01+0.8i -pv 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi plume.mode -fo plume -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 0
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi plume.hopf -fo plume -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo plume -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -nf 0
```

6. Continue leading critical eigenmodes over parameter space
```sh
export decp=60
export name="fourth"
export omegaguess="1"
for dec in {61..63..1}; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi "$name"_phi_0p"$decp".hopf -fo $name -snes_divergence_tolerance 1e100 -snes_linesearch_type l2 -phi 0."$dec"5
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi "$name".base -fo "$name" -eps_target 0.01+"$omegaguess"i
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi "$name".mode -fo "$name" -param Re -snes_divergence_tolerance 1e100 -nf 0 -mo "$name" -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -nf 0 -anisomax 5 -snes_max_it 100 -phi 0."$decp"5 -bifmode eps -snes_linesearch_type l2
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi "$name".hopf -fo "$name" -param Re -snes_divergence_tolerance 1e100 -mo "$name" -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -nf 0 -anisomax 5 -snes_max_it 100 -phi 0."$decp"5 -bifmode eps -snes_linesearch_type l2
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi "$name".hopf -fo "$name" -snes_divergence_tolerance 1e100 -snes_linesearch_type l2 -phi 0."$dec"
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi "$name".base -fo "$name" -eps_target 0.01+"$omegaguess"i
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi "$name".mode -fo "$name" -param Re -snes_divergence_tolerance 1e100 -nf 0 -mo "$name" -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -nf 0 -anisomax 5 -snes_max_it 100 -phi 0."$dec" -bifmode eps -snes_linesearch_type l2
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi "$name".hopf -fo "$name" -param Re -snes_divergence_tolerance 1e100 -mo "$name" -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -nf 0 -anisomax 5 -snes_max_it 100 -phi 0."$dec" -bifmode eps -snes_linesearch_type l2
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi "$name".hopf -fo "$name"_phi_0p"$dec" -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo "$name"_phi_0p"$dec" -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -snes_max_it 100 -phi 0."$dec" -bifmode eps
export decp=$dec
done

export nump="0"
export decp=68
export num="0"
for dec in {794..700..4}; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi flametip_phi_"$nump"p"$decp".hopf -fo flametip -snes_divergence_tolerance 1e100 -snes_linesearch_type l2 -phi "$num"."$dec"
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi flametip.base -fo flametip -eps_target 0.01+1.7i
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.mode -fo flametip -param Re -snes_divergence_tolerance 1e100 -nf 0 -snes_max_it 100 -phi "$num"."$dec"
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.hopf -fo flametip -param Re -snes_divergence_tolerance 1e100 -mo flametip -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -nf 0 -anisomax 5 -snes_max_it 100 -phi "$num"."$dec"
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.hopf -fo flametip_phi_"$num"p"$dec" -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo flametip_phi_"$num"p"$dec" -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -snes_max_it 100 -snes_rtol 0 -snes_atol 2.22e-14 -phi "$num"."$dec" -nf 0
export decp=$dec
done

export num="1"
for dec in 00; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi flametip_phi_"$nump"p"$decp".hopf -fo flametip -snes_divergence_tolerance 1e100 -snes_linesearch_type l2 -phi "$num"."$dec"
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi flametip.base -fo flametip -eps_target 0.01+2.1i
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.mode -fo flametip -param Re -snes_divergence_tolerance 1e100 -nf 0 -snes_max_it 100
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.hopf -fo flametip -param Re -snes_divergence_tolerance 1e100 -mo flametip -adaptto bda -err 0.01 -hmin 1e-5 -hmax 0.2 -nf 0 -anisomax 5 -snes_max_it 100
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.hopf -fo flametip_phi_"$num"p"$dec" -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo flametip_phi_"$num"p"$dec" -adaptto bda -err 0.01 -hmin 1e-5 -hmax 0.2 -anisomax 5 -snes_max_it 100 -snes_rtol 0 -snes_atol 2.22e-14
export decp=$dec
done


for dec in {78..99..1}; do
ff-mpirun -np $nproc hopfcomputeEVP.md -v 0 -dir $workdir -fi flametip_phi_0p"$dec"0.hopf -fo flametip_phi_0p"$dec" -mo flametip_phi_0p"$dec" -param Re -pv 1 -snes_divergence_tolerance 1e100 -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -snes_max_it 100 -nf 0
done

export name="flametip"
for dec in {76..99}; do
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi "$name"_phi_0p"$dec".hopf -fo "$name"_phi_0p"$dec" -mo "$name"_phi_0p"$dec" -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 0 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5
done
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi "$name"_phi_1p00.hopf -fo "$name"_phi_1p00 -mo "$name"_phi_1p00 -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 0 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5

export name="plume"
for dec in {64..83}; do
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi "$name"_phi_0p"$dec".hopf -fo "$name"_phi_0p"$dec" -mo "$name"_phi_0p"$dec" -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 0 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5
done

export name="flametip"
for dec in {76..99}; do
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi "$name"_phi_0p"$dec".hopf -fo "$name" -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 1 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5
done
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi "$name"_phi_1p00.hopf -fo "$name" -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 1 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5

export name="plume"
for dec in {64..83}; do
ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi "$name"_phi_0p"$dec".hopf -fo "$name" -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 1 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5
done

ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi flametip_plume.hoho -fo flametip -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 1 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5

ff-mpirun -np $nproc hopfcontinue.md -v 0 -dir $workdir -fi flametip_plume.hoho -fo plume -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 1 -snes_rtol 0 -snes_stol 0 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5


export name="plume_phi_0p68"
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi $name.hopf -fo $name -mo $name -param Re -pv 1 -snes_divergence_tolerance 1e100 -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -snes_max_it 100 -bifmode eps -nf 0

export nump="0"
export decp="78"
export num="0"
for dec in 77; do
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip_phi_0p"$dec".hopf -fo flametip -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo flametip -adaptto bd -err 0.007 -hmin 1e-5 -hmax 0.2 -nf 0 -anisomax 5 -snes_max_it 100 -snes_linesearch_type l2
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi flametip.hopf -fo flametip_phi_"$num"p"$dec" -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo flametip_phi_"$num"p"$dec" -adaptto bd -err 0.007 -hmin 1e-5 -hmax 0.2 -anisomax 5 -snes_max_it 100 -snes_rtol 0 -snes_atol 2.22e-14
export decp=$dec
done

for phi in 0p74 0p76 0p78 0p80 0p82 0p84 0p86 0p88 0p90 0p92 0p94 0p96 0p98 1p00; do
ff-mpirun -np $nproc hopfcontinue.md -dir $workdir -v 0 -fi plume_phi_"$phi".hopf -fo flametip
done



export nump="0"
export decp="68"
export num="0"
for dec in {67..60..1}; do
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi plume_phi_"$nump"p"$decp".hopf -fo plume -snes_divergence_tolerance 1e100 -snes_linesearch_type l2 -phi "$num"."$dec"5
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi plume.base -fo plume -eps_target 0.01+0.9i -pv 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi plume.mode -fo plume -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 0 -snes_max_it 100 -phi "$num"."$dec"5 -bifmode eps
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi plume.hopf -fo plume -param Re -pv 1  -snes_divergence_tolerance 1e100 -mo plume -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -nf 0 -anisomax 5 -snes_max_it 100 -phi "$num"."$dec"5 -bifmode eps
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi plume.hopf -fo plumehalf -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo plumehalf -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -snes_max_it 100 -phi "$num"."$dec"5 -bifmode eps
ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi plumehalf.hopf -fo plume -snes_divergence_tolerance 1e100 -snes_linesearch_type l2 -phi "$num"."$dec"
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi plume.base -fo plume -eps_target 0.01+0.9i -pv 1
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi plume.mode -fo plume -param Re -pv 1 -snes_divergence_tolerance 1e100 -nf 0 -snes_max_it 100 -phi "$num"."$dec" -bifmode eps
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi plume.hopf -fo plume -param Re -pv 1  -snes_divergence_tolerance 1e100 -mo plume -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -nf 0 -anisomax 5 -snes_max_it 100 -phi "$num"."$dec" -bifmode eps
ff-mpirun -np $nproc hopfcompute.md -v 0 -dir $workdir -fi plume.hopf -fo plume_phi_"$num"p"$dec" -param Re -pv 1 -snes_divergence_tolerance 1e100 -mo plume_phi_"$num"p"$dec" -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -snes_max_it 100 -phi "$num"."$dec" -bifmode eps
export decp=$dec
done


ff-mpirun -np $nproc basecompute.md -v 0 -dir $workdir -fi plume_phi_0p77.base -fo flametip_plume -snes_divergence_tolerance 1e100 -snes_linesearch_type l2 -phi 0.775
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi flametip_plume.base -fo plume -eps_target 0.01+0.8i -pv 1
ff-mpirun -np $nproc modecompute.md -v 0 -dir $workdir -fi flametip_plume.base -fo flametip -eps_target 0.01+1.5i -pv 1
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi plume.mode -fi2 flametip.mode -fo flametip_plume -param Re -param2 phi -pv 1 -snes_divergence_tolerance 1e100 -nf 0 -snes_max_it 100
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi flametip_plume.hoho -fo flametip_plume -param Re -param2 phi -pv 1 -snes_divergence_tolerance 1e100 -mo flametip_plume -adaptto bda -err 0.01 -hmin 1e-5 -hmax 0.1 -nf 0 -anisomax 5 -snes_max_it 100
ff-mpirun -np $nproc hohocompute.md -v 0 -dir $workdir -fi flametip_plume.hoho -fo flametip_plume -param Re -param2 phi -pv 1 -snes_divergence_tolerance 1e100 -mo flametip_plume -adaptto bd -err 0.01 -hmin 1e-5 -hmax 0.1 -anisomax 5 -snes_max_it 100 -snes_rtol 0 -snes_atol 2.22e-14
```


7. Compute harmonic balance
```sh
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi tip_phi_0p80.hopf -fo tip_porb_phi_0p80 -param Re -adaptto 01 -Nh 2 -mo tip_porb_phi_0p80 -hmin 1e-5 -hmax 0.15 -anisomax 5 -scount 5 -h0 0.0125


ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi hi_phi_0p80.hopf -fo hi_Nh2 -param Re -adaptto 01 -Nh 2 -mo hi_Nh2 -hmin 1e-5 -hmax 0.15 -anisomax 5 -scount 5 -h0 0.25 -kmax 1 -fieldsplit_fieldsplit_mat_mumps_icntl_35 1 -fieldsplit_fieldsplit_mat_mumps_cntl_7 1.0e-8 -blocks 3
ff-mpirun -np $nproc porbcontinue.md -v 0 -dir $workdir -fi lo_phi_0p80.hopf -fo lo_Nh2 -param Re -adaptto 01 -Nh 2 -mo lo_Nh2 -hmin 1e-5 -hmax 0.15 -anisomax 5 -scount 5 -h0 0.25 -kmax 1 -fieldsplit_fieldsplit_mat_mumps_icntl_35 1 -fieldsplit_fieldsplit_mat_mumps_cntl_7 1.0e-8 -blocks 3
```


### Nonlinear analysis
8. Compute nonlinear dynamics in time domain.
```sh
ff-mpirun -np 1 examples/schulke_etal_2026/moderescale.md -v 0 -dir $workdir -fi plume_phi_0p80_Re_2000_Nh3.porb -amp 1.1 -fo plume_phi_0p80_Re_2000_plus
ff-mpirun -np 1 examples/schulke_etal_2026/moderescale.md -v 0 -dir $workdir -fi plume_phi_0p80_Re_2000_Nh3.porb -amp 0.9 -fo plume_phi_0p80_Re_2000_minus
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir -fi plume_phi_0p80_Re_2000_plus.base -fo plume_phi_0p80_Re_2000_plus -ts_time_step 0.005 -ts_type dirk -ts_dirk_type s7511sal -ts_adapt_type basic -scount 5 -maxcount 10000 -mo plume_phi_0p80_Re_2000_plus -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5 -pv 1 -ts_adapt_dt_max 0.1 -ts_adapt_dt_min 0.005
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir  -fi plume_phi_0p80_Re_2000_minus.base -fo plume_phi_0p80_Re_2000_minus -ts_time_step 0.005 -ts_type dirk -ts_dirk_type s7511sal -ts_adapt_type basic -scount 5 -maxcount 10000 -mo plume_phi_0p80_Re_2000_minus -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5 -pv 1 -ts_adapt_dt_max 0.1 -ts_adapt_dt_min 0.005
```

```sh
ff-mpirun -np 1 examples/schulke_etal_2026/moderescale.md -v 0 -dir $workdir -fi plume_phi_0p80_Re_2000_Nh3.porb -amp 1.01 -fo plume_phi_0p80_Re_2000_plus2
ff-mpirun -np 1 examples/schulke_etal_2026/moderescale.md -v 0 -dir $workdir -fi plume_phi_0p80_Re_2000_Nh3.porb -amp 0.99 -fo plume_phi_0p80_Re_2000_minus2
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir -fi plume_phi_0p80_Re_2000_plus2.base -fo plume_phi_0p80_Re_2000_plus2 -ts_time_step 0.005 -ts_type dirk -ts_dirk_type s7511sal -ts_adapt_type basic -scount 5 -maxcount 10000 -mo plume_phi_0p80_Re_2000_plus2 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5 -pv 1 -ts_adapt_dt_max 0.1 -ts_adapt_dt_min 0.005
ff-mpirun -np $nproc tdnscompute.md -v 0 -dir $workdir  -fi plume_phi_0p80_Re_2000_minus2.base -fo plume_phi_0p80_Re_2000_minus2 -ts_time_step 0.005 -ts_type dirk -ts_dirk_type s7511sal -ts_adapt_type basic -scount 5 -maxcount 10000 -mo plume_phi_0p80_Re_2000_minus2 -snes_atol 2.22e-14 -err 0.01 -hmax 0.1 -hmin 1e-5 -anisomax 5 -pv 1 -ts_adapt_dt_max 0.1 -ts_adapt_dt_min 0.005
```

Get Saddle branch
```sh

for num in {0..200..1}; do
for dec in {0..9..1}; do
ff-mpirun -np 1 examples/schulke_etal_2026/moderescale.md -v 0 -dir $workdir -fi plume_phi_0p80_Re_2000_Nh3.porb -amp 1 -timeshift "$num"."$dec" -fo plume_phi_0p80_Re_2000_saddle
done
done
```
