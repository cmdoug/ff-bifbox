# cuspcompute.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This script computes the normal form at a non-degenerate cusp point.

The normal form is written for the real amplitude $`A`$ as:

$$
\frac{dA}{dt} + \alpha_1 \cdot \delta\lambda + \alpha_2 \cdot \delta\lambda A + \beta A^3 = 0
$$

where:
- $`\alpha_1,\alpha_2`$ are the coefficients for the terms from parameter changes,
- $`\delta\lambda`$ are the parameter increments,
- $`\beta`$ is the coefficient for the term from harmonic interactions.

#### RESIDUAL EVALUATION IN MINIMALLY AUGMENTED FORMULATION
We can directly compute the residual using the varf `vR()`.

To build the augmented residual `Ra`, we must additionally compute the fold residual augmentation:

$$
g = \left\langle{}v,\mathcal{J}w\right\rangle{} = v^T\mathcal{J}w
$$

and the cusp residual augmentation:

$$
h = \langle{}v,\mathcal{H}\left(w,w\right)\rangle = v^T\mathcal{H}\left(w,w\right)
$$

where $`g`$ is the fold residual and $`v`$ and $`w`$ are the adjoint and direct eigenvectors, respectively.

$`g`$, $`v`$, and $`w`$ can be found using minimially augmented systems:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J} & \mathcal{M}p_0 \\
\left(\mathcal{M}q_0\right)^T & 0
\end{bmatrix}
\begin{bmatrix}
w \\
g
\end{bmatrix} = \begin{bmatrix}
0 \\
1
\end{bmatrix}
\end{equation}
$$

where $`q_0`$, $`p_0`$ are initial approximations of the direct & adjoint eigenvectors.

This implies:

$$
\mathcal{J}w = \mathcal{M}p_0g\qquad{}\text{and}\qquad{}(\mathcal{M}q_0)^Tw = 1
$$

so

$$
w = \mathcal{J}^{-1}\mathcal{M}p_0g\qquad{}\text{and}\qquad{}g = \frac{1}{(\mathcal{M}q_0)^T\mathcal{J}^{-1}\mathcal{M}p_0}.
$$

Note that, at $`g = 0`$, we have $`\mathcal{J}w = 0`$ and $`\left(\mathcal{M}q_0\right)^Tw = 1`$.

Similarly, we can find the adjoint eigenmode using the related system:

$$
\begin{bmatrix}
v^T & g
\end{bmatrix}\begin{bmatrix}
-\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^T & 0
\end{bmatrix} = \begin{bmatrix}
0 & 1
\end{bmatrix}
$$

This implies:

$$
v^T\mathcal{J} = g(\mathcal{M}q_0)^T\qquad{}\text{and}\qquad{}v^T\mathcal{M}p_0 = 1
$$

or, taking the transpose:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J}^T & \mathcal{M}q_0 \\
(\mathcal{M}p_0)^T & 0
\end{bmatrix}
\begin{bmatrix}
v \\
g
\end{bmatrix} = \begin{bmatrix}
0 \\
1
\end{bmatrix}
\end{equation}
$$

giving, equivalently,

$$
\mathcal{J}^Tv = \mathcal{M}q_0g\qquad{}\text{and}\qquad{}(\mathcal{M}p_0)^Tv = 1
$$

so

$$
v = \mathcal{J}^{-T}\mathcal{M}q_0g\qquad{}\text{and}\qquad{}g = \frac{1}{(\mathcal{M}p_0)^T\mathcal{J}^{-T}\mathcal{M}q_0}
$$

At $`g = 0`$, we have $`\mathcal{J}^Tv = 0`$ and $`(\mathcal{M}p_0)^Tv = 1`$, so $`v^T\mathcal{J} = 0`$ and $`v^T\mathcal{M}p_0 = 1`$.

It can then be confirmed that $`g = v^T\mathcal{J}w = w^T\mathcal{J}^Tv`$.

Finally, the augmented cusp residual can be computed directly:

$$
\begin{equation}
h = \langle{}v,\mathcal{H}\left(w,w\right)\rangle = v^T\mathcal{H}\left(w,w\right)
\end{equation}
$$

#### JACOBIAN CONSTRUCTION IN MINIMALLY AUGMENTED FORMULATION
Having computed the RHS of the augmented system in `funcRa`, we now have to build the augmented Jacobian matrix for the Newton scheme:

$$
\begin{equation}
\begin{bmatrix}
\mathcal{J} & \frac{\partial\mathcal{J}}{\partial \lambda_1} & \frac{\partial\mathcal{J}}{\partial \lambda_2} \\
(\frac{\partial{}g}{\partial q})^T& \frac{\partial{}g}{\partial\lambda_1} & \frac{\partial{}g}{\partial\lambda_2} \\
(\frac{\partial{}h}{\partial q})^T& \frac{\partial{}h}{\partial\lambda_1} & \frac{\partial{}h}{\partial\lambda_2} \\
\end{bmatrix}
\begin{bmatrix}
\delta{}q \\
\delta{}\lambda_1 \\
\delta{}\lambda_2
\end{bmatrix} = \begin{bmatrix}
\mathcal{R} \\
g \\
h
\end{bmatrix},
\end{equation}
$$

where $`g = v^T\mathcal{J}w`$ and $`h = v^T\mathcal{H}\left(w,w\right)`$.

To determine the augmented matrix entries in the second row, we differentiate Eq. (1) along each $`z`$ in $`q, \lambda_1, \lambda_2`$ to find:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^T & 0
\end{bmatrix}
\begin{bmatrix}
\frac{\partial w}{\partial z} \\
\frac{\partial g}{\partial z}
\end{bmatrix} = \begin{bmatrix}
\frac{\partial\mathcal{J}}{\partial z}w \\
0
\end{bmatrix}
\end{equation}
$$

We now left-multiply Eq. (5) by $`\begin{bmatrix}v^T & g\end{bmatrix}`$, finding due to Eq. (2) that:

$$
\frac{\partial g}{\partial z} = v^T\frac{\partial \mathcal{J}}{\partial z}w
$$

This also implies that:

$$
\mathcal{J}\frac{\partial w}{\partial z}=-\frac{\partial\mathcal{J}}{\partial z}w+\mathcal{M}p_0\frac{\partial g}{\partial z}
$$

Similarly differentiating Eq. (2), we find

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J}^T & \mathcal{M}q_0 \\
(\mathcal{M}p_0)^T & 0
\end{bmatrix}
\begin{bmatrix}
\frac{\partial v}{\partial z} \\
\left(\frac{\partial g}{\partial z}\right)^T
\end{bmatrix} = \begin{bmatrix}
\left(
\frac{\partial\mathcal{J}}{\partial z}\right)^Tv \\
0
\end{bmatrix}
\end{equation}
$$

which similarly yields after left-multiplication by $`\begin{bmatrix}w^T & g\end{bmatrix}`$ and application of Eq. (1):

$$
\left(\frac{\partial g}{\partial z}\right)^T = w^T\left(\frac{\partial \mathcal{J}}{\partial z}\right)^Tv
$$

We can also find:

$$
\mathcal{J}^T\frac{\partial v}{\partial z}=-\left(\frac{\partial\mathcal{J}}{\partial z}\right)^Tv+\mathcal{M}q_0\left(\frac{\partial g}{\partial z}\right)^T
$$

To determine the augmented matrix entries in the third row, we differentiate Eq. (3) along each $`z`$ in $`q, \lambda_1, \lambda_2`$ to find:

$$
\begin{equation}
\frac{\partial h}{\partial z} = \left(\frac{\partial v}{\partial z}\right)^T\mathcal{H}\left(w,w\right) + v^T\frac{\partial \mathcal{H}}{\partial z}\left(w,w\right) + 2v^T\mathcal{H}\left(w,\frac{\partial w}{\partial z}\right)
\end{equation}
$$

However, it is not desirable or necessary to ever construct $`\frac{\partial w}{\partial z}`$ or $`\frac{\partial v}{\partial z}`$ explicitly. Instead of computing these dense operators, we focus on their action in the associated inner products.

For the first term in Eq. (7), we have:

$$
\left(\frac{\partial v}{\partial z}\right)^T\mathcal{H}\left(w,w\right)=\left(-v^T\frac{\partial\mathcal{J}}{\partial z}+\frac{\partial g}{\partial z}\left(\mathcal{M}q_0\right)^T\right)\mathcal{J}^{-1}\mathcal{H}\left(w,w\right)=-v^T\frac{\partial\mathcal{J}}{\partial z}\hat{w}+\frac{\partial g}{\partial z}\left(\mathcal{M}q_0\right)^T\hat{w}
$$

where $`\hat{w}`$ solves the non-singular system:

$$
\begin{bmatrix}
\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^T & 0
\end{bmatrix}
\begin{bmatrix}
\hat{w} \\
h
\end{bmatrix} = \begin{bmatrix}
\mathcal{H}\left(w,w\right) \\
0
\end{bmatrix}
$$

giving equivalently,

$$
\mathcal{J}\hat{w}=\mathcal{H}\left(w,w\right)-\mathcal{M}p_0h\qquad{}\text{and}\qquad{}(\mathcal{M}q_0)^T\hat{w}=0
$$

so, using the identities derived above that $`w=\mathcal{J}^{-1}\mathcal{M}p_0g`$ and $`v=\mathcal{J}^{-T}\mathcal{M}q_0g`$,

$$
\hat{w}=\mathcal{J}^{-1}\mathcal{H}\left(w,w\right)-w\frac{h}{g}\qquad{}\text{and}\qquad{}h=v^T\mathcal{H}\left(w,w\right)
$$

Then, similarly, for the last term in Eq. (7), we have:

$$
2v^T\mathcal{H}\left(w,\frac{\partial w}{\partial z}\right)=2v^T\mathcal{H}\left(w, \cdot\right)\mathcal{J}^{-1}\left(-\frac{\partial\mathcal{J}}{\partial z}w+\mathcal{M}p_0\frac{\partial g}{\partial z}\right)=-2\hat{v}^T\frac{\partial\mathcal{J}}{\partial z}w+2\hat{v}^T\mathcal{M}p_0\frac{\partial g}{\partial z}
$$

where $`\hat{v}`$ solves the non-singular system:

$$
\begin{bmatrix}
\mathcal{J}^T & \mathcal{M}q_0 \\
(\mathcal{M}p_0)^T & 0
\end{bmatrix}
\begin{bmatrix}
\hat{v} \\
h
\end{bmatrix} = \begin{bmatrix}
\mathcal{H}\left(w,\cdot\right)^Tv \\
0
\end{bmatrix}
$$

giving equivalently,

$$
\mathcal{J}^T\hat{v}=\mathcal{H}\left(w,\cdot\right)^Tv-\mathcal{M}q_0h\qquad{}\text{and}\qquad{}(\mathcal{M}p_0)^T\hat{v}=0
$$

so, using the identities derived above that $`w=\mathcal{J}^{-1}\mathcal{M}p_0g`$ and $`v=\mathcal{J}^{-T}\mathcal{M}q_0g`$,

$$
\hat{v}=\mathcal{J}^{-T}\left(\mathcal{H}\left(w,\cdot\right)^Tv\right)-v\frac{h}{g}\qquad{}\text{and}\qquad{}h=\mathcal{H}\left(w,w\right)^Tv
$$

So we can write Eq. (4) explicitly as

$$
\begin{bmatrix}
\mathcal{J} & \frac{\partial\mathcal{J}}{\partial \lambda_1} & \frac{\partial\mathcal{J}}{\partial \lambda_2} \\
v^T\frac{\partial \mathcal{J}}{\partial q}w & v^T\frac{\partial \mathcal{J}}{\partial \lambda_1}w & v^T\frac{\partial \mathcal{J}}{\partial \lambda_2}w \\
\frac{\partial h}{\partial q} & \frac{\partial h}{\partial \lambda_1} & \frac{\partial h}{\partial \lambda_2}
\end{bmatrix}
\begin{bmatrix}
\delta{}q \\
\delta\lambda_1 \\
\delta\lambda_2
\end{bmatrix} = \begin{bmatrix}
\mathcal{R} \\
g \\
h
\end{bmatrix}
$$

where

$$
\frac{\partial h}{\partial z} = -v^T\frac{\partial\mathcal{J}}{\partial z}\hat{w} + v^T\frac{\partial \mathcal{H}}{\partial z}\left(w,w\right) - 2\hat{v}^T\frac{\partial\mathcal{J}}{\partial z}w + \left(\left(\mathcal{M}q_0\right)^T\hat{w}+2\hat{v}^T\mathcal{M}p_0\right)\frac{\partial g}{\partial z}
$$

## EXAMPLE USAGE:
### Initialize with fold guess from base file, solve on same mesh
```sh
ff-mpirun -np 4 cuspcompute.md -param <PARAM> -param2 <PARAM2> -fi <FILEIN> -fo <FILEOUT>
```

### Initialize with fold from fold file, solve on same mesh
```sh
ff-mpirun -np 4 cuspcompute.md -param <PARAM> -param2 <PARAM2> -fi <FILEIN> -fo <FILEOUT>
```

### Initialize with fold guess from file on a mesh from file
```sh
ff-mpirun -np 4 cuspcompute.md -param <PARAM> -param2 <PARAM2> -mi <MESHIN> -fi <FILEIN> -fo <FILEOUT>
```

### Initialize with fold from file, adapt mesh/solution
```sh
ff-mpirun -np 4 cuspcompute.md -param <PARAM> -param2 <PARAM2> -fi <FILEIN> -fo <FILEOUT> -mo <MESHOUT>
```

NOTE: This file should not be changed unless you know what you're doing.

SEE ALSO: [modecompute.md](./modecompute.md), [basecontinue.md](./basecontinue.md), [foldcompute.md](./foldcompute.md), [foldcontinue.md](./foldcontinue.md), [hopfcontinue.md](./hopfcontinue.md), [fohocompute.md](./fohocompute.md)

```freefem
load "iovtk"
load "PETSc"
include "settings.idp"
include "macros_bifbox.idp"
// arguments
string meshin = getARGV("-mi", "");
string meshout = getARGV("-mo", "");
string filein = getARGV("-fi", "");
string fileout = getARGV("-fo", "");
bool normalform = getARGV("-nf", 1);
bool wnlsave = getARGV("-wnl", 0);
string param = getARGV("-param", "");
string param2 = getARGV("-param2", "");
string adaptto = getARGV("-adaptto", "b");
real eps = getARGV("-eps", 1e-7);
real eps2 = getARGV("-eps2", 1e-7);
real TGV = getARGV("-tgv", -1);
string sneslinesearchtype = getARGV("-snes_linesearch_type", "none");
real[string] alpha1;
real[string] alpha2;
real beta;

// Load mesh, make FE basis
string fileroot, fileext = parsefilename(filein, fileroot); //extract file name and extension
parsefilename(fileout, fileout); // trim extension from output file, if given
if(fileext == "mode" || fileext == "resp" || fileext == "rslv" || fileext == "tdls" || fileext == "floq") {
  filein = readbasename(workdir + filein);
  fileext = parsefilename(filein, fileroot);
}
if(meshin == "") meshin = readmeshname(workdir + filein); // get mesh file
string meshroot, meshext = parsefilename(meshin, meshroot);
parsefilename(meshout, meshout); // trim extension from output mesh, if given
Th = readmeshN(workdir + meshin);
Thg = Th;
DmeshCreate(Th);
restu = restrict(XMh, XMhg, n2o);
XMh defu(ub), defu(um), defu(uma), defu(um2), defu(um3);
// Initialize solution with guess or file
if(fileext == "base") {
  ub[] = loadbase(fileroot, meshin);
}
else if(fileext == "fold") {
  ub[] = loadfold(fileroot, meshin, um[], uma[], alpha1, beta);
}
else if(fileext == "cusp") {
  ub[] = loadcusp(fileroot, meshin, um[], uma[], alpha1, alpha2, beta);
}
else if(fileext == "foho") {
  real omega, beta23, gamma22, gamma23;
  complex[string] alpha2;
  complex beta1, gamma12, gamma13;
  complex[int] q1m, q1ma;
  ub[] = loadfoho(fileroot, meshin, q1m, q1ma, um[], uma[], sym, omega, alpha2, alpha1, beta1, beta, beta23, gamma12, gamma13, gamma22, gamma23);
}
else if(fileext == "hopf") {
  real omega;
  complex[string] alpha;
  complex beta;
  complex[int] qm, qma;
  ub[] = loadhopf(fileroot, meshin, qm, qma, sym, omega, alpha, beta);
}
else if(fileext == "hoho") {
  real[int] sym1(sym.n), sym2(sym.n);
  real omega1, omega2;
  complex[string] alpha1, alpha2;
  complex beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23;
  complex[int] q1m, q1ma, q2m, q2ma;
  ub[] = loadhoho(fileroot, meshin, q1m, q1ma, q2m, q2ma, sym1, sym2, omega1, omega2, alpha1, alpha2, beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
}
else if(fileext == "tdns") {
  real time;
  ub[] = loadtdns(fileroot, meshin, time);
}
else if(fileext == "porb") {
  int Nh=1;
  real omega;
  complex[int, int] qh(um[].n, Nh);
  ub[] = loadporb(fileroot, meshin, qh, sym, omega, Nh);
}
real[int] paramvals(2);
paramvals(0) = getparam(param);
paramvals(1) = getparam(param2);
// Create distributed Mat
Mat J;
createMatu(Th, J, Pk);
// MESH ADAPTATION
bool adapt = false;
if(meshout == "") meshout = meshin; // if no adaptation
else { // if output meshfile is given, adapt mesh
  adapt = true;
  meshout = meshout + "." + meshext;
  real[int] q;
  ChangeNumbering(J, ub[], q);
  ChangeNumbering(J, ub[], q, inverse = true);
  ChangeNumbering(J, um[], q);
  ChangeNumbering(J, um[], q, inverse = true);
  ChangeNumbering(J, uma[], q);
  ChangeNumbering(J, uma[], q, inverse = true);
  XMhg defu(uG), defu(umG), defu(umaG), defu(tempu), defu(uoG); // create private global FE functions
  tempu[](restu) = ub[]; // populate local portion of global soln
  mpiAllReduce(tempu[], uG[], mpiCommWorld, mpiSUM); //aggregate local solns into global soln
  tempu[](restu) = um[]; // populate local portion of global soln
  mpiAllReduce(tempu[], umG[], mpiCommWorld, mpiSUM); //aggregate local solns into global soln
  tempu[](restu) = uma[]; // populate local portion of global soln
  mpiAllReduce(tempu[], umaG[], mpiCommWorld, mpiSUM); //aggregate local solns into global soln
  if(mpirank == 0) {  // Perform mesh adaptation (serially) on processor 0
    if(adaptto == "bo") {
      defu(tempu) = initu(defu(umaG)'*defu(umaG));
      tempu[] = sqrt(tempu[]);
      uoG[] = (umG[].*umG[]);
      uoG[] = sqrt(uoG[]);
      uoG[] .*= tempu[];
    }
    IFMACRO(dimension,2)
      if(adaptto == "b") Thg = adaptmesh(Thg, adaptu(uG), adaptmeshoptions);
      else if(adaptto == "bd") Thg = adaptmesh(Thg, adaptu(uG), adaptu(umG), adaptmeshoptions);
      else if(adaptto == "ba") Thg = adaptmesh(Thg, adaptu(uG), adaptu(umaG), adaptmeshoptions);
      else if(adaptto == "bda") Thg = adaptmesh(Thg, adaptu(uG), adaptu(umG), adaptu(umaG), adaptmeshoptions);
      else if(adaptto == "bo") Thg = adaptmesh(Thg, adaptu(uG), adaptu(uoG), adaptmeshoptions);
    ENDIFMACRO
    IFMACRO(dimension,3)
      //NOTE: 3D mesh adaptation is still under development.
      load "mshmet"
      load "mmg"
      real anisomax = getARGV("-anisomax",1.0);
      real[int] met((bool(anisomax > 1) ? 6 : 1)*Thg.nv);
      if(adaptto == "b") met = mshmet(Thg, adaptu(uG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
      else if(adaptto == "bd") met = mshmet(Thg, adaptu(uG), adaptu(umG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
      else if(adaptto == "ba") met = mshmet(Thg, adaptu(uG), adaptu(umaG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
      else if(adaptto == "bda") met = mshmet(Thg, adaptu(uG), adaptu(umG), adaptu(umaG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
      else if(adaptto == "bo") met = mshmet(Thg, adaptu(uG), adaptu(uoG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
      if(anisomax > 1.0) {
        load "aniso"
        boundaniso(6, met, anisomax);
      }
      Thg = mmg3d(Thg, metric = met, hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), hgrad = -1, verbose = verbosity-(verbosity==0));
    ENDIFMACRO
  }
  broadcast(processor(0), Thg); // broadcast global mesh to all processors
  defu(uG) = defu(uG); //interpolate global solution from old mesh to new mesh
  defu(umG) = defu(umG); //interpolate global solution from old mesh to new mesh
  defu(umaG) = defu(umaG); //interpolate global solution from old mesh to new mesh
  Th = Thg; //Reinitialize local mesh with global mesh
  Mat Adapt; // Partition new mesh and update the PETSc numbering
  createMatu(Th, Adapt, Pk);
  J = Adapt;
  defu(ub) = initu(0.0); // set local values to zero
  defu(um) = initu(0.0); // set local values to zero
  defu(uma) = initu(0.0); // set local values to zero
  defu(um2) = initu(0.0);
  defu(um3) = initu(0.0);
  restu.resize(ub[].n); // Change size of restriction operator
  restu = restrict(XMh, XMhg, n2o); // Compute new restriction from global mesh to local mesh
  ub[] = uG[](restu); //restrict global solution to each local mesh
  um[] = umG[](restu); //restrict global solution to each local mesh
  uma[] = umaG[](restu); //restrict global solution to each local mesh
}
// Build bordered block matrix from only Mat components
sym = 0;
real[int] ik(sym.n), ik2(sym.n), ik3(sym.n);
real iomega = 0.0, iomega2 = 0.0, iomega3 = 0.0;
include "eqns.idp"
Mat JlPM(J.n, mpirank == 0 ? 2 : 0), gqPM(J.n, mpirank == 0 ? 2 : 0), glPM(mpirank == 0 ? 2 : 0, mpirank == 0 ? 2 : 0); // Initialize Mat objects for bordered matrix
Mat H(J), Ja = [[J, JlPM], [gqPM', glPM]]; // make dummy Jacobian
real[int] R(ub[].n), qm(J.n), qma(J.n), qpm(J.n), qpma(J.n), pP(J.n), qP(J.n);
real h, ginv;
// FUNCTIONS
  func real[int] funcRa(real[int]& qa) {
      ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true); // PETSc to FreeFEM
      if(mpirank == 0) paramvals = qa(J.n:Ja.n-1); // Extract parameter value from state vector on proc 0
      broadcast(processor(0), paramvals);
      updateparam(param, paramvals(0));
      updateparam(param2, paramvals(1));
      R = vR(0, XMh, tgv = TGV);
      real[int] Ra;
      ChangeNumbering(J, R, Ra); // FreeFEM to PETSc
      J = vJ(XMh, XMh, tgv = -2);
      KSPSolve(J, pP, qm);
      KSPSolveTranspose(J, qP, qma);
      real ginvl = (qP'*qm);
      mpiAllReduce(ginvl, ginv, mpiCommWorld, mpiSUM);
      qm /= ginv;
      qma /= ginv;
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ChangeNumbering(J, uma[], qma, inverse = true);
      um2[] = um[];
      um3[] = vH(0, XMh, tgv = -10);
      h = J(uma[], um3[]);      
      Ra.resize(Ja.n); // Append 0 to residual vector on proc 0
      if(mpirank == 0) Ra(J.n:Ja.n-1) = [1.0/ginv, h];
      return Ra;
  }

  func int funcJa(real[int]& qa) {
      ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true); // PETSc to FreeFEM
      if(mpirank == 0) paramvals = qa(J.n:Ja.n-1); // Extract parameter value from state vector on proc 0
      broadcast(processor(0), paramvals);
      real[int] temp1, temp3;
      ChangeNumbering(J, um3[], qpm);
      KSPSolve(J, qpm, qpm);
      qpm -= h*ginv*qm;
      H = vH(XMh, XMh, tgv = 0);
      MatMultTranspose(H, qma, temp3);
      KSPSolveTranspose(J, temp3, qpma);
      qpma -= h*ginv*qma;
      qpma *= 2.0;
      updateparam(param, paramvals(0) + eps);
      updateparam(param2, paramvals(1));
      real[int] Jl = vR(0, XMh, tgv = TGV);
      real[int] Hl1 = vJ(0, XMh, tgv = -10);
      real[int] Tl1 = vH(0, XMh, tgv = -10);
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      Tl1 -= um3[];
      updateparam(param, paramvals(0));
      um3[] = vJ(0, XMh, tgv = -10);
      Tl1 += um3[];
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      Jl -= R;
      Jl /= eps;
      ChangeNumbering(J, Jl, temp1); // FreeFEM to PETSc
      um3[] = vJ(0, XMh, tgv = -10);
      Hl1 -= um3[];
      um3[] = vH(0, XMh, tgv = -10);
      Tl1 -= um3[];
      updateparam(param2, paramvals(1) + eps2);
      Jl = vR(0, XMh, tgv = TGV);
      real[int] Hl2 = vJ(0, XMh, tgv = -10);
      real[int] Tl2 = vH(0, XMh, tgv = -10);
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      Tl2 -= um3[];
      updateparam(param2, paramvals(1));
      um3[] = vJ(0, XMh, tgv = -10);
      Tl2 += um3[];
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      Jl -= R;
      Jl /= eps2;
      ChangeNumbering(J, Jl, qm); // FreeFEM to PETSc
      matrix tempPms = [[temp1, qm]]; // dense array to sparse matrix
      ChangeOperator(JlPM, tempPms, parent = Ja); // send to Mat
      um3[] = vJ(0, XMh, tgv = -10);
      Hl2 -= um3[];
      um3[] = vH(0, XMh, tgv = -10);
      Tl2 -= um3[];
      MatMultTranspose(H, qpma, qm);
      qm *= -1.0;
      IFMACRO(cubic)
      H = vT(XMh, XMh, tgv = 0);
      MatMultTranspose(H, qma, temp1);
      qm += temp1;
      ENDIFMACRO
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      H = vH(XMh, XMh, tgv = 0);
      MatMultTranspose(H, qma, temp1);
      qm -= temp1;
      real tt, ttl = (qP'*qpm) + (qpma'*pP);
      mpiAllReduce(ttl, tt, mpiCommWorld, mpiSUM);
      qm += tt*temp3;
      tempPms = [[temp3, qm]]; // dense array to sparse matrix
      ChangeOperator(gqPM, tempPms, parent = Ja); // send to Mat
      real gl1 = J(uma[], Hl1)/eps;
      real gl2 = J(uma[], Hl2)/eps2;
      ChangeNumbering(J, um3[], qpma, inverse = true);
      real hl1 = (J(uma[], Tl1) - J(um3[], Hl1))/eps + tt*gl1;
      real hl2 = (J(uma[], Tl2) - J(um3[], Hl2))/eps2 + tt*gl2;
      tempPms = [[gl1, gl2],
                 [hl1, hl2]];
      ChangeOperator(glPM, tempPms, parent = Ja); // send to Mat
      J = vJ(XMh, XMh, tgv = TGV);
      return 0;
  }
// set up Mat parameters
IFMACRO(Jprecon) Jprecon(0); ENDIFMACRO
set(Ja, sparams = "-ksp_type preonly -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_precondition full"
                + " -prefix_push fieldsplit_1_ -ksp_type preonly -pc_type redundant -redundant_pc_type lu -prefix_pop"
                + " -prefix_push fieldsplit_0_ " + KSPparams + " -prefix_pop", setup = 1);
set(J, IFMACRO(Jsetargs) Jsetargs, ENDIFMACRO prefix = "fieldsplit_0_", parent = Ja);
// Initialize
real[int] qa;
ChangeNumbering(J, ub[], qa);
qa.resize(Ja.n);
if(mpirank == 0) qa(J.n:Ja.n-1) = paramvals;
if (fileext != "fold" && fileext != "foho"){
  updateparam(param, paramvals(0) + eps);
  um2[] = vR(0, XMh);
  updateparam(param, paramvals(0));
  R = vR(0, XMh);
  um2[] -= R;
  um2[] /= eps;
  J = vJ(XMh, XMh);
  um[] = J^-1*um2[];
  uma[] = J'^-1*um2[];
}
um2[] = vM(0, XMh, tgv = 0);
ChangeNumbering(J, um[], qm);
ChangeNumbering(J, um[], qm, inverse = true);
real Mnorm = sqrt(J(um[], um2[]));
um2[] /= Mnorm;
ChangeNumbering(J, um2[], qP);
ChangeNumbering(J, uma[], qma);
ChangeNumbering(J, um[], qma, inverse = true, exchange = true);
um2[] = vM(0, XMh, tgv = 0);
ChangeNumbering(J, um[], qm, inverse = true);
um2[] *= (Mnorm/J(um[], um2[])); // so that <uma[],M*um[]> = 1
ChangeNumbering(J, um2[], pP);
// solve nonlinear problem with SNES
int ret;
SNESSolve(Ja, funcJa, funcRa, qa, reason = ret,
          sparams = "-snes_linesearch_type " + sneslinesearchtype + " -options_left no -snes_monitor -snes_converged_reason");
if (ret > 0) { // Save solution if solver converged and output file is given
  ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true);
  if(mpirank == 0) paramvals = qa(J.n:Ja.n-1);
  broadcast(processor(0), paramvals);
  updateparam(param, paramvals(0));
  updateparam(param2, paramvals(1));
  ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
  um2[] = vM(0, XMh, tgv = 0);
  ChangeNumbering(J, um[], qm, inverse = true);
  ChangeNumbering(J, uma[], qma, inverse = true);
  Mnorm = sqrt(J(um[], um2[]));
  um[] /= Mnorm; // so that <um[],M*um[]> = 1
  uma[] *= (Mnorm/J(uma[], um2[])); // so that <uma[],M*um[]> = 1
  ChangeNumbering(J, um[], qm);
  ChangeNumbering(J, uma[], qma);
  if (normalform){
    real[int,int] qDa(paramnames.n, J.n);
    Mat qPM(J.n, mpirank == 0 ? 1 : 0), pPM(J.n, mpirank == 0 ? 1 : 0); // Initialize Mat objects for bordered matrix
    Ja = [[J, qPM], [pPM', 0]]; // make dummy Jacobian
    set(Ja, sparams = "-ksp_type preonly -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_precondition full"
                    + " -prefix_push fieldsplit_1_ -ksp_type preonly -pc_type redundant -redundant_pc_type lu -prefix_pop"
                    + " -prefix_push fieldsplit_0_ " + KSPparams + " -prefix_pop", setup = 1);
    // 2nd-order
    //  A: base modification due to parameter changes
    ChangeNumbering(J, um[], qma, inverse = true, exchange = true);
    um2[] = vM(0, XMh, tgv = -10);
    ChangeNumbering(J, um2[], pP);
    matrix tempPms = [[pP]]; // dense array to sparse matrix
    ChangeOperator(pPM, tempPms, parent = Ja); // send to Mat
    ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
    um2[] = vM(0, XMh, tgv = -10);
    ChangeNumbering(J, um2[], qP);
    tempPms = [[qP]]; // dense array to sparse matrix
    ChangeOperator(qPM, tempPms, parent = Ja); // send to Mat
    J = vJ(XMh, XMh, tgv = TGV);
    ChangeNumbering(J, uma[], qma, inverse = true);
    if(paramnames[0] != ""){
      for (int k = 0; k < paramnames.n; ++k){
        real paramval = getparam(paramnames[k]);
        updateparam(paramnames[k], paramval + eps);
        um2[] = vR(0, XMh, tgv = TGV);
        updateparam(paramnames[k], paramval);
        um2[] -= R;
        um2[] /= -eps;
        ChangeNumbering(J, um2[], qP); // FreeFEM to PETSc
        qP.resize(Ja.n);
        if(mpirank == 0) qP(Ja.n-1) = 0.0;
        KSPSolve(Ja, qP, qP);
        if(mpirank == 0) alpha1[paramnames[k]] = -qP(Ja.n-1);
        broadcast(processor(0), alpha1[paramnames[k]]);
        qDa(k, :) = qP(0:J.n-1);
      }
    }
    //  B: base modifications due to quadratic nonlinear interactions
    um2[] = -0.5*um[];
    um3[] = vH(0, XMh, tgv = -10);
    ChangeNumbering(J, um3[], pP); // FreeFEM to PETSc
    pP.resize(Ja.n);
    if(mpirank == 0) pP(Ja.n-1) = 0.0;
    KSPSolve(Ja, pP, pP);
    pP.resize(J.n);
    // 3rd-order
    //  A: base modification due to parameter changes
    if(paramnames[0] != ""){
      R = vJ(0, XMh, tgv = -10);
      for (int k = 0; k < paramnames.n; ++k){
        ChangeNumbering(J, um2[], qDa(k, :), inverse = true, exchange = true);
        um3[] = vH(0, XMh, tgv = -10);
        real paramval = getparam(paramnames[k]);
        updateparam(paramnames[k], paramval + eps);
        um2[] = vJ(0, XMh, tgv = -10);
        updateparam(paramnames[k], paramval);
        um2[] -= R;
        um3[] += um2[]/eps;
        alpha2[paramnames[k]] = J(uma[], um3[]);
      }
    }
    //  B: base modification due to quadratic nonlinear interaction
    IFMACRO(cubic)
    um2[] = um[];
    um3[] = um[]/3.0;
    R = vT(0, XMh, tgv = -10);
    ENDIFMACRO
    //  C: fundamental modification due to quadratic interaction of fundamental with 2nd order modification B
    ChangeNumbering(J, um2[], pP, inverse = true, exchange = true); // FreeFEM to PETSc
    IFMACRO(cubic)
    um3[] = vH(0, XMh, tgv = -10);
    R += um3[];
    ENDIFMACRO
    IFMACRO(!cubic)
    R = vH(0, XMh, tgv = -10);
    ENDIFMACRO
    beta = 0.5*J(uma[], R);
    if(wnlsave){
      complex[int] val(1);
      XMh<complex>[int] defu(vec)(1);
      XMh<complex> defu(um);
      if(paramnames[0] != ""){
        for (int k = 0; k < paramnames.n; ++k){
          ChangeNumbering(J, um2[], qDa(k, :), inverse = true);
          vec[0][].re = um2[];
          savemode(fileout + "_wnl_param" + k, "", fileout + ".cusp", meshout, vec, val, sym, true);
        }
      }
      ChangeNumbering(J, um2[], pP, inverse = true);
      vec[0][].re = um2[];
      savemode(fileout + "_wnl_AA", "", fileout + ".cusp", meshout, vec, val, sym, true);
    }
  }
  else {
    for (int k = 0; k < paramnames.n; ++k){
      alpha1[paramnames[k]] = 0.0;
      alpha2[paramnames[k]] = 0.0;
    }
    beta = 0.0;
  }
  if(mpirank==0 && adapt) { // Save adapted mesh
    cout << "  Saving adapted mesh '" + meshout + "' in '" + workdir + "'." << endl;
    savemesh(Thg, workdir + meshout);
  }
  ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true);
  ChangeNumbering(J, um[], qm, inverse = true);
  savecusp(fileout, "", meshout, alpha1, alpha2, beta, true, true);
}
```