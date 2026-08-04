# botacompute.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This script computes the normal form at a non-degenerate Bogdanov-Takens point.

The normal form is written for the amplitude $`A`$ as:

$$
\frac{d^2 A}{dt^2}-\alpha_1\cdot\delta\lambda-\alpha_2\cdot \delta\lambda A-\beta_1A^2+\beta_2A\left(\frac{d A}{dt}\right)-\beta_3A^3-\beta_4\left(\frac{d A}{dt}\right)^2=0.
$$

This may also be equivalently expressed as the first order system:

$$
\begin{aligned}
\frac{dA}{dt}+B&=0,\\
\frac{dB}{dt}+\alpha_1\cdot\delta\lambda+\alpha_2\cdot\delta\lambda A+\beta_1\,A^2+\beta_2AB+\beta_3A^3+\beta_4B^2&=0,
\end{aligned}
$$

where:
- $`\alpha_i`$ are the coefficients for the terms from parameter changes,
- $`\delta\lambda`$ are the parameter increments,
- $`\beta_i`$ are the coefficients for the terms from nonlinear interactions.

#### RESIDUAL EVALUATION IN MINIMALLY AUGMENTED FORMULATION
We can directly compute the residual using the varf `vR()`.

To build the augmented residual `Ra`, we augment the residual with two additional functions:

$$
\begin{aligned}
g &= \langle{}v,\mathcal{J}w\rangle = v^T\mathcal{J}w,\\
h &= \langle{}v,\mathcal{M}w\rangle = v^T\mathcal{M}w,
\end{aligned}
$$

where $`g`$ and $`h`$ are scalars and $`v`$ and $`w`$ are the adjoint and direct eigenvectors, respectively.

$`g`$, $`v`$, and $`w`$ can be found using minimially augmented systems:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^T & 0
\end{bmatrix}
\begin{bmatrix}
w \\
g
\end{bmatrix} = \begin{bmatrix}
0 \\
1
\end{bmatrix},
\end{equation}
$$

where $`q_0`$, $`p_0`$ are initial approximations of the direct & adjoint eigenvectors.

This implies:

$$
\mathcal{J}w = \mathcal{M}p_0g\qquad{}\text{and}\qquad{}(\mathcal{M}q_0)^Tw = 1,
$$

so

$$
w = \mathcal{J}^{-1}\mathcal{M}p_0g\qquad{}\text{and}\qquad{}g = \frac{1}{(\mathcal{M}q_0)^T\mathcal{J}^{-1}\mathcal{M}p_0}.
$$

Note that, at $`g = 0`$, we have $`\mathcal{J}w = 0`$ and $`(\mathcal{M}q_0)^Tw = 1`$.

Similarly, we can find the adjoint eigenmode using the related system:

$$
\begin{bmatrix}
v^T & g
\end{bmatrix}\begin{bmatrix}
-\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^T & 0
\end{bmatrix} = \begin{bmatrix}
0 & 1
\end{bmatrix}.
$$

This implies:

$$
v^T\mathcal{J} = g(\mathcal{M}q_0)^T\qquad{}\text{and}\qquad{}v^T\mathcal{M}p_0 = 1,
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
\end{bmatrix},
\end{equation}
$$

giving, equivalently,

$$
\mathcal{J}^Tv = \mathcal{M}q_0g\qquad{}\text{and}\qquad{}(\mathcal{M}p_0)^Tv = 1,
$$

so

$$
v = \mathcal{J}^{-T}\mathcal{M}q_0g\qquad{}\text{and}\qquad{}g = \frac{1}{(\mathcal{M}p_0)^T\mathcal{J}^{-T}\mathcal{M}q_0}.
$$

At $`g = 0`$, we have $`\mathcal{J}^Tv = 0`$ and $`(\mathcal{M}p_0)^Tv = 1`$, so $`v^T\mathcal{J} = 0`$ and $`v^T\mathcal{M}p_0 = 1`$.

Finally, the augmented Bogdanov-Takens residual can be computed directly:

$$
\begin{equation}
h = \langle{}v,\mathcal{M}w\rangle = v^T\mathcal{M}w.
\end{equation}
$$

#### JACOBIAN CONSTRUCTION IN MINIMALLY AUGMENTED FORMULATION
Having computed the RHS of the augmented system in `funcRa`, we now have to build the augmented Jacobian matrix for the Newton scheme:

$$
\begin{equation}
\begin{bmatrix}
\mathcal{J} & \frac{\partial\mathcal{J}}{\partial \lambda_1} & \frac{\partial\mathcal{J}}{\partial \lambda_2} \\
(\frac{\partial{}g}{\partial q})^T& \frac{\partial{}g}{\partial\lambda_1} & \frac{\partial{}g}{\partial \lambda_2} \\
(\frac{\partial{}h}{\partial q})^T& \frac{\partial{}h}{\partial\lambda_1} & \frac{\partial{}h}{\partial \lambda_2}
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

where $`g = v^T\mathcal{J}w`$ and $`h = v^T\mathcal{M}w`$.

To determine the matrix entries, we differentiate Eq. (1) along each $`z`$ in $`q, \lambda_1, \lambda_2`$ to find:

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

We now left-multiply Eq. (4) by $\begin{bmatrix}v^T & g\end{bmatrix}$, finding due to Eq. (2) that:

$$
\frac{\partial g}{\partial z} = v^T\frac{\partial \mathcal{J}}{\partial z}w
$$

This also implies that:

$$
\mathcal{J}\frac{\partial w}{\partial z}=-\frac{\partial\mathcal{J}}{\partial z}w+\mathcal{M}p_0\frac{\partial g}{\partial z}
$$

Similarly differentiating Eq. (2), we find:

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
\left(\frac{\partial\mathcal{J}}{\partial z}\right)^Tv \\
0
\end{bmatrix}
\end{equation}
$$

which similarly yields after left-multiplying by $\left[w^T\quad{}g\right]$ and application of Eq. (1):

$$
\left(\frac{\partial g}{\partial z}\right)^T = w^T\left(\frac{\partial \mathcal{J}}{\partial z}\right)^Tv
$$

This also implies that:

$$
\mathcal{J}^T\frac{\partial v}{\partial z}=-\left(\frac{\partial\mathcal{J}}{\partial z}\right)^Tv+\mathcal{M}q_0\left(\frac{\partial g}{\partial z}\right)^T
$$

To determine the augmented matrix entries in the third row, we differentiate Eq. (3) along each $`z`$ in $`q, \lambda_1, \lambda_2`$ to find:

$$
\begin{equation}
\frac{\partial h}{\partial z} = \left(\frac{\partial v}{\partial z}\right)^T\mathcal{M}w + v^T\frac{\partial \mathcal{M}}{\partial z}w + v^T\mathcal{M}\frac{\partial w}{\partial z}
\end{equation}
$$

However, it is not desirable or necessary to ever construct $`\frac{\partial w}{\partial z}`$ or $`\frac{\partial v}{\partial z}`$ explicitly. Instead of computing these dense operators, we focus on their action in the associated inner products.

For the first term in Eq. (7), we have:

$$
\left(\frac{\partial v}{\partial z}\right)^T\mathcal{M}w=\left(-v^T\frac{\partial\mathcal{J}}{\partial z}+\frac{\partial g}{\partial z}\left(\mathcal{M}q_0\right)^T\right)\mathcal{J}^{-1}\mathcal{M}w=-v^T\frac{\partial\mathcal{J}}{\partial z}\hat{w}+\frac{\partial g}{\partial z}\left(\mathcal{M}q_0\right)^T\hat{w}
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
\mathcal{M}w \\
0
\end{bmatrix}
$$

giving equivalently,

$$
\mathcal{J}\hat{w}=\mathcal{M}w-\mathcal{M}p_0h\qquad{}\text{and}\qquad{}(\mathcal{M}q_0)^T\hat{w}=0
$$

so, using the identities derived above that $`w=\mathcal{J}^{-1}\mathcal{M}p_0g`$ and $`v=\mathcal{J}^{-T}\mathcal{M}q_0g`$,

$$
\hat{w}=\mathcal{J}^{-1}\mathcal{M}w-w\frac{h}{g}\qquad{}\text{and}\qquad{}h=v^T\mathcal{M}w
$$

Then, similarly, for the last term in Eq. (7), we have:

$$
v^T\mathcal{M}\frac{\partial w}{\partial z}=v^T\mathcal{M}\mathcal{J}^{-1}\left(-\frac{\partial\mathcal{J}}{\partial z}w+\mathcal{M}p_0\frac{\partial g}{\partial z}\right)=-\hat{v}^T\frac{\partial\mathcal{J}}{\partial z}w+\hat{v}^T\mathcal{M}p_0\frac{\partial g}{\partial z}
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
\mathcal{M}^Tv \\
0
\end{bmatrix}
$$

giving equivalently,

$$
\mathcal{J}^T\hat{v}=\mathcal{M}^Tv-\mathcal{M}q_0h\qquad{}\text{and}\qquad{}(\mathcal{M}p_0)^T\hat{v}=0
$$

so, using the identities derived above that $`w=\mathcal{J}^{-1}\mathcal{M}p_0g`$ and $`v=\mathcal{J}^{-T}\mathcal{M}q_0g`$,

$$
\hat{v}=\mathcal{J}^{-T}\mathcal{M}^Tv-v\frac{h}{g}\qquad{}\text{and}\qquad{}h=w^T\mathcal{M}^Tv
$$

So we can write Eq. (3) explicitly as

$$
\begin{bmatrix}
\mathcal{J} & \frac{\partial\mathcal{J}}{\partial \lambda_1} & \frac{\partial\mathcal{J}}{\partial \lambda_2} \\
\Re\left(v^T\frac{\partial \mathcal{J}}{\partial q}w\right) & \Re\left(v^T\frac{\partial \mathcal{J}}{\partial \lambda_1}w\right) & \Re\left(v^T\frac{\partial \mathcal{J}}{\partial \lambda_2}w\right) \\
\Re\left(\frac{\partial h}{\partial q}\right) & \Re\left(\frac{\partial h}{\partial \lambda_1}\right) & \Re\left(\frac{\partial h}{\partial \lambda_1}\right)
\end{bmatrix}
\begin{bmatrix}
\delta{}q \\
\delta\lambda_1 \\
\delta\lambda_2
\end{bmatrix} = \begin{bmatrix}
\mathcal{R} \\
\Re(g) \\
\Re(h)
\end{bmatrix}
$$

where

$$
\frac{\partial h}{\partial z} = -v^T\frac{\partial\mathcal{J}}{\partial z}\hat{w} + v^T\frac{\partial \mathcal{M}}{\partial z}w - \hat{v}T\frac{\partial\mathcal{J}}{\partial z}w + \left(\left(\mathcal{M}q_0\right)^T\hat{w}+\hat{v}^T\mathcal{M}p_0\right)\frac{\partial g}{\partial z}
$$

## EXAMPLE USAGE:
### Initialize with Bogdanov-Takens guess from base file, solve on same mesh
```sh
ff-mpirun -np 4 botacompute.md -param <PARAM> -fi <FILEIN> -bfi <BASEFILEIN> -fo <FILEOUT>
```

### Initialize with Bogdanov-Takens from base and mode file, solve on same mesh
```sh
ff-mpirun -np 4 botacompute.md -param <PARAM> -fi <FILEIN> -fo <FILEOUT>
```

### Initialize with Bogdanov-Takens guess from file on a mesh from file
```sh
ff-mpirun -np 4 botacompute.md -param <PARAM> -mi <MESHIN> -bfi <BASEFILEIN> -fi <FILEIN> -fo <FILEOUT>
```

### Initialize with Bogdanov-Takens from file, adapt mesh/solution
```sh
ff-mpirun -np 4 botacompute.md -param <PARAM> -fi <FILEIN> -fo <FILEOUT> -mo <MESHOUT>
```

NOTE: This file should not be changed unless you know what you're doing.

SEE ALSO: [modecompute.md](./modecompute.md), [hopfcompute.md](./hopfcompute.md), [hopfcontinue.md](./hopfcontinue.md), [fohocompute.md](./fohocompute.md), [foldcompute.md](./foldcompute.md), [foldcontinue.md](./foldcontinue.md), [porbcontinue.md](./porbcontinue.md)

```freefem
load "iovtk"
load "PETSc"
include "settings.idp"
include "macros_bifbox.idp"
// arguments
string meshin = getARGV("-mi", "");
string meshout = getARGV("-mo", "");
string filein = getARGV("-fi", "");
string basefilein = getARGV("-bfi", "");
string fileout = getARGV("-fo", "");
bool normalform = getARGV("-nf", 1);
bool wnlsave = getARGV("-wnl", 0);
int select = getARGV("-select", 1);
string param = getARGV("-param", "");
string param2 = getARGV("-param2", "");
string adaptto = getARGV("-adaptto", "b");
real eps = getARGV("-eps", 1e-7);
real eps2 = getARGV("-eps2", 1e-7);
string sneslinesearchtype = getARGV("-snes_linesearch_type", "none");
real TGV = getARGV("-tgv", -1);
real[string] alpha1, alpha2;
real beta1, beta2, beta3, beta4;

// Load mesh, make FE basis
string fileroot, fileext = parsefilename(filein, fileroot); //extract file name and extension
parsefilename(fileout, fileout); // trim extension from output file, if given
if((fileext == "mode" || fileext == "resp" || fileext == "rslv" || fileext == "tdls" || fileext == "floq") && basefilein == "") basefilein = readbasename(workdir + filein);
string basefileroot, basefileext = parsefilename(basefilein, basefileroot);
if(meshin == "") meshin = readmeshname(workdir + filein); // get mesh file
string meshroot, meshext = parsefilename(meshin, meshroot);
parsefilename(meshout, meshout); // trim extension from output mesh, if given
Th = readmeshN(workdir + meshin);
Thg = Th;
DmeshCreate(Th);
restu = restrict(XMh, XMhg, n2o);
XMh defu(ub), defu(um), defu(uma), defu(um2), defu(um3);
if (fileext == "bota") {
  ub[] = loadbota(fileroot, meshin, um[], uma[], alpha1, alpha2, beta1, beta2, beta3, beta4);
}
else if (fileext == "hopf") {
  real omega;
  complex[string] alpha;
  complex beta;
  real[int] sym1(sym.n);
  complex[int] qm(um[].n), qma(um[].n);
  ub[] = loadhopf(fileroot, meshin, qm, qma, sym1, omega, alpha, beta);
  um[] = qm.re;
  uma[] = qma.re;
}
else if (fileext == "fold") {
  real[string] alpha; 
  real beta;
  ub[] = loadfold(fileroot, meshin, um[], uma[], alpha, beta);
}
else if (fileext == "foho") {
  real omega;
  complex[string] alpha1;
  real[string] alpha2;
  real[int] sym1(sym.n);
  real beta22, beta23, gamma22, gamma23;
  complex beta1, gamma12, gamma13;
  complex[int] q1m(um[].n), q1ma(um[].n);
  if(select == 1){
    ub[] = loadfoho(fileroot, meshin, q1m, q1ma, um[], uma[], sym1, omega, alpha1, alpha2, beta1, beta22, beta23, gamma12, gamma13, gamma22, gamma23);
    um[] = q1m.re;
    uma[] = q1ma.re;
  }
  else if(select == 2){
    ub[] = loadfoho(fileroot, meshin, q1m, q1ma, um[], uma[], sym1, omega, alpha1, alpha2, beta1, beta22, beta23, gamma12, gamma13, gamma22, gamma23);
  }
}
else if(fileext == "hoho") {
  real omega, omegaN;
  complex[string] alpha1, alpha2;
  real[int] sym1(sym.n), symN(sym.n);
  complex gamma11, gamma12, gamma13, gamma21, gamma22, gamma23, beta1, beta2;
  complex[int] q1m(um[].n), q1ma(um[].n), qNm, qNma;
  if(select == 1){
    ub[] = loadhoho(fileroot, meshin, q1m, q1ma, qNm, qNma, sym1, symN, omega, omegaN, alpha1, alpha2, beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
  }
  else if(select == 2){
    ub[] = loadhoho(fileroot, meshin, qNm, qNma, q1m, q1ma, symN, sym1, omegaN, omega, alpha2, alpha1, beta2, beta1, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
  }
  um[] = q1m.re;
  uma[] = q1ma.re;
}
else if (fileext == "mode") {
  complex eigenvalue;
  real[int] sym1(sym.n);
  complex[int] q1m = loadmode(fileroot, meshin, sym1, eigenvalue);
  um[] = q1m.re;
}
else if (fileext == "resp") {
  real omega;
  real[int] sym1(sym.n);
  complex[int] q1m = loadresp(fileroot, meshin, sym1, omega);
  um[] = q1m.re;
}
else if (fileext == "rslv") {
  real omega, gain;
  real[int] sym1(sym.n);
  complex[int] fm;
  complex[int] q1m = loadrslv(fileroot, meshin, fm, sym1, omega, gain);
  um[] = q1m.re;
}
else if(fileext == "porb") {
  int Nh=1;
  real omega;
  real[int] sym1(sym.n);
  complex[int, int] qh(um[].n, Nh);
  ub[] = loadporb(fileroot, meshin, qh, sym1, omega, Nh);
  um[] = qh(:, 0).re;
}
else if(fileext == "floq") {
  int Nh=1;
  real omega;
  real[int] sym1(sym.n);
  complex[int, int] qh(um[].n, 2);
  complex eigenvalue;
  real[int] symtemp(sym.n);
  complex[int] q1m = loadfloq(fileroot, meshin, qh, sym1, eigenvalue, symtemp, omega, Nh);
  um[] = q1m.re;
}
else assert(false); // invalid input filetype
if (basefileext == "base") {
  ub[] = loadbase(basefileroot, meshin);
}
else if(basefileext == "fold") {
  real[string] alpha;
  real beta;
  real[int] qm, qma;
  ub[] = loadfold(basefileroot, meshin, qm, qma, alpha, beta);
}
else if(basefileext == "cusp") {
  real[string] alpha, alphaR;
  real beta;
  real[int] qm, qma;
  ub[] = loadcusp(basefileroot, meshin, qm, qma, alpha, alphaR, beta);
}
else if(basefileext == "hopf") {
  real omega;
  complex[string] alpha;
  complex beta;
  complex[int] qm, qma;
  ub[] = loadhopf(basefileroot, meshin, qm, qma, sym, omega, alpha, beta);
}
else if(basefileext == "bota") {
  real[string] alpha1, alpha2;
  real beta1, beta2;
  real[int] qm, qma;
  ub[] = loadbota(basefileroot, meshin, qm, qma, alpha1, alpha2, beta1, beta2, beta3, beta4);
}
else if(basefileext == "foho") {
  real omega;
  complex[string] alpha1;
  complex beta1, gamma12, gamma13;
  real[string] alpha2;
  real beta22, beta23, gamma22, gamma23;
  complex[int] q1m, q1ma;
  real[int] q2m, q2ma;
  ub[] = loadfoho(basefileroot, meshin, q1m, q1ma, q2m, q2ma, sym, omega, alpha1, alpha2, beta1, beta22, beta23, gamma12, gamma13, gamma22, gamma23);
}
else if(basefileext == "hoho") {
  real[int] sym2(sym.n);
  real omega1, omega2;
  complex[string] alpha1, alpha2;
  complex beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23;
  complex[int] q1m, q1ma, q2m, q2ma;
  ub[] = loadhoho(basefileroot, meshin, q1m, q1ma, q2m, q2ma, sym, sym2, omega1, omega2, alpha1, alpha2, beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
}
else if(basefileext == "tdns") {
  real time;
  ub[] = loadtdns(basefileroot, meshin, time);
}
else if(basefileext == "porb") {
  int Nh=1;
  real omega;
  complex[int, int] qh(um[].n, Nh);
  ub[] = loadporb(basefileroot, meshin, qh, sym, omega, Nh);
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
  XMhg defu(uG), defu(umG), defu(umaG), defu(tempu), defu(uoG);
  tempu[](restu) = ub[]; // populate local portion of global soln
  mpiAllReduce(tempu[], uG[], mpiCommWorld, mpiSUM);
  tempu[](restu) = um[]; // populate local portion of global soln
  mpiAllReduce(tempu[], umG[], mpiCommWorld, mpiSUM);
  tempu[](restu) = uma[]; // populate local portion of global soln
  mpiAllReduce(tempu[], umaG[], mpiCommWorld, mpiSUM);
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
  broadcast(processor(0), Thg);
  defu(uG) = defu(uG);
  defu(umG) = defu(umG);
  defu(umaG) = defu(umaG);
  Th = Thg;
  Mat Adapt;
  createMatu(Th, Adapt, Pk);
  J = Adapt;
  defu(ub) = initu(0.0);
  defu(um) = initu(0.0);
  defu(uma) = initu(0.0);
  defu(um2) = initu(0.0);
  defu(um3) = initu(0.0);
  restu.resize(ub[].n); // Change size of restriction operator
  restu = restrict(XMh, XMhg, n2o); // Compute new restriction from global mesh to local mesh
  ub[] = uG[](restu);
  um[] = umG[](restu);
  uma[] = umaG[](restu);
}

// Build bordered block matrix from only Mat components
sym = 0;
real[int] ik(sym.n), ik2(sym.n), ik3(sym.n);
real iomega = 0.0, iomega2 = 0.0, iomega3 = 0.0;
include "eqns.idp"
Mat JlPM(J.n, mpirank == 0 ? 2 : 0), gqPM(J.n, mpirank == 0 ? 2 : 0), glPM(mpirank == 0 ? 2 : 0, mpirank == 0 ? 2 : 0); // Initialize Mat objects for bordered matrix
Mat M(J), Ja = [[J, JlPM], [gqPM', glPM]]; // make dummy Jacobian
real[int] R(ub[].n), qm(J.n), qma(J.n), qpm(J.n), qpma(J.n), pP(J.n), qP(J.n);
// FUNCTIONS
  func real[int] funcRa(real[int]& qa) {
      ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true); // PETSc to FreeFEM
      if(mpirank == 0) paramvals = qa(J.n:Ja.n-1);
      broadcast(processor(0), paramvals);
      updateparam(param, paramvals(0));
      updateparam(param2, paramvals(1));
      R = vR(0, XMh, tgv = TGV);
      real[int] Ra;
      J = vJ(XMh, XMh, tgv = -2);
      KSPSolve(J, pP, qm);
      KSPSolveTranspose(J, qP, qma);
      real h, ginv, ginvl = (qP'*qm);
      mpiAllReduce(ginvl, ginv, mpiCommWorld, mpiSUM);
      qm /= ginv; // rescale direct mode
      qma /= ginv; // rescale adjoint mode
      M = vM(XMh, XMh, tgv = 0);
      MatMult(M, qm, Ra);
      ginvl = (qma'*Ra);
      mpiAllReduce(ginvl, h, mpiCommWorld, mpiSUM);
      KSPSolve(J, Ra, qpm);
      qpm -= h*ginv*qm;
      MatMultTranspose(M, qma, Ra);
      KSPSolveTranspose(J, Ra, qpma);
      qpma -= h*ginv*qma;
      ChangeNumbering(J, R, Ra); // FreeFEM to PETSc
      Ra.resize(Ja.n); // Append 0 to residual vector on proc 0
      if(mpirank == 0) Ra(J.n:Ja.n-1) = [1.0/ginv, h];
      return Ra;
  }

  func int funcJa(real[int]& qa) {
      ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true); // PETSc to FreeFEM
      if(mpirank == 0) paramvals = qa(J.n:Ja.n-1);
      broadcast(processor(0), paramvals);
      real[int] temp1, temp3;
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ChangeNumbering(J, uma[], qma, inverse = true);
      updateparam(param, paramvals(0) + eps);
      um3[] = vR(0, XMh, tgv = TGV);
      um3[] -= R;
      um3[] /= eps;
      ChangeNumbering(J, um3[], temp1);
      real[int] Jl1 = vJ(0, XMh, tgv = -10);
      real[int] Ml1 = vM(0, XMh, tgv = -10);
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      Ml1 -= um3[];
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      updateparam(param, paramvals(0));
      updateparam(param2, paramvals(1) + eps2);
      um3[] = vR(0, XMh, tgv = TGV);
      um3[] -= R;
      um3[] /= eps2;
      ChangeNumbering(J, um3[], temp3);
      matrix tempPms = [[temp1, temp3]];
      ChangeOperator(JlPM, tempPms, parent = Ja);
      R = vJ(0, XMh, tgv = -10);
      um2[] = vM(0, XMh, tgv = -10);
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      um2[] -= um3[];
      updateparam(param2, paramvals(1));
      um3[] = vJ(0, XMh, tgv = -10);
      Ml1 += um3[];
      um2[] += um3[];
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      Jl1 -= um3[];
      R -= um3[];
      um3[] = vM(0, XMh, tgv = -10);
      Ml1 -= um3[];
      um2[] -= um3[];
      J = vH(XMh, XMh, tgv = 0);
      MatMultTranspose(J, qma, temp3);
      MatMultTranspose(J, qpma, qm);
      qm *= -1.0;
      J = vdM(XMh, XMh, tgv = 0);
      MatMultTranspose(J, qma, temp1);
      qm += temp1;
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      J = vH(XMh, XMh, tgv = 0);
      MatMultTranspose(J, qma, temp1);
      qm -= temp1;
      tempPms = [[temp3, qm]];
      ChangeOperator(gqPM, tempPms, parent = Ja);
      J = vJ(XMh, XMh, tgv = TGV);
      real gl1 = J(uma[], Jl1)/eps;
      real gl2 = J(uma[], R)/eps2;
      ChangeNumbering(J, um3[], qpma, inverse = true);
      real hl1 = (J(uma[], Ml1) - J(um3[], Jl1))/eps;
      real hl2 = (J(uma[], um2[]) - J(um3[], R))/eps2;
      tempPms = [[gl1, gl2], [hl1, hl2]];
      ChangeOperator(glPM, tempPms, parent = Ja);
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
ChangeNumbering(J, ub[], qa, exchange = true, inverse = true);
qa.resize(Ja.n);
if(mpirank == 0) qa(J.n:Ja.n-1) = paramvals;
ChangeNumbering(J, um[], qm);
M = vM(XMh, XMh, tgv = 0);
MatMult(M, qm, qP);
real Mnorm, local = (qm'*qP);
mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
qP /= sqrt(Mnorm);
if (fileext == "hopf" || fileext == "hoho" || fileext == "foho" || fileext == "bota" || fileext == "baut") ChangeNumbering(J, uma[], qma);
else {
  J = vJ(XMh, XMh, tgv = -2);
  KSPSolveTranspose(J, qP, qma);
}
MatMultTranspose(M, qma, pP);
local = (qma'*qP);
mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
pP /= Mnorm;
// solve nonlinear problem with SNES
int ret;
SNESSolve(Ja, funcJa, funcRa, qa, reason = ret,
          sparams = "-snes_linesearch_type " + sneslinesearchtype + " -snes_monitor -snes_converged_reason -options_left no");
if (ret > 0) { // Save solution if solver converged and output file is given
  ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true);
  if(mpirank == 0) paramvals = qa(J.n:Ja.n-1);
  broadcast(processor(0), paramvals);
  updateparam(param, paramvals(0));
  updateparam(param2, paramvals(1));
  M = vM(XMh, XMh, tgv = 0);
  MatMult(M, qm, qP);
  local = (qm'*qP);
  mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
  Mnorm = sqrt(Mnorm);
  qm /= Mnorm;
  qpm /= Mnorm;
  local = (qpma'*qP)/Mnorm;
  mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
  qma /= Mnorm;
  qpma /= Mnorm;
  MatMult(M, qpm, qP);
  local = (qm'*qP);
  mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
  qpm -= Mnorm*qm;
  MatMult(M, qpm, qP);
  local = (qpma'*qP);
  mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
  qpma -= Mnorm*qma;
  if (normalform){
    real[int,int] qDa(paramnames.n, J.n);
    Mat qPM(J.n, mpirank == 0 ? 1 : 0), pPM(J.n, mpirank == 0 ? 1 : 0); // Initialize Mat objects for bordered matrix
    Ja = [[J, pPM], [qPM', 0]]; // make dummy Jacobian
    set(Ja, sparams = "-ksp_type preonly -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_precondition full"
                    + " -prefix_push fieldsplit_1_ -ksp_type preonly -pc_type redundant -redundant_pc_type lu -prefix_pop"
                    + " -prefix_push fieldsplit_0_ " + KSPparams + " -prefix_pop", setup = 1);
    MatMultTranspose(M, qma, pP);
    matrix tempPms = [[pP]]; // dense array to sparse matrix
    ChangeOperator(pPM, tempPms, parent = Ja); // send to Mat
    MatMult(M, qm, qP);
    tempPms = [[qP]]; // dense array to sparse matrix
    ChangeOperator(qPM, tempPms, parent = Ja); // send to Mat
    // 2nd-order
    //  A: base modification due to parameter changes
    J = vJ(XMh, XMh, tgv = TGV);
    ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
    ChangeNumbering(J, uma[], qma, inverse = true);
    if(paramnames[0] != ""){
      for (int k = 0; k < paramnames.n; ++k){
        real paramval = getparam(paramnames[k]);
        updateparam(paramnames[k], paramval + eps);
        um2[] = vR(0, XMh, tgv = TGV);
        updateparam(paramnames[k], paramval);
        um2[] -= R;
        um2[] /= -eps;
        alpha1[paramnames[k]] = -J(uma[], um2[]);
        ChangeNumbering(J, um2[], qP); // FreeFEM to PETSc
        MatMult(M, qpm, pP);
        qP += alpha1[paramnames[k]]*pP;
        qP.resize(Ja.n);
        if(mpirank == 0) qP(Ja.n-1) = 0.0;
        KSPSolve(Ja, qP, pP);
        qDa(k, :) = pP(0:J.n-1);
      }
    }
    // particular solution for quadratic interactions
    um2[] = -0.5*um[];
    um3[] = vH(0, XMh, tgv = -10);
    beta1 = -J(uma[], um3[]);
    ChangeNumbering(J, um3[], qP);
    MatMult(M, qpm, pP);
    qP += beta1*pP;
    qP.resize(Ja.n);
    if(mpirank == 0) qP(Ja.n-1) = 0.0;
    KSPSolve(Ja, qP, pP);
    qP.resize(J.n);
    pP.resize(J.n);
    // inner products
    ChangeNumbering(J, R, qpma, inverse = true);
    beta2 = -2.0*J(R, um3[]);
    R = vdM(0, XMh, tgv = -10);
    ChangeNumbering(J, um2[], qpm, inverse = true, exchange = true);
    um3[] = vH(0, XMh, tgv = -10);
    um3[] += 2.0*R;
    beta2 += J(uma[], um3[]);
    MatMult(M, pP, qP);
    ChangeNumbering(J, R, qP, inverse = true, exchange = true);
    um3[] -= 2.0*R;
    ChangeNumbering(J, R, qpma, inverse = true);
    real gamma = J(R, um3[]);
    R = vdM(0, XMh, tgv = -10);
    gamma -= J(uma[], R);
    um[] = um2[];
    um3[] = vH(0, XMh, tgv = -10);
    beta4 = gamma - 0.5*J(uma[], um3[]);
    ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
    IFMACRO(cubic)
    um2[] = um[];
    um3[] = um[]/6.0;
    R = vT(0, XMh, tgv = -10);
    ENDIFMACRO
    ChangeNumbering(J, um2[], pP, inverse = true, exchange = true);
    IFMACRO(cubic)
    um3[] = vH(0, XMh, tgv = -10);
    R += um3[];
    ENDIFMACRO
    IFMACRO(!cubic)
    R = vH(0, XMh, tgv = -10);
    ENDIFMACRO
    beta3 = beta1*gamma + J(uma[], R);
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
        alpha2[paramnames[k]] = alpha1[paramnames[k]]*gamma + J(uma[], um3[]);
      }
    }
    if(wnlsave){
      complex[int] val(1);
      XMh<complex>[int] defu(vec)(1);
      XMh<complex> defu(um);
      ChangeNumbering(J, um3[], qpm, inverse = true);
      vec[0][].re = um3[];
      val(0) = 0.0;
      savemode(fileout + "_wnl_B", "", fileout + ".bota", meshout, vec, val, sym, true);
      ChangeNumbering(J, um3[], qpma, inverse = true);
      vec[0][].re = um3[];
      savemode(fileout + "_wnl_Badj", "", fileout + ".bota", meshout, vec, val, sym, true);
      if(paramnames[0] != ""){
        for (int k = 0; k < paramnames.n; ++k){
          ChangeNumbering(J, um3[], qDa(k, :), inverse = true);
          vec[0][].re = um3[];
          savemode(fileout + "_wnl_param" + k, "", fileout + ".bota", meshout, vec, val, sym, true);
        }
      }
      ChangeNumbering(J, um3[], pP, inverse = true);
      vec[0][].re = um3[];
      savemode(fileout + "_wnl_AA", "", fileout + ".bota", meshout, vec, val, sym, true);
    }
  }
  else {
    if(paramnames[0] != ""){
      for (int k = 0; k < paramnames.n; ++k){
        alpha1[paramnames[k]] = 0.0;
        alpha2[paramnames[k]] = 0.0;
      }
    }
    beta1 = 0.0;
    beta2 = 0.0;
    beta3 = 0.0;
    beta4 = 0.0;
    ChangeNumbering(J, uma[], qma, inverse = true);
  }
  if(mpirank==0 && adapt) { // Save adapted mesh
    cout << "  Saving adapted mesh '" + meshout + "' in '" + workdir + "'." << endl;
    savemesh(Thg, workdir + meshout);
  }
  ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true);
  ChangeNumbering(J, um[], qm, inverse = true);
  savebota(fileout, "", meshout, alpha1, alpha2, beta1, beta2, beta3, beta4, true, true);
}
```