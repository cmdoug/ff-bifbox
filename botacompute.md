# botacompute.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This script computes the normal form at a non-degenerate Bogdanov-Takens point.

The normal form is written for the amplitude $`A`$ as:

$$
\begin{aligned}
\frac{dA}{dt} - B = 0 \\\\
\frac{dB}{dt} + \alpha_1 \cdot \delta\lambda + \alpha_2 \cdot \delta\lambda A + \beta_1 A^2 + \beta_2 AB = 0
\end{aligned}
$$

or, equivalently,

$$
\ddot{A} + \alpha_1 \cdot \delta\lambda + \alpha_2 \cdot \delta\lambda \dot{A} + \beta_1 A^2 + \beta_2 A\dot{A} = 0
$$


where:
- $`\alpha`$ is the coefficient for the term from parameter changes,
- $`\delta\lambda`$ are the parameter increments,
- $`\beta`$ is the coefficient for the term from harmonic interactions.

#### RESIDUAL EVALUATION IN MINIMALLY AUGMENTED FORMULATION
We can directly compute the residual using the varf `vR()`.

To build the augmented residual `Ra`, we must additionally compute the Hopf residual augmentation:

$$
g = \langle{}v,\mathcal{J}w\rangle = v^H\mathcal{J}w
$$

and the Bogdanov-Takens residual augmentation:

$$
h = \langle{}v,\mathcal{M}w\rangle = v^H\mathcal{M}w
$$

where $`g`$ is the Hopf residual and $`v`$ and $`w`$ are the adjoint and direct eigenvectors, respectively.

$`g`$, $`v`$, and $`w`$ can be found using minimially augmented systems (For more details, see Govaerts, (2000), Ch. 4, particularly page 87.):

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^H & 0
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
\mathcal{J}w = \mathcal{M}p_0g\qquad{}\text{and}\qquad{}(\mathcal{M}q_0)^Hw = 1
$$

so

$$
w = \mathcal{J}^{-1}\mathcal{M}p_0g\qquad{}\text{and}\qquad{}g = \frac{1}{(\mathcal{M}q_0)^H\mathcal{J}^{-1}\mathcal{M}p_0}.
$$

Note that, at $`g = 0`$, we have $`\mathcal{J}w = 0`$ and $`(\mathcal{M}q_0)^Hw = 1`$.

Similarly, we can find the adjoint eigenmode using the related system:

$$
\begin{bmatrix}
v^H & g
\end{bmatrix}\begin{bmatrix}
-\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^H & 0
\end{bmatrix} = \begin{bmatrix}
0 & 1
\end{bmatrix}
$$

This implies:

$$
v^H\mathcal{J} = g(\mathcal{M}q_0)^H\qquad{}\text{and}\qquad{}v^H\mathcal{M}p_0 = 1
$$

or, taking the complex conjugate transpose:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J}^H & \mathcal{M}q_0 \\
(\mathcal{M}p_0)^H & 0
\end{bmatrix}
\begin{bmatrix}
v \\
g^\ast
\end{bmatrix} = \begin{bmatrix}
0 \\
1
\end{bmatrix}
\end{equation}
$$

giving, equivalently,

$$
\mathcal{J}^Hv = \mathcal{M}q_0g^\ast\qquad{}\text{and}\qquad{}(\mathcal{M}p_0)^Hv = 1
$$

so

$$
v = \mathcal{J}^{-H}\mathcal{M}q_0g^\ast\qquad{}\text{and}\qquad{}g^\ast = \frac{1}{(\mathcal{M}p_0)^H\mathcal{J}^{-H}\mathcal{M}q_0}
$$

At $`g^\ast = 0`$, we have $`\mathcal{J}^Hv = 0`$ and $`~`(\mathcal{M}p_0)^Hv = 1`$, so $`v^H\mathcal{J} = 0`$ and $`v^H\mathcal{M}p_0 = 1`$.

Finally, the augmented Bogdanov-Takens residual can be computed directly:

$$
\begin{equation}
h = \langle{}v,\mathcal{M}w\rangle = v^H\mathcal{M}w
\end{equation}
$$

#### JACOBIAN CONSTRUCTION IN MINIMALLY AUGMENTED FORMULATION
Having computed the RHS of the augmented system in `funcRa`, we now have to build the complex augmented Jacobian matrix for the Newton scheme:

$$
\begin{equation}
\begin{bmatrix}
\mathcal{J} & \frac{\partial\mathcal{J}}{\partial \lambda_1} & \frac{\partial\mathcal{J}}{\partial \lambda_2} \\
(\frac{\partial{}g}{\partial q})^H& \frac{\partial{}g}{\partial\lambda_1} & \frac{\partial{}g}{\partial \lambda_2} \\
(\frac{\partial{}h}{\partial q})^H& \frac{\partial{}h}{\partial\lambda_1} & \frac{\partial{}h}{\partial \lambda_2}
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

where $`g = v^H\mathcal{J}w`$ and $`h = v^H\mathcal{M}w`$.

To determine the matrix entries, we differentiate Eq. (1) along each $`z`$ in $`q, \lambda_1, \lambda_2`$ to find:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^H & 0
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

We now left-multiply Eq. (4) by $\begin{bmatrix}v^H & g\end{bmatrix}$, finding due to Eq. (2) that:

$$
\frac{\partial g}{\partial z} = v^H\frac{\partial \mathcal{J}}{\partial z}w
$$

This also implies that:

$$
\mathcal{J}\frac{\partial w}{\partial z}=-\frac{\partial\mathcal{J}}{\partial z}w+\mathcal{M}p_0\frac{\partial g}{\partial z}
$$

Similarly differentiating Eq. (2), we find:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{J}^H & \mathcal{M}q_0 \\
(\mathcal{M}p_0)^H & 0
\end{bmatrix}
\begin{bmatrix}
\frac{\partial v}{\partial z} \\
\left(\frac{\partial g}{\partial z}\right)^H
\end{bmatrix} = \begin{bmatrix}
\left(\frac{\partial\mathcal{J}}{\partial z}\right)^Hv \\
0
\end{bmatrix}
\end{equation}
$$

which similarly yields after left-multiplying by $\left[w^H\quad{}g\right]$ and application of Eq. (1):

$$
\left(\frac{\partial g}{\partial z}\right)^H = w^H\left(\frac{\partial \mathcal{J}}{\partial z}\right)^Hv
$$

This also implies that:

$$
\mathcal{J}^H\frac{\partial v}{\partial z}=-\left(\frac{\partial\mathcal{J}}{\partial z}\right)^Hv+\mathcal{M}q_0\left(\frac{\partial g}{\partial z}\right)^H
$$

To determine the augmented matrix entries in the third row, we differentiate Eq. (3) along each $`z`$ in $`q, \lambda_1, \lambda_2`$ to find:

$$
\begin{equation}
\frac{\partial h}{\partial z} = \left(\frac{\partial v}{\partial z}\right)^H\mathcal{M}w + v^H\frac{\partial \mathcal{M}}{\partial z}w + v^H\mathcal{M}\frac{\partial w}{\partial z}
\end{equation}
$$

However, it is not desirable or necessary to ever construct $`\frac{\partial w}{\partial z}`$ or $`\frac{\partial v}{\partial z}`$ explicitly. Instead of computing these dense operators, we focus on their action in the associated inner products.

For the first term in Eq. (7), we have:

$$
\left(\frac{\partial v}{\partial z}\right)^H\mathcal{M}w=\left(-v^H\frac{\partial\mathcal{J}}{\partial z}+\frac{\partial g}{\partial z}\left(\mathcal{M}q_0\right)^H\right)\mathcal{J}^{-1}\mathcal{M}w=-v^H\frac{\partial\mathcal{J}}{\partial z}\hat{w}+\frac{\partial g}{\partial z}\left(\mathcal{M}q_0\right)^H\hat{w}
$$

where $`\hat{w}`$ solves the non-singular system:

$$
\begin{bmatrix}
\mathcal{J} & \mathcal{M}p_0 \\
(\mathcal{M}q_0)^H & 0
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
\mathcal{J}\hat{w}=\mathcal{M}w-\mathcal{M}p_0h\qquad{}\text{and}\qquad{}(\mathcal{M}q_0)^H\hat{w}=0
$$

so, using the identities derived above that $`w=\mathcal{J}^{-1}\mathcal{M}p_0g`$ and $`v=\mathcal{J}^{-H}\mathcal{M}q_0g^\ast`$,

$$
\hat{w}=\mathcal{J}^{-1}\mathcal{M}w-w\frac{h}{g}\qquad{}\text{and}\qquad{}h=v^H\mathcal{M}w
$$

Then, similarly, for the last term in Eq. (7), we have:

$$
v^H\mathcal{M}\frac{\partial w}{\partial z}=v^H\mathcal{M}\mathcal{J}^{-1}\left(-\frac{\partial\mathcal{J}}{\partial z}w+\mathcal{M}p_0\frac{\partial g}{\partial z}\right)=-\hat{v}^H\frac{\partial\mathcal{J}}{\partial z}w+\hat{v}^H\mathcal{M}p_0\frac{\partial g}{\partial z}
$$

where $`\hat{v}`$ solves the non-singular system:

$$
\begin{bmatrix}
\mathcal{J}^H & \mathcal{M}q_0 \\
(\mathcal{M}p_0)^H & 0
\end{bmatrix}
\begin{bmatrix}
\hat{v} \\
h^\ast
\end{bmatrix} = \begin{bmatrix}
\mathcal{M}^Hv \\
0
\end{bmatrix}
$$

giving equivalently,

$$
\mathcal{J}^H\hat{v}=\mathcal{M}^Hv-\mathcal{M}q_0h^\ast\qquad{}\text{and}\qquad{}(\mathcal{M}p_0)^H\hat{v}=0
$$

so, using the identities derived above that $`w=\mathcal{J}^{-1}\mathcal{M}p_0g`$ and $`v=\mathcal{J}^{-H}\mathcal{M}q_0g^\ast`$,

$$
\hat{v}=\mathcal{J}^{-H}\mathcal{M}^Hv-v\frac{h^\ast}{g^\ast}\qquad{}\text{and}\qquad{}h^\ast=w^H\mathcal{M}^Hv
$$

So we can write Eq. (3) explicitly as

$$
\begin{bmatrix}
\mathcal{J} & \frac{\partial\mathcal{J}}{\partial \lambda_1} & \frac{\partial\mathcal{J}}{\partial \lambda_2} \\
\Re\left(v^H\frac{\partial \mathcal{J}}{\partial q}w\right) & \Re\left(v^H\frac{\partial \mathcal{J}}{\partial \lambda_1}w\right) & \Re\left(v^H\frac{\partial \mathcal{J}}{\partial \lambda_2}w\right) \\
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
\frac{\partial h}{\partial z} = -v^H\frac{\partial\mathcal{J}}{\partial z}\hat{w} + v^H\frac{\partial \mathcal{M}}{\partial z}w - \hat{v}^H\frac{\partial\mathcal{J}}{\partial z}w + \left(\left(\mathcal{M}q_0\right)^H\hat{w}+\hat{v}^H\mathcal{M}p_0\right)\frac{\partial g}{\partial z}
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

SEE ALSO: [modecompute.md](./modecompute.md), [hopfcompute.md](./hopfcompute.md), [hopfcontinue.md](./hopfcontinue.md), [fohocompute.md](./fohocompute.md), [hohocompute.md](./hohocompute.md), [porbcontinue.md](./porbcontinue.md)

```freefem
load "iovtk"
load "PETSc-complex"
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
real[int] sym1(sym.n);
complex[string] alpha1, alpha2;
complex beta1, beta2;

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
XMh<complex> defu(ub), defu(um), defu(uma), defu(um2), defu(um3);
if (fileext == "bota") {
  ub[].re = loadbota(fileroot, meshin, um[], uma[], sym1, alpha1, alpha2, beta1, beta2);
}
else if (fileext == "hopf") {
  real omega;
  ub[].re = loadhopf(fileroot, meshin, um[], uma[], sym1, omega, alpha1, beta1);
}
else if (fileext == "fold") {
  real[int] qm(um[].n), qma(um[].n);
  real[string] alpha; 
  real beta;
  ub[].re = loadfold(fileroot, meshin, qm, qma, alpha, beta);
  um[].re = qm;
  uma[].re = qma;
}
else if (fileext == "foho") {
  real omega;
  real[string] alpha2;
  real beta22, beta23, gamma22, gamma23;
  complex gamma12, gamma13;
  real[int] q2m, q2ma;
  ub[].re = loadfoho(fileroot, meshin, um[], uma[], q2m, q2ma, sym1, omega, alpha1, alpha2, beta1, beta22, beta23, gamma12, gamma13, gamma22, gamma23);
}
else if(fileext == "hoho") {
  real omega, omegaN;
  complex gamma11, gamma12, gamma13, gamma21, gamma22, gamma23;
  complex[int] qNm, qNma;
  if(select == 1){
    ub[].re = loadhoho(fileroot, meshin, um[], uma[], qNm, qNma, sym1, sym, omega, omegaN, alpha1, alpha2, beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
  }
  else if(select == 2){
    ub[].re = loadhoho(fileroot, meshin, qNm, qNma, um[], uma[], sym, sym1, omegaN, omega, alpha2, alpha1, beta2, beta1, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
  }
}
else if (fileext == "mode") {
  complex eigenvalue;
  um[] = loadmode(fileroot, meshin, sym1, eigenvalue);
}
else if (fileext == "resp") {
  real omega;
  um[] = loadresp(fileroot, meshin, sym1, omega);
}
else if (fileext == "rslv") {
  real omega, gain;
  complex[int] fm;
  um[] = loadrslv(fileroot, meshin, fm, sym1, omega, gain);
}
else if(fileext == "porb") {
  int Nh=1;
  real omega;
  complex[int, int] qh(um[].n, Nh);
  ub[].re = loadporb(fileroot, meshin, qh, sym1, omega, Nh);
  um[] = qh(:, 0);
}
else if(fileext == "floq") {
  int Nh=1;
  real omega;
  complex[int, int] qh(um[].n, 2);
  complex eigenvalue;
  real[int] symtemp(sym.n);
  um[] = loadfloq(fileroot, meshin, qh, sym1, eigenvalue, symtemp, omega, Nh);
}
else assert(false); // invalid input filetype
if (basefileext == "base") {
  ub[].re = loadbase(basefileroot, meshin);
}
else if(basefileext == "fold") {
  real[string] alpha;
  real beta;
  real[int] qm, qma;
  ub[].re = loadfold(basefileroot, meshin, qm, qma, alpha, beta);
}
else if(basefileext == "cusp") {
  real[string] alpha, alphaR;
  real beta;
  real[int] qm, qma;
  ub[].re = loadcusp(basefileroot, meshin, qm, qma, alpha, alphaR, beta);
}
else if(basefileext == "hopf") {
  real omega;
  complex[string] alpha;
  complex beta;
  complex[int] qm, qma;
  ub[].re = loadhopf(basefileroot, meshin, qm, qma, sym, omega, alpha, beta);
}
else if(basefileext == "bota") {
  complex[string] alpha1, alpha2;
  complex beta1, beta2;
  complex[int] qm, qma;
  ub[].re = loadbota(basefileroot, meshin, qm, qma, sym, alpha1, alpha2, beta1, beta2);
}
else if(basefileext == "foho") {
  real omega;
  complex[string] alpha1;
  complex beta1, gamma12, gamma13;
  real[string] alpha2;
  real beta22, beta23, gamma22, gamma23;
  complex[int] q1m, q1ma;
  real[int] q2m, q2ma;
  ub[].re = loadfoho(basefileroot, meshin, q1m, q1ma, q2m, q2ma, sym, omega, alpha1, alpha2, beta1, beta22, beta23, gamma12, gamma13, gamma22, gamma23);
}
else if(basefileext == "hoho") {
  real[int] sym2(sym.n);
  real omega1, omega2;
  complex[string] alpha1, alpha2;
  complex beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23;
  complex[int] q1m, q1ma, q2m, q2ma;
  ub[].re = loadhoho(basefileroot, meshin, q1m, q1ma, q2m, q2ma, sym, sym2, omega1, omega2, alpha1, alpha2, beta1, beta2, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
}
else if(basefileext == "tdns") {
  real time;
  ub[].re = loadtdns(basefileroot, meshin, time);
}
else if(basefileext == "porb") {
  int Nh=1;
  real omega;
  complex[int, int] qh(um[].n, Nh);
  ub[].re = loadporb(basefileroot, meshin, qh, sym, omega, Nh);
}
real[int] paramvals(2);
paramvals(0) = getparam(param);
paramvals(1) = getparam(param2);
// Create distributed Mat
Mat<complex> J;
createMatu(Th, J, Pk);
// MESH ADAPTATION
bool adapt = false;
if(meshout == "") meshout = meshin; // if no adaptation
else { // if output meshfile is given, adapt mesh
  adapt = true;
  meshout = meshout + "." + meshext;
  complex[int] q;
  ChangeNumbering(J, ub[], q);
  ChangeNumbering(J, ub[], q, inverse = true);
  ChangeNumbering(J, um[], q);
  ChangeNumbering(J, um[], q, inverse = true);
  ChangeNumbering(J, uma[], q);
  ChangeNumbering(J, uma[], q, inverse = true);
  XMhg defu(uG), defu(umrG), defu(umiG), defu(umarG), defu(umaiG), defu(tempu), defu(uoG);
  tempu[](restu) = ub[].re; // populate local portion of global soln
  mpiAllReduce(tempu[], uG[], mpiCommWorld, mpiSUM);
  tempu[](restu) = um[].re; // populate local portion of global soln
  mpiAllReduce(tempu[], umrG[], mpiCommWorld, mpiSUM);
  tempu[](restu) = um[].im; // populate local portion of global soln
  mpiAllReduce(tempu[], umiG[], mpiCommWorld, mpiSUM);
  tempu[](restu) = uma[].re; // populate local portion of global soln
  mpiAllReduce(tempu[], umarG[], mpiCommWorld, mpiSUM);
  tempu[](restu) = uma[].im; // populate local portion of global soln
  mpiAllReduce(tempu[], umaiG[], mpiCommWorld, mpiSUM);
  if(mpirank == 0) {  // Perform mesh adaptation (serially) on processor 0
    if(adaptto == "bo") {
      defu(uoG) = initu(defu(umarG)'*defu(umarG));
      defu(tempu) = initu(defu(umaiG)'*defu(umaiG));
      tempu[] += uoG[];
      tempu[] = sqrt(tempu[]);
      uoG[] = (umrG[].*umrG[]);
      uoG[] += (umiG[].*umiG[]);
      uoG[] = sqrt(uoG[]);
      uoG[] .*= tempu[];
    }
    IFMACRO(dimension,2)
      if(adaptto == "b") Thg = adaptmesh(Thg, adaptu(uG), adaptmeshoptions);
      else if(adaptto == "bd") Thg = adaptmesh(Thg, adaptu(uG), adaptu(umrG), adaptu(umiG), adaptmeshoptions);
      else if(adaptto == "ba") Thg = adaptmesh(Thg, adaptu(uG), adaptu(umarG), adaptu(umaiG), adaptmeshoptions);
      else if(adaptto == "bda") Thg = adaptmesh(Thg, adaptu(uG), adaptu(umrG), adaptu(umiG), adaptu(umarG), adaptu(umaiG), adaptmeshoptions);
      else if(adaptto == "bo") Thg = adaptmesh(Thg, adaptu(uG), adaptu(uoG), adaptmeshoptions);
    ENDIFMACRO
    IFMACRO(dimension,3)
      //NOTE: 3D mesh adaptation is still under development.
      load "mshmet"
      load "mmg"
      real anisomax = getARGV("-anisomax",1.0);
      real[int] met((bool(anisomax > 1) ? 6 : 1)*Thg.nv);
      if(adaptto == "b") met = mshmet(Thg, adaptu(uG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
      else if(adaptto == "bd") met = mshmet(Thg, adaptu(uG), adaptu(umrG), adaptu(umiG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
      else if(adaptto == "ba") met = mshmet(Thg, adaptu(uG), adaptu(umarG), adaptu(umaiG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
      else if(adaptto == "bda") met = mshmet(Thg, adaptu(uG), adaptu(umrG), adaptu(umiG), adaptu(umarG), adaptu(umaiG), normalization = getARGV("-normalization",1), aniso = bool(anisomax > 1.0),hmin = getARGV("-hmin", 1.0e-6), hmax = getARGV("-hmax", 1.0e+2), err = getARGV("-err", 1.0e-2));
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
  defu(umrG) = defu(umrG);
  defu(umiG) = defu(umiG);
  defu(umarG) = defu(umarG);
  defu(umaiG) = defu(umaiG);
  Th = Thg;
  Mat<complex> Adapt;
  createMatu(Th, Adapt, Pk);
  J = Adapt;
  defu(ub) = initu(0.0);
  defu(um) = initu(0.0);
  defu(uma) = initu(0.0);
  defu(um2) = initu(0.0);
  defu(um3) = initu(0.0);
  restu.resize(ub[].n); // Change size of restriction operator
  restu = restrict(XMh, XMhg, n2o); // Compute new restriction from global mesh to local mesh
  ub[].re = uG[](restu);
  um[].re = umrG[](restu);
  um[].im = umiG[](restu);
  uma[].re = umarG[](restu);
  uma[].im = umaiG[](restu);
}
// Build bordered block matrix from only Mat components
complex[int] ik(sym.n), ik2(sym.n), ik3(sym.n);
complex iomega = 0.0, iomega2 = 0.0, iomega3 = 0.0;
include "eqns.idp"
Mat<complex> JlPM(J.n, mpirank == 0 ? 2 : 0), gqPM(J.n, mpirank == 0 ? 2 : 0), glPM(mpirank == 0 ? 2 : 0, mpirank == 0 ? 2 : 0); // Initialize Mat objects for bordered matrix
Mat<complex> Ja = [[J, JlPM], [gqPM', glPM]]; // make dummy Jacobian
complex[int] R(ub[].n), qm(J.n), qma(J.n), qpm(J.n), qpma(J.n), pP(J.n), qP(J.n);
PetscScalar ginv, h;
// FUNCTIONS
  func PetscScalar[int] funcRa(PetscScalar[int]& qa) {
      ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true); // PETSc to FreeFEM
      if(mpirank == 0) paramvals = qa(J.n:Ja.n-1).re;
      broadcast(processor(0), paramvals);
      updateparam(param, paramvals(0));
      updateparam(param2, paramvals(1));
      sym = 0;
      ik = 0.0;
      R = vR(0, XMh, tgv = TGV);
      PetscScalar[int] Ra;
      ChangeNumbering(J, R, Ra); // FreeFEM to PETSc
      sym = sym1;
      ik.im = sym1;
      J = vJ(XMh, XMh, tgv = -2);
      KSPSolve(J, pP, qm);
      KSPSolveHermitianTranspose(J, qP, qma);
      PetscScalar ginvl = (qP'*qm);
      mpiAllReduce(ginvl, ginv, mpiCommWorld, mpiSUM);
      qm /= ginv; // rescale direct mode
      qma /= conj(ginv); // rescale adjoint mode
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ChangeNumbering(J, uma[], qma, inverse = true);
      um3[] = vM(0, XMh, tgv = -10);
      h = J(uma[], um3[]);
      Ra.resize(Ja.n); // Append 0 to residual vector on proc 0
      if(mpirank == 0) {
        Ra(J.n) = real(1.0/ginv);
        Ra(Ja.n-1) = real(h);
        cout << 1.0/ginv << " " << h << endl;
      }
      return Ra;
  }

  func int funcJa(PetscScalar[int]& qa) {
      ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true); // PETSc to FreeFEM
      if(mpirank == 0) paramvals = qa(J.n:Ja.n-1).re;
      broadcast(processor(0), paramvals);
      complex[int] temp1, temp3;
      ChangeNumbering(J, um3[], qpm);
      KSPSolve(J, qpm, qpm);
      qpm -= h*ginv*qm;
      ChangeNumbering(J, um[], qma, inverse = true, exchange = true);
      um2[] = vM(0, XMh, tgv = -10);
      ChangeNumbering(J, um2[], temp3);
      KSPSolveHermitianTranspose(J, temp3, qpma);
      qpma -= conj(h)*conj(ginv)*qma;
      qpma *= -1.0;
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      updateparam(param, paramvals(0) + eps);
      sym = 0;
      ik = 0.0;
      um3[] = vR(0, XMh, tgv = TGV);
      um3[] -= R;
      um3[] /= eps;
      ChangeNumbering(J, um3[], temp1);
      sym = sym1;
      ik.im = sym1;
      PetscScalar[int] Jl1 = vJ(0, XMh, tgv = -10);
      PetscScalar[int] Ml1 = vM(0, XMh, tgv = -10);
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      Ml1 -= um3[];
      updateparam(param, paramvals(0));
      um3[] = vJ(0, XMh, tgv = -10);
      Ml1 += um3[];
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      Jl1 -= um3[];
      um3[] = vM(0, XMh, tgv = -10);
      Ml1 -= um3[];
      updateparam(param2, paramvals(1) + eps2);
      sym = 0;
      ik = 0.0;
      um3[] = vR(0, XMh, tgv = TGV);
      um3[] -= R;
      um3[] /= eps2;
      ChangeNumbering(J, um3[], temp3);
      matrix<PetscScalar> tempPms = [[temp1, temp3]];
      ChangeOperator(JlPM, tempPms, parent = Ja);
      sym = sym1;
      ik.im = sym1;
      R = vJ(0, XMh, tgv = -10);
      um2[] = vM(0, XMh, tgv = -10);
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      um2[] -= um3[];
      updateparam(param2, paramvals(1));
      um3[] = vJ(0, XMh, tgv = -10);
      um2[] += um3[];
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      um3[] = vJ(0, XMh, tgv = -10);
      R -= um3[];
      um3[] = vM(0, XMh, tgv = -10);
      um2[] -= um3[];
      sym = 0;
      J = vH(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(J, qma, temp3);
      MatMultHermitianTranspose(J, qpma, qm);
      J = vdM(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(J, qma, temp1);
      qm += temp1;
      ChangeNumbering(J, um[], qpm, inverse = true, exchange = true);
      J = vH(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(J, qma, temp1);
      qm -= temp1;
      temp3.im = 0.0;
      qm.im = 0.0;
      tempPms = [[temp3, qm]];
      ChangeOperator(gqPM, tempPms, parent = Ja);
      PetscScalar gl1 = J(uma[], Jl1)/eps;
      PetscScalar gl2 = J(uma[], R)/eps2;
      ChangeNumbering(J, um3[], qpma, inverse = true);
      PetscScalar hl1 = (J(uma[], Ml1) + J(um3[], Jl1))/eps;
      PetscScalar hl2 = (J(uma[], um2[]) + J(um3[], R))/eps2;
      tempPms = [[real(gl1), real(gl2)], [real(hl1), real(hl2)]];
      ChangeOperator(glPM, tempPms, parent = Ja); 
      ik = 0.0;
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
complex[int] qa;
ChangeNumbering(J, ub[], qa);
qa.resize(Ja.n);
if(mpirank == 0) qa(J.n:Ja.n-1).re = paramvals;
sym = sym1;
ik.im = sym1;
um2[] = vM(0, XMh, tgv = -10);
complex phaseref, phaserefl = um2[].sum;
mpiAllReduce(phaserefl, phaseref, mpiCommWorld, mpiSUM);
um[] /= phaseref;
um2[] /= phaseref;
ChangeNumbering(J, um[], qm);
ChangeNumbering(J, um[], qm, inverse = true);
real Mnorm = sqrt(real(J(um[], um2[])));
um2[] /= Mnorm;
ChangeNumbering(J, um2[], qP);
if (fileext == "hopf" || fileext == "hoho" || fileext == "foho" || fileext == "bota" || fileext == "baut") um[] = uma[];
else {
  J = vJ(XMh, XMh, tgv = -2);
  KSPSolveHermitianTranspose(J, qP, qma);
  ChangeNumbering(J, um[], qma, inverse = true, exchange = true);
}
um2[] = vM(0, XMh, tgv = 0);
ChangeNumbering(J, um[], qm, inverse = true);
um2[] *= (Mnorm/J(um[], um2[])); // so that <uma[],M*um[]> = 1
ChangeNumbering(J, um2[], pP);
// solve nonlinear problem with SNES
int ret;
SNESSolve(Ja, funcJa, funcRa, qa, reason = ret,
          sparams = "-snes_linesearch_type " + sneslinesearchtype + " -snes_monitor -snes_converged_reason -options_left no");
if (ret > 0) { // Save solution if solver converged and output file is given
  ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true);
  if(mpirank == 0) paramvals = qa(J.n:Ja.n-1).re;
  broadcast(processor(0), paramvals);
  updateparam(param, paramvals(0));
  updateparam(param2, paramvals(1));
  ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
  sym = sym1;
  ik.im = sym1;
  um2[] = vM(0, XMh, tgv = 0);
  phaserefl = um2[].sum;
  mpiAllReduce(phaserefl, phaseref, mpiCommWorld, mpiSUM);
  ChangeNumbering(J, um[], qm, inverse = true);
  ChangeNumbering(J, uma[], qma, inverse = true);
  um[] /= phaseref;
  um2[] /= phaseref;
  Mnorm = sqrt(real(J(um[], um2[])));
  um[] /= Mnorm; // so that <um[],M*um[]> = 1
  uma[] *= (Mnorm/J(um2[], uma[])); // so that <uma[],M*um[]> = 1
  ChangeNumbering(J, um[], qm);
  ChangeNumbering(J, uma[], qma);
  if (normalform){
    complex[int,int] qDa(paramnames.n, J.n);
/*    // 2nd-order
    //  A: base modification due to parameter changes
    ik = 0.0;
    ik2 = 0.0;
    sym = 0;
    J = vJ(XMh, XMh, tgv = TGV);
    if(paramnames[0] != ""){
      for (int k = 0; k < paramnames.n; ++k){
        real paramval = getparam(paramnames[k]);
        updateparam(paramnames[k], paramval + eps);
        um[] = vR(0, XMh, tgv = TGV);
        updateparam(paramnames[k], paramval);
        um[] -= R;
        um[] /= -eps;
        ChangeNumbering(J, um[], qP);
        KSPSolve(J, qP, qP);
        qDa(k, :) = qP;
      }
    }
    //  B: base modification due to quadratic nonlinear interaction
    ik.im = sym1;
    ik2.im = -sym1;
    ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
    um2[] = conj(um[]);
    um3[] = vH(0, XMh, tgv = -10);
    um3[].re *= -1.0; // -2.0/2.0
    um3[].im = 0.0;
    ChangeNumbering(J, um3[], qP);
    KSPSolve(J, qP, qP);
    //  C: harmonic generation due to quadratic nonlinear interaction
    ik2.im = sym1;
    um2[] = -0.5*um[];
    sym = 2*sym1;
    um3[] = vH(0, XMh, tgv = -10);
    ChangeNumbering(J, um3[], pP);
    ik.im = 2*sym1;
    J = vJ(XMh, XMh, tgv = TGV);
    KSPSolve(J, pP, pP);
    // 3rd-order
    //  A: fundamental modification due to parameter change and quadratic interaction of fundamental with 2nd order modification A.
    sym = sym1;
    ik.im = sym1;
    ik2 = 0.0;
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
        alpha[paramnames[k]] = J(uma[], um3[]);
      }
    }
    IFMACRO(cubic)
    //  B: fundamental modification due to cubic self-interaction of fundamental
    ik2.im = sym1;
    ik3.im = -sym1;
    um2[] = 0.5*um[];
    um3[] = conj(um[]);
    R = vT(0, XMh, tgv = -10);
    ENDIFMACRO
    //  C: fundamental modification due to quadratic interaction of fundamental with 2nd order modification B
    ik2 = 0.0;
    ChangeNumbering(J, um2[], qP, inverse = true, exchange = true);
    IFMACRO(cubic)
    um3[] = vH(0, XMh, tgv = -10);
    R += um3[];
    ENDIFMACRO
    IFMACRO(!cubic)
    R = vH(0, XMh, tgv = -10);
    ENDIFMACRO
    //  D: fundamental modification due to quadratic interaction of fundamental with 2nd order modification C
    ik.im = -sym1;
    ik2.im = 2*sym1;
    um[] = conj(um[]);
    ChangeNumbering(J, um2[], pP, inverse = true, exchange = true);
    um3[] = vH(0, XMh, tgv = -10);
    R += um3[];
    beta = J(uma[], R);
    if(wnlsave){
      complex[int] val(1);
      XMh<complex>[int] defu(vec)(1);
      sym = 0;
      val = 0.0;
      if(paramnames[0] != ""){
        for (int k = 0; k < paramnames.n; ++k){
          ChangeNumbering(J, vec[0][], qDa(k, :), inverse = true);
          savemode(fileout + "_wnl_param" + k, "", fileout + ".hopf", meshout, vec, val, sym, true);
        }
      }
      ChangeNumbering(J, vec[0][], qP, inverse = true);
      savemode(fileout + "_wnl_AAs", "", fileout + ".hopf", meshout, vec, val, sym, true);
      ChangeNumbering(J, vec[0][], pP, inverse = true);
      val = 2i*omega;
      sym = 2*sym1;
      savemode(fileout + "_wnl_AA", "", fileout + ".hopf", meshout, vec, val, sym, true);
    }*/
  } else {
    if(paramnames[0] != ""){
      for (int k = 0; k < paramnames.n; ++k){
        alpha1[paramnames[k]] = 0.0;
        alpha2[paramnames[k]] = 0.0;
      }
    }
    beta1 = 0.0;
    beta2 = 0.0;
  }
  if(mpirank==0 && adapt) { // Save adapted mesh
    cout << "  Saving adapted mesh '" + meshout + "' in '" + workdir + "'." << endl;
    savemesh(Thg, workdir + meshout);
  }
  ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true);
  ChangeNumbering(J, um[], qm, inverse = true);
  savebota(fileout, "", meshout, sym1, alpha1, alpha2, beta1, beta2, true, true);
}
```