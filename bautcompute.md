# bautcompute.md
Author: Chris Douglas ([@cmdoug](https://github.com/cmdoug)) [christopher.douglas@duke.edu](mailto:christopher.douglas@duke.edu)

This script computes the normal form at a Bautin (generalized Hopf) point.

The normal form is written for the complex amplitude $Z = A \exp(\mathrm{i} \omega t)$ as:

$$
\frac{dA}{dt} + \alpha \cdot \delta\lambda A + \beta A |A|^2 = 0
$$

where:
- $\alpha$ is the coefficient for the term from parameter changes,
- $\delta\lambda$ are the parameter increments,
- $\beta_i$ are the coefficients for the term from harmonic interactions.

#### RESIDUAL EVALUATION IN MINIMALLY AUGMENTED FORMULATION
We can directly compute the residual using the varf `vR()`.

To build the augmented residual `Ra`, we augment the residual with two additional functions:

$$
\begin{aligned}
g &= \left\langle{}v,\left(\mathrm{i}\omega{}\mathcal{M} + \mathcal{J}\right)w\right\rangle = v^H\mathcal{L}_{\omega}w \\
h &= \frac{\left\langle{}v,\mathcal{F}\left(w,w,w^\ast\right)\right\rangle}{\left\langle{}v,\mathcal{M}w\right\rangle} = \frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{v^H\mathcal{M}w}
\end{aligned}
$$

where $g$ is the Hopf residual, $h$ is the Bautin residual, and $v$ and $w$ are the adjoint and direct eigenvectors, respectively. Here, $`\mathcal{F}`$ is

$$
\mathcal{F}(a,b,c) = \tfrac{1}{2}\mathcal{T}(a,b,c) - \tfrac{1}{2}\mathcal{H}\left(\mathcal{L}_{2\omega}^{-1}\mathcal{H}(a,b),c\right) - \tfrac{1}{2}\mathcal{H}\left(a,\mathcal{J}^{-1}\mathcal{H}(b,c)\right) - \tfrac{1}{2}\mathcal{H}\left(b,\mathcal{J}^{-1}\mathcal{H}(a,c)\right)
$$

$g$, $v$, and $w$ can be found using minimially augmented systems (For more details, see Govaerts, (2000), Ch. 4, particularly page 87.):

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{L}_\omega & \mathcal{M}q_0 \\
q_0^H\mathcal{M}^H & 0
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

where $q_0$ is an initial approximation of the direct eigenvector.

This implies:

$$
\mathcal{L}_{\omega}w = \mathcal{M}q_0g\qquad{}\text{and}\qquad{}q_0^H\mathcal{M}^Hw = 1
$$

so

$$
w = \mathcal{L}_{\omega}^{-1}\mathcal{M}q_0g\qquad{}\text{and}\qquad{}g = \frac{1}{q_0^H\mathcal{M}^H\mathcal{L}_{\omega}^{-1}\mathcal{M}q_0}.
$$

Note that, at $g = 0$, we have $\mathcal{L}_{\omega}w = 0$ and $q_0^H\mathcal{M}^Hw = 1$.

Similarly, we can find the adjoint eigenmode using the related system:

$$
\begin{bmatrix}
v^H & g
\end{bmatrix}\begin{bmatrix}
-\mathcal{L}_{\omega} & \mathcal{M}q_0 \\
q_0^H\mathcal{M}^H & 0
\end{bmatrix} = \begin{bmatrix}
0 & 1
\end{bmatrix}
$$

This implies:

$$
v^H\mathcal{L}_{\omega} = q_0^H\mathcal{M}^Hg\qquad{}\text{and}\qquad{}v^H\mathcal{M}q_0 = 1
$$

or, taking the complex conjugate transpose:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{L}_{\omega}^H & \mathcal{M}q_0 \\
q_0^H\mathcal{M}^H & 0
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
\mathcal{L}_{\omega}^Hv = \mathcal{M}q_0g^\ast\qquad{}\text{and}\qquad{}q_0^H\mathcal{M}^Hv = 1
$$

so

$$
v = \mathcal{L}_{\omega}^{-H}\mathcal{M}q_0g^\ast\qquad{}\text{and}\qquad{}g^\ast = \frac{1}{q_0^H\mathcal{M}^H\mathcal{L}_{\omega}^{-H}\mathcal{M}q_0}
$$

At $g^\ast = 0$, we have $\mathcal{L}_{\omega}^Hv = 0$ and $q_0^H\mathcal{M}^Hv = 1$, so $v^H\mathcal{L}_{\omega} = 0$ and $v^H\mathcal{M}q_0 = 1$.

Finally, the augmented Bautin residual can be computed directly:

$$
\begin{equation}
h = \frac{\left\langle{}v,\mathcal{F}\left(w,w,w^\ast\right)\right\rangle}{\left\langle{}v,\mathcal{M}w\right\rangle}=\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{v^H\mathcal{M}w}.
\end{equation}
$$

Note that only the real part, $`\Re\lbrace{}h\rbrace{}`$, is evaluated as a residual.

#### JACOBIAN CONSTRUCTION IN MINIMALLY AUGMENTED FORMULATION
Having computed the RHS of the augmented system in `funcRa`, we now have to build the augmented Jacobian matrix for the Newton scheme:

$$
\begin{equation}
\begin{bmatrix}
\mathcal{J} & \frac{\partial\mathcal{R}}{\partial \lambda_1} & 0 & \frac{\partial\mathcal{R}}{\partial \lambda_2} \\
\Re\left(v^H\frac{\partial \mathcal{L}_{\omega}}{\partial q}w\right)^T & \Re\left(v^H\frac{\partial \mathcal{L}_{\omega}}{\partial \lambda_1}w\right) & \Re\left(v^H\frac{\partial \mathcal{L}_{\omega}}{\partial \omega}w\right) & \Re\left(v^H\frac{\partial \mathcal{L}_{\omega}}{\partial \lambda_2}w\right) \\
\Im\left(v^H\frac{\partial \mathcal{L}_{\omega}}{\partial q}w\right)^T & \Im\left(v^H\frac{\partial \mathcal{L}_{\omega}}{\partial \lambda_1}w\right) & \Im\left(v^H\frac{\partial \mathcal{L}_{\omega}}{\partial \omega}w\right) & \Im\left(v^H\frac{\partial \mathcal{L}_{\omega}}{\partial \lambda_2}w\right) \\
\Re\left(\frac{\partial h}{\partial q}\right)^T & \Re\left(\frac{\partial h}{\partial \lambda_1}\right) & \Re\left(\frac{\partial h}{\partial \omega}\right) & \Re\left(\frac{\partial h}{\partial \lambda_2}\right) 
\end{bmatrix}
\begin{bmatrix}
\delta{}q \\
\delta\lambda_1 \\
\delta\omega \\
\delta\lambda_2
\end{bmatrix} = \begin{bmatrix}
\mathcal{R} \\
\Re(g) \\
\Im(g) \\
\Re(h)
\end{bmatrix}
\end{equation}
$$

where $g = v^H\mathcal{L}_{\omega}w$.

To determine the matrix entries in the second and third rows, we differentiate Eq. (1) along each $z$ in $q, \lambda, \omega$ to find:

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{L}_{\omega} & \mathcal{M}q_0 \\
q_0^H\mathcal{M}^H & 0
\end{bmatrix}
\begin{bmatrix}
\frac{\partial w}{\partial z} \\
\frac{\partial g}{\partial z}
\end{bmatrix} = \begin{bmatrix}
\frac{\partial\mathcal{L}_{\omega}}{\partial z}w \\
0
\end{bmatrix}
\end{equation}
$$

We now left-multiply Eq. (4) by $\begin{bmatrix}v^H & g\end{bmatrix}$, finding due to Eq. (2) that:

$$
\frac{\partial g}{\partial z} = v^H\frac{\partial \mathcal{L}_{\omega}}{\partial z}w
$$


This also implies that:

$$
\mathcal{L}_{\omega}\frac{\partial w}{\partial z}=-\frac{\partial\mathcal{L}_{\omega}}{\partial z}w+\mathcal{M}q_0\frac{\partial g}{\partial z}
$$

Similarly differentiating Eq. (2), we find

$$
\begin{equation}
\begin{bmatrix}
-\mathcal{L}_{\omega}^H & \mathcal{M}q_0 \\
q_0^H\mathcal{M}^H & 0
\end{bmatrix}
\begin{bmatrix}
\frac{\partial v}{\partial z} \\
\left(\frac{\partial g}{\partial z}\right)^H
\end{bmatrix} = \begin{bmatrix}
\left(
\frac{\partial\mathcal{L}_{\omega}}{\partial z}\right)^Hv \\
0
\end{bmatrix}
\end{equation}
$$

which similarly yields after left-multiplication by $`\begin{bmatrix}w^H & g^\ast\end{bmatrix}`$ and application of Eq. (1):

$$
\left(\frac{\partial g}{\partial z}\right)^H = w^H\left(\frac{\partial \mathcal{L}_{\omega}}{\partial z}\right)^Hv
$$

We can also find:

$$
\mathcal{L}_{\omega}^H\frac{\partial v}{\partial z}=-\left(\frac{\partial\mathcal{L}_{\omega}}{\partial z}\right)^Hv+\mathcal{M}q_0\left(\frac{\partial g}{\partial z}\right)^H
$$

To determine the augmented matrix entries in the final row, we differentiate Eq. (3) along each $z$ in $q,\lambda_1,\omega,\lambda_2$ to find:

$$
\begin{equation}
\begin{aligned}
\Re\left\lbrace \frac{\partial h}{\partial z}\right\rbrace &= \Re\left\lbrace \left(\frac{\partial v}{\partial z}\right)^H\left[\frac{\mathcal{F}(w,w,w^\ast)}{v^H\mathcal{M}w}-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\mathcal{M}w\right]\right\rbrace \\
&\quad + \Re\left\lbrace \frac{v^H}{v^H\mathcal{M}w}\left[2\mathcal{F}\left(w,\frac{\partial w}{\partial z},w^\ast\right) + \mathcal{F}\left(w,w,\frac{\partial w^\ast}{\partial z}\right)-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{v^H\mathcal{M}w}\mathcal{M}\frac{\partial w}{\partial z}\right]\right\rbrace \\
&\quad + \Re\left\lbrace \frac{v^H}{v^H\mathcal{M}w}\left[\frac{\partial \mathcal{F}}{\partial z}(w,w,w^\ast)-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{v^H\mathcal{M}w}\frac{\partial\mathcal{M}}{\partial z}w\right]\right\rbrace.
\end{aligned}
\end{equation}
$$

However, it is not desirable or necessary to ever construct these dense $`\frac{\partial (\cdot)}{\partial z}`$ operators explicitly. Instead, we focus on their action in the associated inner products.

For the first term, we have:

$$
\begin{aligned}
&\Re\left\lbrace \left(\frac{\partial v}{\partial z}\right)^H\left[\frac{\mathcal{F}(w,w,w^\ast)}{v^H\mathcal{M}w}-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\mathcal{M}w\right]\right\rbrace \\
&\qquad = \Re\left\lbrace\left[-v^H\frac{\partial\mathcal{L}_{\omega}}{\partial z}+\frac{\partial g}{\partial z}q_0^H\mathcal{M}^H\right]\mathcal{L}_{\omega}^{-1}\left[\frac{\mathcal{F}(w,w,w^\ast)}{v^H\mathcal{M}w}-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\mathcal{M}w\right]\right\rbrace \\
&\qquad = \Re\left\lbrace v^H\frac{\partial\mathcal{L}_{\omega}}{\partial z}\hat{w}\right\rbrace
\end{aligned}
$$

where $`\hat{w}`$ solves the non-singular system:

$$
\begin{bmatrix}
-\mathcal{L}_{\omega} & \mathcal{M}q_0 \\
q_0^H\mathcal{M}^H & 0
\end{bmatrix}
\begin{bmatrix}
\hat{w} \\
\mu_1
\end{bmatrix} = \begin{bmatrix}
\frac{\mathcal{F}(w,w,w^\ast)}{v^H\mathcal{M}w}-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\mathcal{M}w \\
0
\end{bmatrix}
$$

giving equivalently,

$$
-\mathcal{L}_{\omega}\hat{w}=\frac{\mathcal{F}(w,w,w^\ast)}{v^H\mathcal{M}w}-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\mathcal{M}w - \mathcal{M}q_0\mu_1\qquad{}\text{and}\qquad{}q_0^H\mathcal{M}^H\hat{w}=0
$$

so, using the identities derived above that $`w=\mathcal{L}_{\omega}^{-1}\mathcal{M}q_0g`$ and $`v=\mathcal{L}_{\omega}^{-H}\mathcal{M}q_0g^\ast`$,

$$
\hat{w}=-\mathcal{L}_{\omega}^{-1}\left[\frac{\mathcal{F}(w,w,w^\ast)}{v^H\mathcal{M}w} - \frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\mathcal{M}w\right] + w\frac{\mu_1}{g}\qquad{}\text{and}\qquad{}\mu_1=0
$$

For the second term in Eq. (7), we have:

$$
\begin{aligned}
&\Re\left\lbrace \frac{v^H}{v^H\mathcal{M}w}\left[2\mathcal{F}\left(w,\frac{\partial w}{\partial z},w^\ast\right) + \mathcal{F}\left(w,w,\frac{\partial w^\ast}{\partial z}\right) - \frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{v^H\mathcal{M}w}\mathcal{M}\frac{\partial w}{\partial z}\right]\right\rbrace \\
&\qquad = \Re\left\lbrace \left[2v^H\frac{\mathcal{F}(w,\cdot,w^\ast)}{v^H\mathcal{M}w} - \frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}v^H\mathcal{M} + v^T\left(\frac{\mathcal{F}(w,w,\cdot)}{v^H\mathcal{M}w}\right)^\ast\right]\frac{\partial w}{\partial z}\right\rbrace \\
&\qquad=\Re\left\lbrace \left[2v^H\frac{\mathcal{F}(w,\cdot,w^\ast)}{v^H\mathcal{M}w} - \frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}v^H\mathcal{M} + v^T\left(\frac{\mathcal{F}(w,w,\cdot)}{v^H\mathcal{M}w}\right)^\ast\right]\mathcal{L}_{\omega}^{-1}\left(-\frac{\partial\mathcal{L}_{\omega}}{\partial z}w+\mathcal{M}q_0\frac{\partial g}{\partial z}\right)\right\rbrace \\
&\qquad=\Re\left\lbrace \hat{v}^H\frac{\partial\mathcal{L}_{\omega}}{\partial z}w\right\rbrace
\end{aligned}
$$

where $`\hat{v}`$ solves the non-singular system:

$$
\begin{bmatrix}
-\mathcal{L}_{\omega}^H & \mathcal{M}q_0 \\
q_0^H\mathcal{M}^H & 0
\end{bmatrix}
\begin{bmatrix}
\hat{v} \\
\mu_2
\end{bmatrix} = \begin{bmatrix}
\frac{2\mathcal{F}(w,\cdot,w^\ast)^Hv}{\left(v^H\mathcal{M}w\right)^\ast} - \left(\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\right)^\ast\mathcal{M}^Hv + \frac{\mathcal{F}(w,w,\cdot)^Tv^\ast}{v^H\mathcal{M}w}\\
0
\end{bmatrix}
$$

giving equivalently,

$$
-\mathcal{L}_{\omega}^H\hat{v}=\frac{2\mathcal{F}(w,\cdot,w^\ast)^Hv}{\left(v^H\mathcal{M}w\right)^\ast} - \left(\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\right)^\ast\mathcal{M}^Hv + \frac{\mathcal{F}(w,w,\cdot)^Tv^\ast}{v^H\mathcal{M}w} - \mathcal{M}q_0\mu_2\qquad{}\text{and}\qquad{}q_0^H\mathcal{M}^H\hat{v}=0
$$

so, using the identities derived above that $`w=\mathcal{L}_{\omega}^{-1}\mathcal{M}q_0g`$ and $`v=\mathcal{L}_{\omega}^{-H}\mathcal{M}q_0g^\ast`$,

$$
\hat{v} = -\mathcal{L}_{\omega}^{-H}\left[\frac{2\mathcal{F}(w,\cdot,w^\ast)^Hv}{\left(v^H\mathcal{M}w\right)^\ast} - \left(\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{(v^H\mathcal{M}w)^2}\right)^\ast\mathcal{M}^Hv + \frac{\mathcal{F}(w,w,\cdot)^Tv^\ast}{v^H\mathcal{M}w}\right] + v\frac{\mu_2}{g^\ast} \qquad{} \text{and} \qquad{} \mu_2 = \frac{\mathcal{F}\left(w,w,w^\ast\right)^Hv}{\left(v^H\mathcal{M}w\right)^\ast} + \frac{\mathcal{F}\left(w,w,w^\ast\right)^Tv^\ast}{v^H\mathcal{M}w}=h^\ast+h
$$

Finally, for the third term, we have:

$$
\begin{aligned}
&\Re\left\lbrace \frac{v^H}{v^H\mathcal{M}w}\left[\frac{\partial \mathcal{F}}{\partial z}(w,w,w^\ast)-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{v^H\mathcal{M}w}\frac{\partial\mathcal{M}}{\partial z}w\right]\right\rbrace \\
&\qquad = \Re\bigg\lbrace\tfrac{1}{2}\frac{v^H}{v^H\mathcal{M}w}\frac{\partial \mathcal{T}}{\partial z}(w,w,w^\ast) - \tfrac{1}{2}\frac{v^H}{v^H\mathcal{M}w}\frac{\partial \mathcal{H}}{\partial z}\left(\mathcal{L}_{2\omega}^{-1}\mathcal{H}(w,w),w^\ast\right) - \tfrac{1}{2}\frac{v^H}{v^H\mathcal{M}w}\mathcal{H}\left(\frac{\partial(\mathcal{L}_{2\omega}^{-1}\mathcal{H})}{\partial z}(w,w),w^\ast\right) \\
&\qquad\qquad - \frac{v^H}{v^H\mathcal{M}w}\frac{\partial\mathcal{H}}{\partial z}\left(w,\mathcal{J}^{-1}\mathcal{H}(w,w^\ast)\right) - \frac{v^H}{v^H\mathcal{M}w}\mathcal{H}\left(w,\frac{\partial(\mathcal{J}^{-1}\mathcal{H})}{\partial z}(w,w^\ast)\right) - \frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{\left(v^H\mathcal{M}w\right)^2}v^H\frac{\partial\mathcal{M}}{\partial z}w\bigg\rbrace
\end{aligned}
$$

The third and fifth terms just above can be simplified using the relation: $\frac{\partial A^{-1}}{\partial z}=-A^{-1}\frac{\partial A}{\partial z}A^{-1}$. This yields:

$$
\begin{aligned}
\frac{\partial(\mathcal{L}_{2\omega}^{-1}\mathcal{H})}{\partial z}(w,w) &= \mathcal{L}_{2\omega}^{-1}\left[\frac{\partial \mathcal{H}}{\partial z}(w,w) - \frac{\partial \mathcal{L}_{2\omega}}{\partial z}\mathcal{L}_{2\omega}^{-1}\mathcal{H}(w,w)\right]\\
\frac{\partial(\mathcal{J}^{-1}\mathcal{H})}{\partial z}(w,w^\ast) &= \mathcal{J}^{-1}\left[\frac{\partial \mathcal{H}}{\partial z}(w,w^\ast) - \frac{\partial \mathcal{J}}{\partial z}\mathcal{J}^{-1}\mathcal{H}(w,w^\ast)\right].
\end{aligned}
$$

We can then define $\hat{x}$ and $\hat{y}$ to satisfy:

$$
-\mathcal{L}_{2\omega}^H\hat{x}=\frac{\mathcal{H}(\cdot, w^\star)^Hv}{\left(v^H\mathcal{M}w\right)^\ast} \qquad\text{and}\qquad -\mathcal{J}^H\hat{y}=\frac{\mathcal{H}(w,\cdot)^Hv}{\left(v^H\mathcal{M}w\right)^\ast}
$$

such that

$$
\begin{aligned}
-\tfrac{1}{2}v^H\mathcal{H}\left(\frac{\partial(\mathcal{L}_{2\omega}^{-1}\mathcal{H})}{\partial z}(w,w),w^\ast\right) &= \tfrac{1}{2}\hat{x}^H\left[\frac{\partial \mathcal{H}}{\partial z}(w,w) - \frac{\partial \mathcal{L}_{2\omega}}{\partial z}\mathcal{L}_{2\omega}^{-1}\mathcal{H}(w,w)\right] \\
-v^H\mathcal{H}\left(w,\frac{\partial(\mathcal{J}^{-1}\mathcal{H})}{\partial z}(w,w^\ast)\right) &= \hat{y}^H\left[\frac{\partial \mathcal{H}}{\partial z}(w,w^\ast) - \frac{\partial \mathcal{J}}{\partial z}\mathcal{J}^{-1}\mathcal{H}(w,w^\ast)\right].
\end{aligned}
$$

Therefore, we have:

$$
\begin{aligned}
&\Re\left\lbrace \frac{v^H}{v^H\mathcal{M}w}\left[\frac{\partial \mathcal{F}}{\partial z}(w,w,w^\ast)-\frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{v^H\mathcal{M}w}\frac{\partial\mathcal{M}}{\partial z}w\right]\right\rbrace \\
&\qquad = \Re\bigg\lbrace \frac{v^H}{v^H\mathcal{M}w}\left[\tfrac{1}{2}\frac{\partial \mathcal{T}}{\partial z}(w,w,w^\ast) - \frac{\partial\mathcal{H}}{\partial z}\left(w,\mathcal{J}^{-1}\mathcal{H}(w,w^\ast)\right) - \tfrac{1}{2}\frac{\partial \mathcal{H}}{\partial z}\left(\mathcal{L}_{2\omega}^{-1}\mathcal{H}(w,w),w^\ast\right)\right] \\
&\qquad\qquad + \tfrac{1}{2}\hat{x}^H\left[\frac{\partial \mathcal{H}}{\partial z}(w,w) - \frac{\partial \mathcal{L}_{2\omega}}{\partial z}\mathcal{L}_{2\omega}^{-1}\mathcal{H}(w,w)\right] + \hat{y}^H\left[\frac{\partial \mathcal{H}}{\partial z}(w,w^\ast) - \frac{\partial \mathcal{J}}{\partial z}\mathcal{J}^{-1}\mathcal{H}(w,w^\ast)\right] - \frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{\left(v^H\mathcal{M}w\right)^2}v^H\frac{\partial\mathcal{M}}{\partial z}w\bigg\rbrace
\end{aligned}
$$

So, putting it all together, we get:

$$
\begin{aligned}
\Re\left\lbrace \frac{\partial h}{\partial z}\right\rbrace &= \Re\left\lbrace \frac{v^H}{v^H\mathcal{M}w}\left[\tfrac{1}{2}\frac{\partial \mathcal{T}}{\partial z}(w,w,w^\ast) + \frac{\partial\mathcal{H}}{\partial z}\left(w,\hat{q}_{AA^\ast}\right) + \frac{\partial \mathcal{H}}{\partial z}\left(\hat{q}_{AA},w^\ast\right) - \frac{v^H\mathcal{F}\left(w,w,w^\ast\right)}{v^H\mathcal{M}w}\frac{\partial\mathcal{M}}{\partial z}w\right]\right\rbrace \\
&\quad + \Re\left\lbrace v^H\frac{\partial\mathcal{L}_{\omega}}{\partial z}\hat{w}\right\rbrace + \Re\left\lbrace\hat{v}^H\frac{\partial\mathcal{L}_{\omega}}{\partial z}w\right\rbrace + \Re\left\lbrace \hat{x}^H\left[\tfrac{1}{2}\frac{\partial \mathcal{H}}{\partial z}(w,w) + \frac{\partial \mathcal{L}_{2\omega}}{\partial z}\hat{q}_{AA}\right]\right\rbrace + \Re\left\lbrace\hat{y}^H\left[\frac{\partial \mathcal{H}}{\partial z}(w,w^\ast) + \frac{\partial \mathcal{J}}{\partial z}\hat{q}_{AA^\ast}\right] \right\rbrace
\end{aligned}
$$

where $`\hat{q}_{AA} = -\tfrac{1}{2}\mathcal{L}_{2\omega}^{-1}\mathcal{H}(w,w)`$ and $`\hat{q}_{AA^\ast} = -\mathcal{J}^{-1}\mathcal{H}(w,w^\ast)`$.

## EXAMPLE USAGE:
### Initialize with Bautin guess from base file, solve on same mesh
```sh
ff-mpirun -np 4 bautcompute.md -param <PARAM> -fi <FILEIN> -bfi <BASEFILEIN> -fo <FILEOUT>
```

### Initialize with Bautin from base and mode file, solve on same mesh
```sh
ff-mpirun -np 4 bautcompute.md -param <PARAM> -fi <FILEIN> -fo <FILEOUT>
```

### Initialize with Bautin guess from file on a mesh from file
```sh
ff-mpirun -np 4 bautcompute.md -param <PARAM> -mi <MESHIN> -bfi <BASEFILEIN> -fi <FILEIN> -fo <FILEOUT>
```

### Initialize with Bautin from file, adapt mesh/solution
```sh
ff-mpirun -np 4 bautcompute.md -param <PARAM> -fi <FILEIN> -fo <FILEOUT> -mo <MESHOUT>
```

NOTE: This file should not be changed unless you know what you're doing.

SEE ALSO: [modecompute.md](./modecompute.md), [hopfcompute.md](./hopfcompute.md), [hopfcontinue.md](./hopfcontinue.md), [fohocompute.md](./fohocompute.md), [./botacompute.md](./botacompute.md), [hohocompute.md](./hohocompute.md), [porbcontinue.md](./porbcontinue.md)

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
bool zerofreq = getARGV("-zero", 0);
string param = getARGV("-param", "");
string param2 = getARGV("-param2", "");
string adaptto = getARGV("-adaptto", "b");
real eps = getARGV("-eps", 1e-7);
real eps2 = getARGV("-eps2", 1e-7);
string sneslinesearchtype = getARGV("-snes_linesearch_type", "none");
real TGV = getARGV("-tgv", -1);
real[int] sym1(sym.n);
real omega;
complex[string] alpha;
complex beta;

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
if (fileext == "baut") {
  ub[].re = loadbaut(fileroot, meshin, um[], uma[], sym1, omega, alpha, beta);
}
else if (fileext == "hopf") {
  ub[].re = loadhopf(fileroot, meshin, um[], uma[], sym1, omega, alpha, beta);
}
else if(fileext == "bota") {
  real[string] alpha1, alpha2;
  real beta1, beta2, beta3, beta4;
  ub[].re = loadbota(fileroot, meshin, um[].re, uma[].re, alpha1, alpha2, beta1, beta2, beta3, beta4);
}
else if (fileext == "foho") {
  real[string] alpha2;
  real beta22, beta23, gamma22, gamma23;
  complex gamma12, gamma13;
  real[int] q2m, q2ma;
  ub[].re = loadfoho(fileroot, meshin, um[], uma[], q2m, q2ma, sym1, omega, alpha, alpha2, beta, beta22, beta23, gamma12, gamma13, gamma22, gamma23);
}
else if(fileext == "hoho") {
  real omegaN;
  complex[string] alphaN;
  complex betaN, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23;
  complex[int] qNm, qNma;
  if(select == 1){
    ub[].re = loadhoho(fileroot, meshin, um[], uma[], qNm, qNma, sym1, sym, omega, omegaN, alpha, alphaN, beta, betaN, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
  }
  else if(select == 2){
    ub[].re = loadhoho(fileroot, meshin, qNm, qNma, um[], uma[], sym, sym1, omegaN, omega, alphaN, alpha, betaN, beta, gamma11, gamma12, gamma13, gamma21, gamma22, gamma23);
  }
}
else if (fileext == "mode") {
  complex eigenvalue;
  um[] = loadmode(fileroot, meshin, sym1, eigenvalue);
  omega = imag(eigenvalue);
}
else if (fileext == "resp") {
  um[] = loadresp(fileroot, meshin, sym1, omega);
}
else if (fileext == "rslv") {
  real gain;
  complex[int] fm;
  um[] = loadrslv(fileroot, meshin, fm, sym1, omega, gain);
}
else if(fileext == "porb") {
  int Nh=1;
  complex[int, int] qh(um[].n, Nh);
  ub[].re = loadporb(fileroot, meshin, qh, sym1, omega, Nh);
  um[] = qh(:, 0);
}
else if(fileext == "floq") {
  int Nh=1;
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
else if(basefileext == "baut") {
  real omega;
  complex[string] alpha;
  complex beta;
  complex[int] qm, qma;
  ub[].re = loadbaut(basefileroot, meshin, qm, qma, sym, omega, alpha, beta);
}
else if(basefileext == "bota") {
  real[string] alpha1, alpha2;
  real beta1, beta2, beta3, beta4;
  real[int] qm, qma;
  ub[].re = loadbota(basefileroot, meshin, qm, qma, alpha1, alpha2, beta1, beta2, beta3, beta4);
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
real[int] paramvals(3-zerofreq);
paramvals(1-zerofreq) = zerofreq ? 0.0 : omega;
paramvals(0) = getparam(param);
paramvals(2-zerofreq) = getparam(param2);
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
complex iomega, iomega2 = 0.0, iomega3 = 0.0;
include "eqns.idp"
Mat<complex> JlPM(J.n, mpirank == 0 ? (3-zerofreq) : 0), gqPM(J.n, mpirank == 0 ? (3-zerofreq) : 0), glPM(mpirank == 0 ? (3-zerofreq) : 0, mpirank == 0 ? (3-zerofreq) : 0); // Initialize Mat objects for bordered matrix
Mat<complex> H(J), Ja = [[J, JlPM], [gqPM', glPM]]; // make dummy Jacobian
complex[int] R(ub[].n), F(J.n), qm(J.n), qma(J.n), qP(J.n), qAA(J.n), qAAs(J.n);
complex ginv, h, s;

func PetscScalar[int] buildJvec(PetscScalar[int] & qm1, real n1){
      ChangeNumbering(J, um[], qm1, inverse = true, exchange = true); // PETSc to FreeFEM
      if (n1 < 0.0) um[] = conj(um[]);
      ik.im = n1*sym1;
      iomega = 1i*n1*omega;
      sym = n1*sym1;
      PetscScalar[int] out = vJ(0, XMh, tgv = -10);
      return out;
}
func int buildJmat(real n1, real tgv){
      ik.im = n1*sym1;
      iomega = 1i*n1*omega;
      sym = n1*sym1;
      J = vJ(XMh, XMh, tgv = tgv);
      return 0;
}
func PetscScalar[int] buildHvec(PetscScalar[int] & qm1, real n1, PetscScalar[int] & qm2, real n2){
      ChangeNumbering(J, um[], qm1, inverse = true, exchange = true); // PETSc to FreeFEM
      if (n1 < 0.0) um[] = conj(um[]);
      ik.im = n1*sym1;
      iomega = 1i*n1*omega;
      ChangeNumbering(J, um2[], qm2, inverse = true, exchange = true); // PETSc to FreeFEM
      if (n2 < 0.0) um2[] = conj(um2[]);
      ik2.im = n2*sym1;
      iomega2 = 1i*n2*omega;
      sym = (n1 + n2)*sym1;
      PetscScalar[int] out = vH(0, XMh, tgv = -10);
      return out;
}
func int buildHmat(PetscScalar[int] & qm1, real n1, real n2){
      ChangeNumbering(J, um[], qm1, inverse = true, exchange = true); // PETSc to FreeFEM
      if (n1 < 0.0) um[] = conj(um[]);
      ik.im = n1*sym1;
      iomega = 1i*n1*omega;
      ik2.im = n2*sym1;
      iomega2 = 1i*n2*omega;
      sym = (n1 + n2)*sym1;
      H = vH(XMh, XMh, tgv = -10);
      return 0;
}
IFMACRO(cubic)
func PetscScalar[int] buildTvec(PetscScalar[int] & qm1){
      ChangeNumbering(J, um[], qm1, inverse = true, exchange = true); // PETSc to FreeFEM
      ik.im = sym1;
      iomega = 1i*omega;
      um2[] = um[];
      ik2.im = sym1;
      iomega2 = 1i*omega;
      um3[] = conj(um[]);
      ik3.im = -sym1;
      iomega3 = -1i*omega;
      sym = sym1;
      PetscScalar[int] out = vT(0, XMh, tgv = -10);
      return out;
}
func int buildTmat(PetscScalar[int] & qm1, real n1, PetscScalar[int] & qm2, real n2, real n3){
      ChangeNumbering(J, um[], qm1, inverse = true, exchange = true); // PETSc to FreeFEM
      if (n1 < 0.0) um[] = conj(um[]);
      ik.im = n1*sym1;
      iomega = 1i*n1*omega;
      ChangeNumbering(J, um2[], qm2, inverse = true, exchange = true); // PETSc to FreeFEM
      if (n2 < 0.0) um2[] = conj(um2[]);
      ik2.im = n2*sym1;
      iomega2 = 1i*n2*omega;
      ik3.im = n3*sym1;
      iomega3 = 1i*n3*omega;
      sym = (n1 + n2 + n3)*sym1;
      H = vT(XMh, XMh, tgv = -10);
      return 0;
}
ENDIFMACRO

// FUNCTIONS
  func PetscScalar[int] funcRa(PetscScalar[int]& qa) {
      ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true);
      if(mpirank == 0) paramvals = qa(J.n:Ja.n-1).re;
      broadcast(processor(0), paramvals);
      updateparam(param, paramvals(0));
      omega = zerofreq ? 0.0 : paramvals(1-zerofreq);
      updateparam(param2, paramvals(2-zerofreq));
      // Hopf part
      buildJmat(1.0, -2);
      KSPSolve(J, qP, qm);
      KSPSolveHermitianTranspose(J, qP, qma);
      h = (qP'*qm);
      mpiAllReduce(h, ginv, mpiCommWorld, mpiSUM);
      qm /= ginv; // rescale direct mode
      qma /= conj(ginv); // rescale adjoint mode
      // Bautin part
      um3[] = -buildHvec(qm, 1.0, qm, -1.0);
    	ChangeNumbering(J, um3[], qAAs);
      buildJmat(0.0, TGV);
    	KSPSolve(J, qAAs, qAAs);
      um3[] = -0.5*buildHvec(qm, 1.0, qm, 1.0);
    	ChangeNumbering(J, um3[], qAA);
      buildJmat(2.0, TGV);
    	KSPSolve(J, qAA, qAA);
      R = buildHvec(qm, 1.0, qAAs, 0.0);
      R += buildHvec(qm, -1.0, qAA, 2.0);
	    IFMACRO(cubic)
      R += 0.5*buildTvec(qm);
      ENDIFMACRO
    	ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ChangeNumbering(J, uma[], qma, inverse = true);
      ik.im = sym1;
      um2[] = vM(0, XMh, tgv = -10);
      s = J(uma[], um2[]);
    	h = J(uma[], R)/s;
      R -= h*um2[];
      R /= -s;
      ChangeNumbering(J, R, F);
      sym = 0;
      R = vR(0, XMh, tgv = TGV);
      PetscScalar[int] Ra;
      ChangeNumbering(J, R, Ra); // FreeFEM to PETSc
      Ra.resize(Ja.n);
      if(mpirank == 0) {
        if (zerofreq) Ra(J.n:Ja.n-1) = [real(1.0/ginv), real(h)];
        else Ra(J.n:Ja.n-1) = [real(1.0/ginv), imag(1.0/ginv), real(h)];
      }
      return Ra;
  }

  func int funcJa(PetscScalar[int]& qa) {
      ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true);
      if(mpirank == 0) paramvals = qa(J.n:Ja.n-1).re;
      broadcast(processor(0), paramvals);
      omega = zerofreq ? 0.0 : paramvals(1-zerofreq);
      broadcast(processor(2-zerofreq), paramvals);
      updateparam(param, paramvals(0) + eps);
      sym = 0;
      um2[] = vR(0, XMh, tgv = TGV);
      um2[] -= R;
      um2[] /= eps;
      PetscScalar[int] temp(J.n), temp2(J.n), temp3(J.n);
      ChangeNumbering(J, um2[], temp);
      updateparam(param, paramvals(0));
      updateparam(param2, paramvals(2-zerofreq) + eps2);
      um2[] = vR(0, XMh, tgv = TGV);
      um2[] -= R;
      um2[] /= eps2;
      ChangeNumbering(J, um2[], temp3);
      updateparam(param2, paramvals(2-zerofreq));
      matrix<PetscScalar> tempPms;
      if(zerofreq) tempPms = [[temp, temp3]];
      else tempPms = [[temp, 0, temp3]]; // dense array to sparse matrix
      ChangeOperator(JlPM, tempPms, parent = Ja); // send to Mat
      PetscScalar[int] what(J.n), vhat(J.n), xhat(J.n), yhat(J.n);
      temp2 = 0.0;
      // Build vhat
      // v^H*T(w, ., w*) + 1/2*conj(v^H*T(w, w, .))
      IFMACRO(cubic)
      buildTmat(qm, 1.0, qm, -1.0, 1.0);
      MatMultHermitianTranspose(H, qma, temp2);
      temp2 /= -conj(s);
      buildTmat(qm, 1.0, qm, 1.0, -1.0);
      MatMultHermitianTranspose(H, qma, temp3);
      temp3 *= 0.5/conj(s);
      temp2 -= conj(temp3);
      ENDIFMACRO
      // build xhat
      buildHmat(qm, -1.0, 2.0);
      MatMultHermitianTranspose(H, qma, temp);
      temp /= -conj(s);
      buildJmat(2.0, TGV);
      KSPSolveHermitianTranspose(J, temp, xhat);
      // v^H*H(L^-1*H(w, .), w*) + 1/2*conj(v^H*H(L^-1*H(w, w), .))
      buildHmat(qm, 1.0, 1.0);
      MatMultHermitianTranspose(H, xhat, temp);
      temp2 -= temp;
      buildHmat(qAA, 2.0, -1.0);
      MatMultHermitianTranspose(H, qma, temp);
      temp /= -conj(s);
      temp2 += conj(temp);
      // build yhat
      buildHmat(qm, 1.0, 0.0);
      MatMultHermitianTranspose(H, qma, temp);
      temp /= -conj(s);
      buildJmat(0.0, TGV);
      KSPSolveHermitianTranspose(J, temp, yhat);
      // v^H*H(w, J^-1*H(., w*)) + conj(v^H*H(w, J^-1*H(w, .)))
      buildHmat(qm, -1.0, 1.0);
      temp3 = yhat;
      temp3 += conj(yhat);
      MatMultHermitianTranspose(H, temp3, temp);
      temp2 -= temp;
      // v^H*H(., J^-1*H(w, w*))
      buildHmat(qAAs, 0.0, 1.0);
      MatMultHermitianTranspose(H, qma, temp);
      temp2 -= temp/conj(s);
      // - conj(h/s)*v^H*M
      ik.im = sym1;
      H = vM(XMh, XMh, tgv = -10);
      MatMultHermitianTranspose(H, qma, temp);
      temp2 += conj(h/s)*temp;

      buildJmat(1.0, -2);
      KSPSolveHermitianTranspose(J, temp2, vhat);
      vhat += (conj(h) + h)*conj(ginv)*qma;

      // Build what
      KSPSolve(J, F, what);

      // dh/dq
      temp3 = 0.0;
      IFMACRO(cubic)
      // + <v, dH/dz(w,qAAs)>
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ChangeNumbering(J, um2[], qAAs, inverse = true, exchange = true);
      ik.im = sym1;
      iomega = 1i*omega;
      ik2.im = 0.0;
      iomega2 = 0.0;
      ik3.im = 0.0;
      iomega3 = 0.0;
      sym = 0;
      H = vT(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, qma, temp3);
      temp3 /= conj(s);
      IFMACRO(quartic)
      // + 1/2*<v, dT/dz(w,w,w*)>
      um2[] = um[];
      um3[] = conj(um[]);
      ik2.im = ik.im;
      iomega2 = iomega;
      ik3.im = -ik.im;
      iomega3 = -iomega;
      H = vQ(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, qma, temp2);
      temp3 += 0.5/conj(s)*temp2;
      ENDIFMACRO
      // + 1/2*<xhat, dHdz(w,w)>
      ik3.im = 0.0;
      iomega3 = 0.0;
      H = vT(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, xhat, temp2);
      temp3 += 0.5*temp2;
      // + <yhat, dHdz(w,w*)>
      ik2.im = -ik.im;
      iomega2 = -iomega;
      H = vT(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, yhat, temp2);
      temp3 += temp2;
      // + <v, dH/dz(qAA,w*)>
      um[] = conj(um[]);
      ik.im = -sym1;
      iomega = -1i*omega;
      ChangeNumbering(J, um2[], qAA, inverse = true, exchange = true);
      ik2.im = 2.0*sym1;
      iomega2 = 2i*omega;
      H = vT(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, qma, temp2);
      temp3 += temp2/conj(s);
      ENDIFMACRO
      // + <v, dLw/dz(what)>
      ChangeNumbering(J, um[], what, inverse = true, exchange = true);
      ik.im = sym1;
      iomega = 1i*omega;
      ik2.im = 0.0;
      iomega2 = 0.0;
      H = vH(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, qma, temp2);
      temp3 += temp2;
      // - h*<v, dM/dz(w)>
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      H = vdM(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, qma, temp2);
      temp3 -= conj(h/s)*temp2;
      // + <xhat, dL2w/dz(qAA)>
      ChangeNumbering(J, um[], qAA, inverse = true, exchange = true);
      ik.im = 2.0*sym1;
      iomega = 2i*omega;
      H = vH(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, xhat, temp2); 
      temp3 += temp2;
      // + <yhat, dJ/dz(qAAs)>
      ChangeNumbering(J, um[], qAAs, inverse = true, exchange = true);
      ik.im = 0.0;
      iomega = 0.0;
      H = vH(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, yhat, temp2); 
      temp3 += temp2;
      // + <vhat, dLw/dz(w)>
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ik.im = sym1;
      iomega = 1i*omega;
      H = vH(XMh, XMh, tgv = 0);
      MatMultHermitianTranspose(H, vhat, temp2);
      temp3 += temp2;
      MatMultHermitianTranspose(H, qma, temp); // compute (dL/dq*w)'*v
      if(!zerofreq) temp2.re = -temp.im;
      temp.im = 0.0;
      temp2.im = 0.0;
      temp3.im = 0.0;
      if(zerofreq) tempPms = [[temp, temp3]];
      else tempPms = [[temp, temp2, temp3]]; // dense array to sparse matrix
      ChangeOperator(gqPM, tempPms, parent = Ja); // send to Mat
      
      // dh/dl1
      updateparam(param, paramvals(0) + eps);
      // + <v, dH/dz(w,qAAs)>
      PetscScalar[int] vterml1 = buildHvec(qm, 1.0, qAAs, 0.0);
      IFMACRO(cubic)
      // + 1/2*<v, dT/dz(w,w,w*)>
      vterml1 += 0.5*buildTvec(qm);
      ENDIFMACRO
      // + <v, dH/dz(qAA,w*)>
      vterml1 += buildHvec(qm, -1.0, qAA, 2.0);
      // - h*<v, dM/dz(w)>
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ik.im = sym1;
      um3[] = vM(0, XMh, tgv = -10);
      vterml1 -= h*um3[];
      vterml1 /= s;
      // + <v, dLw/dz(what)>
      vterml1 += buildJvec(what, 1.0);
      // + 1/2*<xhat, dHdz(w,w)>
      PetscScalar[int] xhterml1 = 0.5*buildHvec(qm, 1.0, qm, 1.0);
      // + <xhat, dL2w/dz(qAA)>
      xhterml1 += buildJvec(qAA, 2.0);
      // + <yhat, dHdz(w,w*)>
      PetscScalar[int] yhterml1 = buildHvec(qm, 1.0, qm, -1.0);
      // + <yhat, dJ/dz(qAAs)>
      yhterml1 += buildJvec(qAAs, 0.0);
      // + <vhat, dLw/dz(w)>
      PetscScalar[int] vhterml1 = buildJvec(qm, 1.0);
      
      updateparam(param, paramvals(0));
	    updateparam(param2, paramvals(2-zerofreq) + eps2);
      // dh/dl2
      // + <v, dH/dz(w,qAAs)>
      PetscScalar[int] vterml2 = buildHvec(qm, 1.0, qAAs, 0.0);
      IFMACRO(cubic)
      // + 1/2*<v, dT/dz(w,w,w*)>
      vterml2 += 0.5*buildTvec(qm);
      ENDIFMACRO
      // <v, dH/dz(qAA,w*)>
      vterml2 += buildHvec(qm, -1.0, qAA, 2.0);
      // - h*<v, dM/dz(w)>
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ik.im = sym1;
      um3[] = vM(0, XMh, tgv = -10);
      vterml2 -= h*um3[];
      vterml2 /= s;
      // + <v, dLw/dz(what)>
      vterml2 += buildJvec(what, 1.0);
      // + 1/2*<xhat, dHdz(w,w)>
      PetscScalar[int] xhterml2 = 0.5*buildHvec(qm, 1.0, qm, 1.0);
      // + <xhat, dL2w/dz(qAA)>
      xhterml2 += buildJvec(qAA, 2.0);
      // + <yhat, dHdz(w,w*)>
      PetscScalar[int] yhterml2 = buildHvec(qm, 1.0, qm, -1.0);
      // + <yhat, dJ/dz(qAAs)>
      yhterml2 += buildJvec(qAAs, 0.0);
      // + <vhat, dLw/dz(w)>
      PetscScalar[int] vhterml2 = buildJvec(qm, 1.0);

	    updateparam(param2, paramvals(2-zerofreq));
      
      // reference point
      // + <v, dH/dz(w,qAAs)>
      R = buildHvec(qm, 1.0, qAAs, 0.0);
      IFMACRO(cubic)
      // + 1/2*<v, dT/dz(w,w,w*)>
      R += 0.5*buildTvec(qm);
      ENDIFMACRO
      // + <v, dH/dz(qAA,w*)>
      R += buildHvec(qm, -1.0, qAA, 2.0);
      // - h*<v, dM/dz(w)>
      ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
      ik.im = sym1;
      um3[] = vM(0, XMh, tgv = -10);
      R -= h*um3[];
      R /= s;
      // + <v, dLw/dz(what)>
      R += buildJvec(what, 1.0);
      vterml1 -= R;
      vterml2 -= R;
      // + 1/2*<xhat, dHdz(w,w)>
      R = 0.5*buildHvec(qm, 1.0, qm, 1.0);
      // + <xhat, dL2w/dz(qAA)>
      R += buildJvec(qAA, 2.0);
      xhterml1 -= R;
      xhterml2 -= R;
      // + <yhat, dHdz(w,w*)>
      R = buildHvec(qm, 1.0, qm, -1.0);
      // + <yhat, dJ/dz(qAAs)>
      R += buildJvec(qAAs, 0.0);
      yhterml1 -= R;
      yhterml2 -= R;
      // + <vhat, dLw/dz(w)>
      R = buildJvec(qm, 1.0);
      vhterml1 -= R;
      vhterml2 -= R;

      PetscScalar gl = J(uma[], vhterml1)/eps;
      PetscScalar gl2 = J(uma[], vhterml2)/eps2;

      PetscScalar hl = J(uma[], vterml1)/eps;
      PetscScalar hl2 = J(uma[], vterml2)/eps2;
      ChangeNumbering(J, um3[], xhat, inverse = true);
      hl += J(um3[], xhterml1)/eps;
      hl2 += J(um3[], xhterml2)/eps2;
      ChangeNumbering(J, um3[], yhat, inverse = true);
      hl += J(um3[], yhterml1)/eps;
      hl2 += J(um3[], yhterml2)/eps2;
      ChangeNumbering(J, um3[], vhat, inverse = true);
      hl += J(um3[], vhterml1)/eps;
      hl2 += J(um3[], vhterml2)/eps2;

      if(zerofreq) tempPms = [[real(gl), real(gl2)],
	  		            				  [real(hl), real(hl2)]];
      else {
        ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
        ik.im = sym1;
        sym = sym1;
        um2[] = vM(0, XMh, tgv = -10);
        PetscScalar gw = J(uma[], um2[]);
        // dh/domega
        // + <vhat, dLw/dz(w)>
        PetscScalar hw = J(um3[], um2[]);

        R = 0.0;
        ChangeNumbering(J, um[], qm, inverse = true, exchange = true); // PETSc to FreeFEM
        ik.im = sym1;
        sym = sym1;
        IFMACRO(cubic)
        // + 1/2*<v, dT/dz(w,w,w*)>
        um2[] = um[];
        ik2.im = sym1;
        um3[] = conj(um[]);
        ik3.im = -sym1;
        R = vddM(0, XMh, tgv = -10);
        R *= 0.5;
        ENDIFMACRO
        // + <v, dH/dz(w,qAAs)>
        ChangeNumbering(J, um2[], qAAs, inverse = true, exchange = true); // PETSc to FreeFEM
        ik2.im = 0.0;
        um3[] = vdM(0, XMh, tgv = -10);
        R += um3[];
        // + <v, dH/dz(qAA,w*)>
        um[] = conj(um[]);
        ik.im = -sym1;
        ChangeNumbering(J, um2[], qAA, inverse = true, exchange = true); // PETSc to FreeFEM
        ik2.im = 2*sym1;
        um3[] = vdM(0, XMh, tgv = -10);
        R += 2.0*um3[];
        um3[] = um[];
        um2[] = um[];
        ik2.im = -sym1;
        um[] = um3[];
        ik.im = 2*sym1;
        um3[] = vdM(0, XMh, tgv = -10);
        R -= um3[];
        R /= s;
        // + <v, dLw/dz(what)>
        ChangeNumbering(J, um[], what, inverse = true, exchange = true); // PETSc to FreeFEM
        ik.im = sym1;
        um3[] = vM(0, XMh, tgv = -10);
        R += um3[];
      
        hw += J(uma[], R);
        // + 1/2*<xhat, dHdz(w,w)>
        ChangeNumbering(J, um[], qm, inverse = true, exchange = true); // PETSc to FreeFEM
        ik.im = sym1;
        um2[] = um[];
        ik2.im = sym1;
        sym = 2*sym1;
        R = vdM(0, XMh, tgv = -10);
        // + <xhat, dL2w/dz(qAA)>
        ChangeNumbering(J, um[], qAA, inverse = true, exchange = true); // PETSc to FreeFEM
        ik.im = 2*sym1;
        um3[] = vM(0, XMh, tgv = -10);
        R += 2.0*um3[];
        
        ChangeNumbering(J, um3[], xhat, inverse = true);
        hw += J(um3[], R);
        tempPms = [[real(gl), -imag(gw), real(gl2)],
				           [imag(gl),  real(gw), imag(gl2)],
				           [real(hl), -imag(hw), real(hl2)]];
      }
      ChangeOperator(glPM, tempPms, parent = Ja); // send to Mat
  	  buildJmat(0.0, TGV);
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
um2[] = vM(0, XMh, tgv = 0);
ChangeNumbering(J, um[], qm);
ChangeNumbering(J, um2[], qP);
complex Mnorm, local = qP.sum;
mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
qm /= Mnorm;
qP /= Mnorm;
local = (qm'*qP);
mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
qP /= sqrt(abs(Mnorm));
// solve nonlinear problem with SNES
int ret;
SNESSolve(Ja, funcJa, funcRa, qa, reason = ret,
          sparams = "-snes_linesearch_type " + sneslinesearchtype + " -snes_monitor -snes_converged_reason -options_left no");
if (ret > 0) { // Save solution if solver converged and output file is given
  ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true, exchange = true);
  if(mpirank == 0) paramvals = qa(J.n:Ja.n-1).re;
  broadcast(processor(0), paramvals);
  updateparam(param, paramvals(0));
  omega = zerofreq ? 0.0 : paramvals(1-zerofreq);
  updateparam(param2, paramvals(2-zerofreq));
  sym = sym1;
  ik.im = sym1;
  J = vM(XMh, XMh, tgv = 0);
  MatMult(J, qm, qP);
  local = qP.sum;
  mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
  qm /= Mnorm;
  qP /= Mnorm;
  local = (qm'*qP);
  mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
  Mnorm = sqrt(abs(Mnorm));
  qP /= Mnorm;
  qm /= Mnorm;
  local = (qP'*qma);
  mpiAllReduce(local, Mnorm, mpiCommWorld, mpiSUM);
  qma /= Mnorm;
  ChangeNumbering(J, uma[], qma, inverse = true);
  if (normalform){
    complex[int,int] qDa(paramnames.n, J.n);
    // 2nd-order
    //  A: base modification due to parameter changes
    ik = 0.0;
    iomega = 0.0;
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
    iomega = 1i*omega;
    iomega2 = -iomega;
    ChangeNumbering(J, um[], qm, inverse = true, exchange = true);
    um2[] = conj(um[]);
    um3[] = vH(0, XMh, tgv = -10);
    um3[].re *= -1.0; // -2.0/2.0
    um3[].im = 0.0;
    ChangeNumbering(J, um3[], qP);
    KSPSolve(J, qP, qP);
    //  C: harmonic generation due to quadratic nonlinear interaction
    ik2.im = sym1;
    iomega2 = iomega;
    um2[] = -0.5*um[];
    sym = 2*sym1;
    um3[] = vH(0, XMh, tgv = -10);
    ChangeNumbering(J, um3[], F);
    ik.im = 2*sym1;
    iomega = 2i*omega;
    J = vJ(XMh, XMh, tgv = TGV);
    KSPSolve(J, F, F);
    // 3rd-order
    //  A: fundamental modification due to parameter change and quadratic interaction of fundamental with 2nd order modification A.
    sym = sym1;
    ik.im = sym1;
    ik2 = 0.0;
    iomega = 1i*omega;
    iomega2 = 0.0;
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
    iomega2 = iomega;
    iomega3 = -iomega;
    um2[] = 0.5*um[];
    um3[] = conj(um[]);
    R = vT(0, XMh, tgv = -10);
    ENDIFMACRO
    //  C: fundamental modification due to quadratic interaction of fundamental with 2nd order modification B
    ik2 = 0.0;
    iomega2 = 0.0;
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
    iomega = -iomega;
    iomega2 = 2i*omega;
    um[] = conj(um[]);
    ChangeNumbering(J, um2[], F, inverse = true, exchange = true);
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
          savemode(fileout + "_wnl_param" + k, "", fileout + ".baut", meshout, vec, val, sym, true);
        }
      }
      ChangeNumbering(J, vec[0][], qP, inverse = true);
      savemode(fileout + "_wnl_AAs", "", fileout + ".baut", meshout, vec, val, sym, true);
      ChangeNumbering(J, vec[0][], F, inverse = true);
      val = 2i*omega;
      sym = 2*sym1;
      savemode(fileout + "_wnl_AA", "", fileout + ".baut", meshout, vec, val, sym, true);
    }
  } else {
    if(paramnames[0] != ""){
      for (int k = 0; k < paramnames.n; ++k){
        alpha[paramnames[k]] = 0.0;
      }
    }
    beta = 0.0;
  }
  if(mpirank==0 && adapt) { // Save adapted mesh
    cout << "  Saving adapted mesh '" + meshout + "' in '" + workdir + "'." << endl;
    savemesh(Thg, workdir + meshout);
  }
  ChangeNumbering(J, ub[], qa(0:J.n-1), inverse = true);
  ChangeNumbering(J, um[], qm, inverse = true);
  savebaut(fileout, "", meshout, sym1, omega, alpha, beta, true, true);
}
```