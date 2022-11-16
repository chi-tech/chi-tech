/**\page CodeTut3 Coding Tutorial 3 - Poisson's problem with Finite Volume

## Table of contents
- \ref CodeTut3Sec1
 - \ref CodeTut3Sec1_1
 - \ref CodeTut3Sec1_2
 - \ref CodeTut3Sec1_3

\section CodeTut3Sec1 1 Poisson's equation
<a href="https://en.wikipedia.org/wiki/Poisson%27s_equation">The Poisson's equation</a> states the following, for
\f$ \phi, q \in \mathcal{R} \f$,
\f{eqnarray*}{
-\boldsymbol{\nabla} \boldsymbol{\cdot} \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr) = q(\mathbf{x}), \quad \quad \quad (\mathbf{x})\in\mathcal{D} \\
\phi(\mathbf{x}) = 0, \quad \quad \quad \mathbf{x}\in \partial \mathcal{D}
\f}
where \f$ \boldsymbol{\nabla} \boldsymbol{\cdot} \bigr( \bigr) \f$ denotes the divergence-operator
and \f$ \boldsymbol{\nabla} \f$ denotes the gradient-operator, i.e.
\f$ \boldsymbol{\nabla} \phi \f$ denotes the gradient of \f$ \phi \f$. The
boundary conditions here state that \f$ \phi=0 \f$ on the boundary.




\subsection CodeTut3Sec1_1 1.1 Our specific problem
For our specific problem we will choose \f$ q(\mathbf{x})=1 \f$ and \f$ \mathcal{D} \f$ a cartesian domain,
either 1D, 2D or 3D, with each dimension always between \f$ -1,+1 \f$. We can generate the mesh for this
problem using an input file
\code
--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=20
L=2
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

--chiMeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,mesh)
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
--chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)
\endcode
This code can be used to generate any of the following meshes,
\image html CodingTutorials/OrthoMesh_1D_2D_3D.png "[From left to right] 1D, 2D, 3D orthogonal mesh" width=1200px




\subsection CodeTut3Sec1_2 1.2 The Weak-form of Poisson's equation
When using the finite element method we develop the so-called weak-form by
weighting the Poisson equation with a trial space function, \f$ t_i(\mathbf{x}) \f$, where \f$ i \f$
is a unique node in the problem, and then integrate over the volume of the
domain
\f[
-\int_V t_i \boldsymbol{\nabla} \boldsymbol{\cdot} \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr) dV
 =
\int_V t_i q(\mathbf{x}) dV.
\f]

Next we apply integration by parts to the left term, for which,
\f[
t_i \boldsymbol{\nabla} \boldsymbol{\cdot} \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr)
 =
\boldsymbol{\nabla} \cdot (t_i \boldsymbol{\nabla} \phi)
-
(\boldsymbol{\nabla} t_i ) \boldsymbol{\cdot} (\boldsymbol{\nabla} \phi)
\f]
which gives
\f[
\int_V (\boldsymbol{\nabla} t_i ) \boldsymbol{\cdot} (\boldsymbol{\nabla} \phi) dV
-\int_V \boldsymbol{\nabla} \boldsymbol{\cdot} (t_i \boldsymbol{\nabla} \phi) dV
 =
\int_V t_i q(\mathbf{x}) dV.
\f]

We then apply Gauss's Divergence Theorem to the second term from the left, to get
\f[
\int_V (\boldsymbol{\nabla} t_i ) \boldsymbol{\cdot} (\boldsymbol{\nabla} \phi) dV
-\int_S \mathbf{n} \boldsymbol{\cdot} (t_i \boldsymbol{\nabla} \phi) dA
 =
\int_V t_i q(\mathbf{x}) dV.
\f]

This is essentially the weak-form.


\subsection CodeTut3Sec1_3 1.3 The discretized equation with the Continuous Galerkin approach
To fully discretize this equation we next approximate the continuous nature of
\f$ \phi(\mathbf{x}) \f$ using essentially a very fancy interpolation scheme.
With this scheme we approximate \f$ \phi(\mathbf{x}) \f$ as
\f$ \phi_h(\mathbf{x}) \f$ where
\f[
\phi(\mathbf{x}) \approx \phi_h(\mathbf{x}) = \sum_j \phi_j b_j(\mathbf{x})
\f]
where the coefficients \f$ \phi_j \f$ do not dependent on space and
\f$ b_j(\mathbf{x}) \f$ are called shape functions. For this tutorial we
will be utilizing the Piecewise Linear Continuous (PWLC) shape functions which
are generally specially in the fact that they can support arbitrary polygons and
polyhedra.

We will also follow a Galerking approach, resulting in the trial-functions
\f$ t_i \f$ being the same as the shape functions, i.e., \f$ t_n(\mathbf{x}) = b_n(\mathbf{x}) \f$.

With this approximation defined, our weak form becomes
\f{align*}{
\int_V (\boldsymbol{\nabla} b_i ) \boldsymbol{\cdot} (\boldsymbol{\nabla} \phi) dV
-\int_S \mathbf{n} \boldsymbol{\cdot} (b_i \boldsymbol{\nabla} \phi) dA
 &=
\int_V b_i q(\mathbf{x}) dV
\\
\sum_j \phi_j \int_V (\boldsymbol{\nabla} b_i ) \boldsymbol{\cdot} (\boldsymbol{\nabla} b_j) dV
-\sum_j \phi_j \int_S \mathbf{n} \boldsymbol{\cdot} (b_i \boldsymbol{\nabla} b_j) dA
 &=
\int_V b_i q(\mathbf{x}) dV
\f}

Now we are faced with the dilemma of how to compute the integral terms in this
equation. The first thing we do here is to observe that the shape functions
\f$ b_n \f$ are defined cell-by-cell, therefore we can drill a little deeper
by defining
\f{align*}{
\int_V f(\mathbf{x}) dV &= \sum_c \int_{V_c} f_c(\mathbf{x}) dV \\
\int_S g(\mathbf{x}) dA &= \sum_c \sum_{\substack{f\\ f\in \partial \mathcal{D}}}
 \int_{S_{c,f}} g_c(\mathbf{x}) dA \\
\f}
where \f$ V_c \f$ is the volume of cell \f$ c \f$, \f$ S_{c,f} \f$ is the surface
of cell \f$ c \f$ face \f$ f \f$. The functions \f$ f_c \f$ and \f$ g_c \f$ are
function uniquely defined within cell \f$ c \f$'s volume and surface respectively.

These integrals are often difficult to perform in real-world coordinates and are
therefore carried out on a reference element in a different coordinate system. This
necessitates the use of a transformation using the magnitude of the Jacobian that
 transforms the above integrals into the form
\f{align*}{
\int_{V_c} f_c(\mathbf{x}) dV
&=
\int_{V_e} f_e(\tilde{\mathbf{x}}) |J(\tilde{\mathbf{x}})| dV
\\
\int_{S_{c,f}} g_c(\mathbf{x}) dA
&=
\int_{S_{c,f}} g_e(\tilde{\mathbf{x}}) |J(\tilde{\mathbf{x}})| dA
\f}
Now, the big advantage of doing this is that we can use quadrature rules for
these integrals
\f{align*}{
\int_{V_e} f_e(\tilde{\mathbf{x}}) |J(\tilde{\mathbf{x}})| dV
&=
\sum_n w_n f_e(\tilde{\mathbf{x}}_n) |J(\tilde{\mathbf{x}}_n)|
\\
\int_{S_{c,f}} g_e(\tilde{\mathbf{x}}) |J(\tilde{\mathbf{x}})| dA
&=
\sum_n w_n g_e(\tilde{\mathbf{x}}_n) |J(\tilde{\mathbf{x}}_n)|
\f}



*/