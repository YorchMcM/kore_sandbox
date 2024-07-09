# Numerical scheme

In order to solve the non-dimensional problem described in Physical Definition and Non-dimensionalization, the implementation of Kore is quite the opposite of straightforward, but the methods are traditional. Below, the translation between the final equations in Non-dimensionalization and the software implementation described in Software Implementation will be presented. An attempt has been made at dividing this section and its steps in a way in which the relationship with the different parts of Kore described in Software Implementation are as clear as possible. The momentum equation without the Lorentz force and the rotation time scale is used to illustrate the procedure, which is extendable to all other equations and terms.

In short, Kore writes $\mathbf u$ and $\mathbf b$ in their poloidal and toroidal components. Discretization of the problem is made in terms of frequency, with spherical harmonics serving for angular discretization and CHebyshev polynomials serving for radial discretization. The differential system in the unknowns then turns into a large algebraic linear system in their coefficients, which is solved as an eigenvalue problem if a forcing is not specified (or solved directly otherwise).

## Poloidal-toroidal decomposition and the *u* and *v* sections

Possibly the best place to start could be the way in which the vector equations are *disected*. This decomposition is not motivated by software, but rather it has been a traditional manner in which to study rotating flows (see e.g. Tilgner 1999). It should be noted that, although it hasn't been explicitely mentioned in Physical Definition, a (quai-) incompressible fluid can only produce a divergence-less velocity field $\mathbf u$ by virtue of the continuity equation. Similarly, Gauss' law for a magnetic field imposes a divergence-less induced magnetic field $\mathbf b$. These conditions impose further constraints on the relationship between the components of each vector field, reducing the number of actual free parameters from three per field to just two.

Rather than using the actual components of the (e.g. velocity) field, it has been a traditional approach to express the vector by mean of two scalar quantites, $\mathcal P$ and $\mathcal T$, called the *poloidal* and *toroidal* potentials respectively, as follows:

$$
\mathbf u = \nabla\times\nabla\times(\mathcal P\mathbf r) + \nabla\times(\mathcal T\mathbf r)
$$

The magnetic field $\mathbf b$ is written in a similar way, but the poloidal and toroidal potentials are denoted $\mathcal F$ and $\mathcal G$ respectively. Note that, because $\mathbf u$ is given as the *curls* of vectors, the condition of zero divergence is automatically satisfied. However, because we now have two unknowns rather than three, the original momentum (or induction) equation cannot be used. Instead, the *u* and *v* sections are used. The *u* section is the radial projection of the second curl of the momentum equation ($\hat{\mathbf r}\cdot\nabla\times\nabla\times$), while the *v* section is the radial projection of the first curl of the momentum equation ($\hat{\mathbf r}\cdot\nabla\times$). For simplicity reasons, however, the projection will be performed here by dot-multiplication with $\mathbf r$ rather than $\hat{\mathbf r}$. The same thing is done for the induction equation, where Kore now calls the sections *f* and *g* rather than *u* and *v* respectively. Each section is now one scalar equation, and the two sections together allow to solve for the two scalar unknowns.

Finally, it will be explicitely noted that the poloidal and toroidal potentials are time-independent complex functions, by virtue of $\mathbf u$ (or $\mathbf b$) having these same properties.

## Angular discretization, operators and symmetry

The introduction of the poloidal and toroidal potentials leaves us with a set of *explicitely* scalar equations and unknowns. Each of this unknowns is a complex number with a purely spatial dependance, which is advantegously expressed in terms of the radial $r$, colatitude $\theta$ and longitude $\varphi$ spherical coordinates given the (quasi-)spherical symmetry of the problem at hand. It is thus reasonable to expand the unknowns in spherical harmonics $Y_l^m(\theta,\varphi) = e^{i\varphi}P_l^m(\cos\theta)$, with coefficients e.g. $\mathcal P_{l,m}$ that only depend on $r$:

$$
\mathcal P = \sum_{l=0}^{L}\sum_{m = -l}^lP_{l,m}(r)Y_l^m(\theta,\varphi)
$$

This greatly simplifies the expressions for all the terms in the *u* and *v* sections. These can be built term by term, in what Kore calls *operators*. Consider, for example, the first term in the momentum equation, $\lambda\mathbf u$ (see Non-dimensionalization). When writing the *v* section, this term will turn into $\mathbf r\cdot\nabla\times\mathbf u$. If one does the (very extensive) math, this is simply $-r^2\nabla^2_s\mathcal P$, where $\nabla^2_s$ represents the surface laplacian, i.e. the standard laplacian with all radial derivatives removed. By virtue of the properties of the spherical harmonics, such a quantiy can be written as:

$$
\mathbf r\cdot\nabla\times\mathbf u = -r^2\nabla^2_s\mathcal P = \sum_{l=0}^{L}\sum_{m = -l}^ll(l+1)P_{l,m}(r)Y_l^m(\theta,\varphi)
$$

This is done with all the terms in the equation. Some need more algebra than others to be written in a useful manner. The Coriolis term is particularly cumbersome. Once a series expansion for the e.g. *v* section of each term is obtained, the equation represents an equality between two series that are imposed to match coefficients. These results in a system of 
