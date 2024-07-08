# Physical definition

Kore can solve the flow inside near-spherical rotating bodies in a wide range of scenarios in the context of planetary liquid cores. The basic physical problem considers the fluid body rotating to a large extent as a rigid body with angular velocity $\mathbf\Omega$, and small velocity perturbations $\vec u$ exist over the background field $\mathbf\Omega\times\mathbf r$. Kore assumes this perturbation to be oscillatory and very small compared to the maximum value of $\mathbf\Omega\times\mathbf r$ attained in the volume's domain. In the most basic of scenarios, Kore solves the linearized and incompressible Navier-Stokes momentum equation and finds $\mathbf u$. In the absence of external forcings, the equation can be solved as an eigenvalue problem, for which the frequencies $\omega$ of oscillation of $\math u$ are obtained together with their respective eigenmodes - characteristic velocity profiles, see Numerical Scheme for more information.

Here, the basic equations defining each problem will be provided, together with their respective boundary conditions. All are given in their homogeneous version - i.e. for the resolution of the egienvalue problem. Forcings will be covered in Numerical Scheme.

## Momentum equation

The most basic of problems poses the Navier-Stokes momentum equation. Under the assumption of oscillatory velocity perturbations, one can write $\vec u \sim e^{\lambda t}$, where $\lambda\in\mathbb C$ encompasses both the damping and oscillatory characteristics of the motion. Time derivatives of $\vec u$ then turn into products by $\lambda$. With this notation, the momentum equation in a reference frame rotating with the container/planet at angular speed $\vec\Omega$ - the \textit{mantle frame} - is written as:

$$
\rho\partial_t \mathbf u +2\rho\mathbf\Omega\times\mathbf u = -\nabla p + \rho\vec g + \rho\nu\nabla^2\mathbf u
$$

Here, $\rho$ and $\nu > 0$ are the fluid's density and kinematic viscosity respectively, while $p$ is the reduce pressure - sum of physical pressure and centrifugal potential - and $\vec g$ is the acceleration of gravity.

In the general case, the momentum equation is imposed in a (near-)spherical shell of inner radius $r_{icb}$ and outer radius $r_{cmb}$ - the \textit{inner core boundary} and the \textit{core mantle boundary}. The boundary condition in this scenario then requires the no-penetration and no-slip boundary condition at both surfaces. Under the assumption that the inner core and mantle rotate at the same angular speef $\vec\Omega$, this boundary condition is trivially written as $\vec u_{icb} = \vec u_{cmb} = \vec 0$. At the moment, Kore cannot handle differential rotation between the inner core and mantle. Alternatively, Kore allows to set stress-free boundary conditions, so that the radial derivative f $\vec u$ vanishes at the boundaries. The BCs at the ICB and the CMB need not be the same.

If the physical model lacks an inner core ($r_{icb} = 0$) the inviscid momentum equation can be considered, i.e. $\nu = 0$. This only leaves one boundary, naemly the CBM, where the only feasible boundary condition is now the stress-free BC. The second boundary condition required to complete the problem is that of regularity at the origin.

## Induction equation

On Earth, the liquid core is conductive iron and its motion produces a magnetic field. In the presence of an imposed background magnetic field $\vec B_o$, the small oscilations $\vec u$ will induce a similarly small and oscillatiory magnetic field perturbation $\vec b \sim e^{\lambda t}$, with magnitude much smaller than $\vec B_o$. The evlution of this magnetic perturbation is given by the electromagnetic induction equation as (see Triana et al. 2021a):

$$
\lambda\vec b = \nabla\times(\vec u\times\vec B_o) + \eta\nabla^2\vec b
$$

Here, $eta$ is the fluid's magnetic diffusivity.

The magnetic boundary conditions are many. TODO: Explicarlas.

In the presence of magnetic fields, the conductive iron will be subject to an electromagnetic force $\vec F_{em}$ (the Lorentz force), which should be added to the right hand side of the momentum equaion. This force is written as:

$$
\vec F_{em} = \frac{1}{\mu_o}(\nabla\times\vec b)\times\vec B_o
$$

Here, $\mu_o$ is the magnetic permeability. The momentum equation - or its boundary conditions - are not affected any further.

## Heat equation
