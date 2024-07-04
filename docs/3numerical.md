# Numerical scheme

## Physical-mathematical problem

Kore solves the linearized oscillatory flow inside a near-spherical rotating container for a wide range of physical scenarios. Depending on the scenario, different equations with different terms need to be solved simultanoeusly: the momentum equation, the induction equation and the heat equation. The relevant unknowns or dependent variables in these equations are the flow's velocity $\mathbf u$ and temperature $\theta$, and the magnetic field $\mathbf b$ induced by it. The flow's density can be considered constant in time and space (incompressible flow), or quasi-constant by the Bousinesq approximation, in which case compressibility effects only play a role in some parts of the problem. On the other hand, the flow's pressure does not appear in the final equations, due to the way in which the final equations are written (this is the topic of this section). # Furthermore the flow's velocity and induced magnetic field are both divergenceless (the Bousinesq approximation does not change this), conditions that bring down the number of free components in these vector fields from three to two.

This sections follows the mathematical procedure that transforms the equations from their original form into the form that is implemented in Kore. For the sake of simplificity and without any effect in the process, the equations will be written in their non-dimensional form with $\tau = 1/\Omega$ as the unit of time (see section on non-dimensionalization). Because the problem is oscillatory, all dependent variables can be written as multiples of $e^{\lambda t}$, with $\lambda \in \mathbb C$ being the frequency of oscillation (see e.g. Eq.(3) of Triana et al., 2022). It is the (complex) constants of proportionality that represent the unknowns, while the time derivatives of the original quantities become products of $\lambda$ with these proportionality constants.

On the other hand, both the flow's velocity and induced magnetic field are divergenceless (even in the quasi-incompressible Bousinesq approximation), so they accept a so-called \textit{poloidal-toroidal} decomposition as shown below, where $\mathcal P$ and $\mathcal F$ are the poloidal potentials of $\mathbf u$ and $\mathbf b$, while $\mathcal T$ and $\mathcal G$ are thir toroidal potentials.

$$
\mathcal u = \nabla\times\nabla\times(\mathcal P\mathbf r) + \nabla\times(\mathcal T\mathbf r)\ \ \ \ \ \ \ \ \ \ \ \mathcal b = \nabla\times\nabla\times(\mathcal F\mathbf r) + \nabla\times(\mathcal G\mathbf r)
$$


