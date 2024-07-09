# Non-dimensionalization

In solving the equations described in Physical Definition, Kore uses *dimensionaless* variables. A non-dimensionalization procedure is always possible thanks to the Buckingham-$\pi$ theorem. For this, the following will be used as units of length, mass, time, temperature respectively: $L$, $\tau$, $\rho L^3$, $\theta^*$. On the other hand, several vectors can be written as the product of a typical value of their magnitude times a non-dimensional vector to provide direction. Thus, one can write $\mathbf g \rightarrow g^*\mathbf g$, $\mathbf B_o \rightarrow B^*\mathbf B_o$ and $\mathbf b \rightarrow B^*\mathbf b$, where the vectors on the left hand sides are dimensional and those on the right hand side are non-dimensional, while the starred quantities are their characteristic values and provide the units. Similarly one can also write $\mathbf\Omega = \Omega\hat{\mathbf z}$ (see Physical Definition). Finally, and because it represents an angular frequency, the non-dimensionalization $\lambda \rightarrow \lambda/\tau$ follows.


the background and induced magnetic fields can be written as $B^*\mathbf B_o$ and $B^*\mathbf b$, where $B^*$ is a typical value of the magnitude of the magnetic field while $\mathbf B_o$ and $\mathbf b$ are now non-dimensional vector quantites. Similarly, one can write $\mathbf\Omega = \Omega\hat{\mathbf z}$ (see Physical Definition).

After dividing the momentum equation by $\rho$, applying all appropriate non-dimensionalizations and cancelling some common terms, the non-dimensional equations read:

$$ \displaylines{
\lambda\mathbf u + 2\Omega\tau\hat{\mathbf z}\times\mathbf u = -\nabla p-\frac{\theta^* g\^ast\alpha\tau^2}{L}\theta\mathbf g + \frac{\nu\tau}{L^2}\nabla^2\mathbf u + \frac{1}{\rho\mu_o}\bigg(\frac{B^*\tau}{L}\bigg)^2(\nabla\times\mathbf b)\times\mathbf B_o \\
\lambda\mathbf b = \nabla\times(\mathbf u\times\mathbf B_o) + \frac{\eta\tau}{L^2}\nabla^2\mathbf b \\
\lambda\theta = -\mathbf u\cdot\nabla T_o + \frac{\kappa\tau}{L^2}\nabla^2\theta
} $$

Here, all vector quantites, derivatives, dependent variables as well as the background temperature profile $T_o$ and the eigenfrequency $\lambda$ are non-dimensional. From here, different standard non-dimensionalizations can take place depending on the choice of time scale. Each of them has their particular non-dimensional numbers of interest.

## Rotation time scale

This non-dimensionalization is characterized by the time scale

$$
\tau = \frac{1}{\Omega}
$$

The relevant non-dimensional numbers are the:

$$ \displaylines{
\text{Ekmann number:}\ \ E = \frac{\nu}{\Omega L^2} \\
\text{Lenhert number:}\ \ Le = \frac{B^*}{\Omega L\sqrt{\rho\mu_o}} \\
\text{Rayleigh number:}\ \ Ra = \alpha g^*\theta^* L^3 \\
\text{Prandtl number:}\ \ Pr = \frac{\nu}{\kappa}
} $$
