.. math::
    \newcommand{\b}[1]{\boldsymbol{#1}}
    \newcommand{\r}[1]{\mathrm{#1}}
    \newcommand{\bz}{\b{z}}
    \newcommand{\bu}{\b{u}}
    \newcommand{\bcdot}{\b{\cdot}}
    \newcommand{\d}{\partial}

    \newcommand{\p}{\, .}
    \newcommand{\c}{\, ,}


Anisotropic minimum dissipation
===============================

The anisotropic minimum dissipation (AMD) model, like :ref:`constant Smagorinsky`,
is an eddy viscosity and eddy diffusivity model for subgrid stress and
subgrid tracer flux. In its simplest implementation, the AMD eddy viscosity is

.. math::

    \nu_e^\dagger = - \frac{ C_k^2 \delta_k^2 u_{i,k} u_{j,k} S_{ij}}{ u_{m, \ell}^2 } \c

where :math:`C_i` and :math:`\delta_i` are the Poincaré constant and
grid spacing in the :math:`i^{\r{th}}` direction.

Unlike the Smagorinsky models, the minimum-dissipation framework
also supplies a minimum variance-destroying eddy diffusivity. For 
a quantity :math:`\theta` the AMD eddy diffusivity is

.. math::

    \kappa_e^\dagger = 
        - \frac{ C_k^2 \delta_k^2 u_{i,k} \theta_{,k} \theta_{,i}}{ \theta_{,\ell}^2 } 

    
In practice, it is necessary to ensure that neither :math:`\nu_e` nor :math:`\kappa_e`
become negative, so that the implemented eddy viscosity is

.. math::
    
    \nu_e = \r{max} \left [ 0, \nu_e^\dagger \right ] \c

and similarly for :math:`\kappa_e`.

Default parameters
------------------

The default Poincaré constant is :math:`C_i = 1/\sqrt{12}` (see `Abkar et al 2016`_).

.. _Abkar et al 2016: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.1.041701
