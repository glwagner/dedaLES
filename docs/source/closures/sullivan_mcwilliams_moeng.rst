.. math::
    \newcommand{\b}[1]{\boldsymbol{#1}}
    \newcommand{\r}[1]{\mathrm{#1}}
    \newcommand{\bz}{\b{z}}
    \newcommand{\bu}{\b{u}}
    \newcommand{\bcdot}{\b{\cdot}}
    \newcommand{\d}{\partial}

    \newcommand{\p}{\, .}
    \newcommand{\c}{\, ,}

    \newcommand{\bnabla}{\b{\nabla}}

.. _Sullivan-McWilliams-Moeng:

Sullivan-McWilliams-Moeng
=========================

In the second-order 'Sullivan-McWilliams-Moeng' (SMM) turbulence closure, 
a prognostic equation for the subfilter energy :math:`e` is solved in 
addition to equations for resolved velocity and tracer fields.
The subfilter energy obeys

.. math::

    e_t + \bu \bcdot \bnabla e = P + B + \epsilon - D \c

where :math:`P`, :math:`B`, :math:`\epsilon`, and :math:`D` are the
production of subfilter energy by resolved shear, buoyancy conversion 
from subfilter potential energy, dissipation of subfilter energy, and 
turbulent diffusion of filter energy by subfilter turbulent velocity.

The subfilter stress :math:`F^\bu_{ij}` is defined in terms of the 
resolved rate of strain tensor
:math:`S_{ij} = \tfrac{1}{2} \left ( \d_i u_j + \d_j u_i \right )` ,
via

.. math::

    F^\bu_{ij} = 2 \gamma \nu_i S_{ij} + 2 \nu_w \bar{S}_{ij} \c

where :math:`\bar{S}` is the horizonal mean of the strain tensor,
:math:`\nu_i` is the 'interior' eddy viscosity, 
:math:`\nu_w`: is the 'near-wall' eddy viscosity, and :math:`\gamma`
is the 'isoptropy factor':

.. math::

    \gamma = \frac{S'}{S' + \bar S} \c

and :math:`S'` is defined in terms of the strain
perturbation tensor :math:`S'_{ij} = S_{ij} - \bar{S}_{ij}` via

.. math::

    S' = \sqrt{ 2 \overline{S'_{ij} S'_{ij}}} \p
