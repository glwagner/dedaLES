.. math::

    \newcommand{\b}[1]{\boldsymbol{#1}}
    \newcommand{\r}[1]{\mathrm{#1}}
    \newcommand{\bz}{\b{z}}
    \newcommand{\bu}{\b{u}}
    \newcommand{\bcdot}{\b{\cdot}}
    \newcommand{\d}{\partial}

    \newcommand{\p}{\, .}
    \newcommand{\c}{\, ,}


.. _Fluid models:

Fluid models
============

dedaLES provides solvers for

* Boussinesq flow in a channel

.. _Boussinesq channel flow:

Boussinesq Channel flow
-----------------------

The rotating, stratified Boussinesq equations are

.. math::

    \d_t \bu + \left ( \bu \cdot \nabla \right ) \bu - f \bz \times \bu + \nabla p = b \bz
        + \nu \nabla^2 \bu + \nabla \cdot \b{F}^{\bu} \c \\
     
    \d_t b + \bu \cdot \nabla b + w N^2 =
        + \kappa \nabla^2 b + \nabla \cdot \b{F}^b \c

where :math:`\bu = (u, v, w)` is the velocity field, :math:`p` is pressure, 
:math:`b` is buoyancy, :math:`f` is the Coriolis frequency, :math:`\nu` is viscosity,
:math:`\kappa` is diffusivity, :math:`N^2` is the background buoyancy gradient
and squared buoyancy frequency, :math:`F^{\bu}_{ij}` is the subgrid stress and 
:math:`\b{F}^b` is the subgrid buoyancy flux. The subgrid stress and buoyancy flux
are determined by the subgrid turbulence closure. When these terms are zero, 
the simulation is a direct numerical simulation.
