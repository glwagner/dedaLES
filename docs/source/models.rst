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
    \newcommand{\bcdot}{\b{\cdot}}


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

    \d_t \bu + \left ( \bu \bcdot \bnabla \right ) \bu - f \bz \times \bu + \bnabla p = b \bz
        + \nu \bnabla^2 \bu + \bnabla \bcdot \b{F}^{\bu} \c \\
     
    \d_t b + \bu \bcdot \bnabla b + w N^2 =
        \kappa \bnabla^2 b + \bnabla \bcdot \b{F}^b \c

where :math:`\bu = (u, v, w)` is the velocity field, :math:`p` is pressure, 
:math:`b` is buoyancy, :math:`f` is the Coriolis frequency, :math:`\nu` is viscosity,
:math:`\kappa` is diffusivity, :math:`N^2` is the background buoyancy gradient
and squared buoyancy frequency, :math:`F^{\bu}_{ij}` is the subgrid stress and 
:math:`\b{F}^b` is the subgrid buoyancy flux. The subgrid stress and buoyancy flux
are determined by the subgrid turbulence closure. When these terms are zero, 
the simulation is a direct numerical simulation.
