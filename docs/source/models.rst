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

* Rotating Boussinesq channel flow


Rotating Stratified Boussinesq equations
----------------------------------------

The rotating, stratified Boussinesq equations are

.. math::

    \d_t \bu + \left ( \bu \cdot \nabla \right ) \bu - f \bz \times \bu + \nabla p = b \bz
        + \nu \nabla^2 \bu + \nabla \cdot \b{F}^{\bu} \c \\
     
    \d_t b + \bu \cdot \nabla b + w N^2 =
        + \kappa \nabla^2 b + \nabla \cdot \b{F}^b \p

The subgrid stress :math:`F^{\bu}_{ij}` and subgrid buoyancy flux :math:`\b{F}^b`
are determined by the subgrid turbulence closure. When these terms are zero, 
the simulation is a direct numerical simulation.
