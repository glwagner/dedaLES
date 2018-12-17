.. math::
    \newcommand{\b}[1]{\boldsymbol{#1}}
    \newcommand{\r}[1]{\mathrm{#1}}
    \newcommand{\bz}{\b{z}}
    \newcommand{\bu}{\b{u}}
    \newcommand{\bcdot}{\b{\cdot}}
    \newcommand{\d}{\partial}

    \newcommand{\p}{\, .}
    \newcommand{\c}{\, ,}


Modified constant Smagorinsky
=============================

The modified constant Smagorinsky closure attempts to improve
the :ref:`constant Smagorinsky` closure in the presence of buoyancy
gradients by multiplying the Smagorinsky subgrid stress by a 
'buoyancy factor' :math:`\lambda` such that 

.. math::

    F^\bu_{ij} = \lambda \, \nu_e S_{ij} \p

The eddy viscosity :math:`\nu_e` and strain tensor :math:`S_{ij}`
are defined as in :ref:`constant Smagorinsky`. The buoyancy factor is

.. math::

    \lambda = \left \{ \begin{matrix}
        1 & \quad \r{for} \quad N^2 \le 0 \\
        \max \left [ 0, \sqrt{ 1 - N^2 / Pr | S|^2 } \right ] & \quad \r{for} \quad N^2 > 0
        \end{matrix} \right . \c

where :math:`Pr` is the turbulent Prandtl number.
The subgrid tracer flux in :ref:`constant Smagorinsky` is not modified.

When a flow is affected by stratification such that
the Richardson-like number :math:`N^2/|S|^2` is greater 
than zero, the buoyancy factor :math:`\lambda` acts to reduce 
the subgrid stress.

References
----------

- `Pressel et al 2015`_

.. _Pressel et al 2015: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015MS000496
