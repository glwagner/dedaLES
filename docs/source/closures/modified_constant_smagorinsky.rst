.. math::
    \newcommand{\b}[1]{\boldsymbol{#1}}
    \newcommand{\r}[1]{\mathrm{#1}}
    \newcommand{\bz}{\b{z}}
    \newcommand{\bu}{\b{u}}
    \newcommand{\bcdot}{\b{\cdot}}
    \newcommand{\d}{\partial}

    \newcommand{\p}{\, .}
    \newcommand{\c}{\, ,}


Modified Constant Smagorinsky
=============================

In modified constant Smagorinsky, the subgrid stress
defined in :ref:`constant Smagorinsky` is modified by
a 'buoyancy factor' :math:`\lambda` such that 

.. math::

    F^\bu_{ij} = \lambda \, \nu_e S_{ij} \p

The eddy viscosity :math:`\nu_e` and 
eddy diffusivity are defined as in :ref:`constant Smagorinsky`.
The buoyancy factor, which is the only difference between 
ordinary constant Smagorinsky and 'modified' constant Smagorinsky,
is

.. math::

    \lambda = \left \{ \begin{matrix}
        1 & \quad \r{for} \quad N^2 \le 0 \\
        \max \left [ 0, \sqrt{ 1 - N^2 / Pr | S|^2 } \right ] & \quad \r{for} \quad N^2 > 0
        \end{matrix} \right . \c

where :math:`Pr` is the turbulent Prandtl number.
When the flow is affected by stratification such that
the Richardson-like number :math:`N^2/|S|^2` is greater 
than zero the buoyancy factor reduces the effective subgrid stress.

References
----------

- `Pressel et al 2015`_

.. _Pressel et al 2015: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015MS000496
