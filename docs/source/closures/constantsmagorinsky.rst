.. math::
    \newcommand{\b}[1]{\boldsymbol{#1}}
    \newcommand{\r}[1]{\mathrm{#1}}
    \newcommand{\bz}{\b{z}}
    \newcommand{\bu}{\b{u}}
    \newcommand{\bcdot}{\b{\cdot}}
    \newcommand{\d}{\partial}

    \newcommand{\p}{\, .}
    \newcommand{\c}{\, ,}

.. _constant Smagorinsky:

Constant Smagorinsky
====================

In the first-order 'constant Smagorinsky' turbulence closure, the subgrid stress
:math:`F^\bu_{ij}` defined in terms of the resolved rate of strain tensor
:math:`S_{ij} = \tfrac{1}{2} \left ( \d_i u_j + \d_j u_i \right )` 
(which we abbreviate as the 'strain tensor', and an eddy viscosity :math:`\nu_e`:

.. math::

    F^\bu_{ij} = 2 \nu_e S_{ij} \, ,

In dedaLES, the eddy viscosity :math:`\nu_e` is defined via a slight
generalization of traditional constant Smagorinsky,

.. math::

    \nu_e = \left [ \delta_c^2 + \left ( C_s \delta_{\r{c}} \right )^2 \right ] | S | \, ,

where :math:`\delta_c` is a constant 'filter width', 
:math:`C_s` is the 'Smagorinsky coefficient', 
and :math:`\delta` is a filter width defined by
some multiple of the grid resolution, and thus dependent on position 
within the chosen grid in general.
The invariant of the resolved strain tensor :math:`|S|` is

.. math::

    | S | \equiv \sqrt{ 2 S_{ij} S_{ji} } \, .

Note that :math:`S_{ij}` is symmetric, so that :math:`S_{ij} = S_{ji}`.
The subgrid buoyancy flux is

.. math::

    \b{F}^b = \kappa_e \nabla b \c

with :math:`\kappa_e = \nu_e / Pr_e` for effective turbulent Prandtl number 
:math:`Pr_e`.
