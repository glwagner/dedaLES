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


.. _Constant Smagorinsky wiki: https://en.wikipedia.org/wiki/Large_eddy_simulation#Smagorinsky–Lilly_model
.. _Moeng and Wyngaard 1988: https://journals.ametsoc.org/doi/pdf/10.1175/1520-0469%281988%29045%3C3573%3ASAOLES%3E2.0.CO%3B2
.. _Sullivan et al 1994: https://link.springer.com/article/10.1007/BF00713741
.. _Pressel et al 2015: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1002/2015MS000496
.. _Rozema et al 2015: https://aip.scitation.org/doi/abs/10.1063/1.4928700
.. _Vreugdenhil and Taylor 2018: https://aip.scitation.org/doi/pdf/10.1063/1.5037039?class=pdf
.. _Abkar et al 2016: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.1.041701

.. _Subgrid closures:

Sub-grid models
===============

dedaLES impements a variety of models that approximate the impact of
unresolved, turbulent, 'subgrid' stress and tracer flux on the evolution 
of the resolved momentum and tracer fields. These models are often called
'subgrid closures'.

The models that dedaLES implements are 'eddy viscosity' and 'eddy diffusivity'
models, because the subgrid stress and tracer flux are assumed 
*proportional* to the resolved rate of strain and tracer gradients.
The constants of proportionality are the eddy viscosity and diffusivity.

The subgrid stress tensor :math:`F^{\bu}_{ij}` is thus written

.. math::

    F^{\bu}_{ij} = 2 \nu_e S_{ij} \c

where :math:`\nu_e` is the eddy viscosity, and

.. math::

    S_{ij} = \tfrac{1}{2} \left ( \d_i u_j + \d_j u_i \right ) \c

is the rate of strain tensor, often abbreviated as the 'strain tensor'.
The subgrid flux of a tracer :math:`\theta` is similarly

.. math::

    \b{F}^\theta = -\kappa_e \bnabla \theta \c

where :math:`\bnabla \theta` is the resolved tracer gradient and
:math:`\kappa_e` is the eddy diffusivity.
The differences between models lies entirely in how :math:`\nu_e` is 
calculated.

Implemented closures
--------------------

The following closure schemes for the subgrid-scale turbulent stress
and tracer flux are implemented in dedaLES:

.. toctree::
    :maxdepth: 2

    closures/constant_smagorinsky
    closures/modified_constant_smagorinsky
    closures/anisotropic_minimum_dissipation
