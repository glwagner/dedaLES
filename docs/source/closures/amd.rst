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
models unresolved turbulent subgrid stress :math:`F^\bu_{ij}` as the product of 
the resolved rate of strain  :math:`S_{ij} = \tfrac{1}{2} \left ( \d_i u_j + \d_j u_i \right )` 
and an eddy viscosity :math:`\nu_e`. The subgrid tracer flux, :math:`F^b`, is
modeled similarly in terms of the resolved tracer gradient :math:`\d_i \theta` 
and an eddy viscosity :math:`\kappa_e`. Unlike :ref:`constant Smagorinsky`, however, 
which introduces a turbulent Prandtl number (typically :math:`O(1)`) to relate 
eddy viscosity to eddy diffusivity, the AMD eddy diffusivity :math:`\kappa_e` 
is determined by an analog of the method used to determine eddy visosity.

The relationship between subgrid stress and rate of strain in the AMD model is 

.. math::

    F^\bu_{ij} = 2 \nu_e S_{ij} \, ,

The eddy viscosity :math:`\nu_e` is determined by the direction-dependent formula

.. math::

    \nu_e^\dagger = - \frac{ C_k^2 \delta_k^2 u_{i,k} u_{j,k} S_{ij}}{ u_{m, \ell}^2 } \c

where :math:`C_i` and :math:`\delta_i` are the Poincaré constant and
grid spacing in the :math:`i^{\r{th}}` direction. In long form, this becomes

.. math::

    \begin{equation}
    \begin{split}
    \, & \nu_e^\dagger = \big [ \\
    \,& \quad C_1 \delta_1 \left (u_x^2 S_{11} + v_x^2 S_{22} + w_x^2 S_{33} + 2 u_x v_x S_{12} + 2 u_x w_x S_{13} + 2 v_x w_x S_{23} \right ) +\\
    \,& \quad C_2 \delta_2 \left (u_y^2 S_{11} + v_y^2 S_{22} + w_y^2 S_{33} + 2 u_y v_y S_{12} + 2 u_y w_y S_{13} + 2 v_y w_y S_{23} \right ) +\\
    \,& \quad C_3 \delta_3 \left (u_z^2 S_{11} + v_z^2 S_{22} + w_z^2 S_{33} + 2 u_z v_z S_{12} + 2 u_z w_z S_{13} + 2 v_z w_z S_{23} \right ) \\ 
    \,& \big ] / \left ( u_x^2 + u_y^2 + u_z^2 + v_x^2 + v_y^2 + v_z^2 + w_x^2 + w_y^2 + w_z^2 \right ) \p
    \end{split}
    \end{equation}
       

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


References
==========

- `Rozema et al 2015`_
- `Abkar et al 2016`_
- `Vreugdenhil and Taylor 2018`_

.. _Abkar et al 2016: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.1.041701
.. _Rozema et al 2015: https://aip.scitation.org/doi/pdf/10.1063/1.4928700
.. _Vreugdenhil and Taylor 2018: https://aip.scitation.org/doi/abs/10.1063/1.5037039
