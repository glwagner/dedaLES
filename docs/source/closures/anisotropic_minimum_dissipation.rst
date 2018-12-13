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

The AMD model is an eddy viscosity model in that the relationship between subgrid 
stress and rate of strain is

.. math::

    F^\bu_{ij} = 2 \nu_e S_{ij} \c

while the relationship between subgrid tracer gradients and subgrid tracer flux is 

.. math::

    F^\theta = -\kappa_e \bnabla \theta \c

where :math:`\nu_e` and :math:`\kappa_e` are the eddy viscosity and eddy diffusivity.
The eddy viscosity and diffusivity are defined in terms of the 'predictor' eddy viscosity 
and diffusivities :math:`\nu_e^\dagger` and :math:`\kappa_e^\dagger` such that

.. math::

    \nu_e = \r{max} \left [ \nu_e^\dagger, 0 \right ] \c

and

.. math::

    \kappa_e = \r{max} \left [ \kappa_e^\dagger, 0 \right ] \c

Ensure that :math:`\nu_e` and :math:`\kappa_e` are always greater than 0.

The buoyancy-modified predictor eddy viscosity :math:`\nu_e^\dagger` 
is determined via the direction-dependent formula (`Abkar et al 2017`_) 

.. math::

    \nu_e^\dagger = - \frac{ \left ( \hat{\d}_k  u_i \right ) \left ( \hat{\d}_k  u_j \right )  S_{ij}
                                - \left ( \hat{\d}_k  w \right ) \hat{\d}_k  b}
                           {\left ( \d_{\ell}  u_m\right )^2} \c


where :math:`\hat{\d}_k` is the anisotropic scaled gradient operator

.. math::

    \hat{\d}_i = \sqrt{C_i} \Delta_i \d_i

for 'filter width' :math:`\Delta_i` and Poincaré constant :math:`C_i`. 
A key feature of AMD is inclusion of a direction-dependent, anisotropic filter
width :math:`\Delta_i`. Nominally, the filter width is some constant factor of the
grid spacing in the :math:`i^{\r{th}}` direction. In principle, the Poincaré constant 
:math:`C_i` depends on the discretization method used in each direction.
Because dedaLES is spectral in all directions, we use a single Poincaré constant
constant :math:`C_i = C`, so that :math:`\nu_e^\dagger` becomes

.. math::

    \nu_e^\dagger = - C^2 \frac{ \Delta_k^2 u_{i,k} u_{j,k} S_{ij} 
                        - \Delta_k^2 w_{,k} b_{,k}}{\r{tr}(\bnabla \bu)} \c

where the tracer, or first invariant of the tensor :math:`\bnabla \bu` is

.. math::

    \r{tr}(\bnabla \bu) = u_x^2 + v_x^2 + w_x^2 + u_y^2 + v_y^2 + w_y^2 + u_z^2 + v_z^2 + w_z^2

Explicitly, :math:`\Delta_k u_{i,k} u_{j,k} S_{ij}` is

.. math::

    \begin{split}
    \Delta_k u_{i,k} u_{j,k} S_{ij} &= 
    \,     \Delta_1^2 \left (u_x^2 S_{11} + v_x^2 S_{22} + w_x^2 S_{33} + 2 u_x v_x S_{12} + 2 u_x w_x S_{13} + 2 v_x w_x S_{23} \right ) \\
    \, & + \Delta_2^2 \left (u_y^2 S_{11} + v_y^2 S_{22} + w_y^2 S_{33} + 2 u_y v_y S_{12} + 2 u_y w_y S_{13} + 2 v_y w_y S_{23} \right ) \\
    \, & + \Delta_3^2 \left (u_z^2 S_{11} + v_z^2 S_{22} + w_z^2 S_{33} + 2 u_z v_z S_{12} + 2 u_z w_z S_{13} + 2 v_z w_z S_{23} \right ) \\ 
    \end{split}
       
while :math:`\Delta_k^2 w_{,k} b_{,k}` is

.. math::

    \Delta_k^2 w_{,k} b_{,k} = \Delta_1^2 w_x b_x + \Delta_2^2 w_y b_y + \Delta_3^2 w_z b_z \p

The predictor eddy diffusivity :math:`\kappa_e` for a quantity :math:`\theta` is

.. math::

    \kappa_e^\dagger = 
        - \frac{ C_k^2 \Delta_k^2 u_{i,k} \theta_{,k} \theta_{,i}}{ \theta_{,\ell}^2 } 

Note that :math:`\Delta_k^2 u_{i,k} \theta_{,k} \theta_{,i}` expands to

.. math::

    \Delta_k^2 u_{i,k} \theta_{,k} \theta_{,i} = \bnabla \theta \bcdot \big ( 
          \Delta_1^2 \theta_x \bu_x 
        + \Delta_2^2 \theta_y \bu_y 
        + \Delta_3^2 \theta_z \bu_z \big ) \p
        

References
----------

- `Rozema et al 2015`_
- `Abkar et al 2016`_
- `Abkar et al 2017`_
- `Vreugdenhil and Taylor 2018`_

.. _Rozema et al 2015: https://aip.scitation.org/doi/pdf/10.1063/1.4928700
.. _Abkar et al 2016: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.1.041701
.. _Abkar et al 2017: https://link.springer.com/article/10.1007/s10546-017-0288-4 
.. _Vreugdenhil and Taylor 2018: https://aip.scitation.org/doi/abs/10.1063/1.5037039
