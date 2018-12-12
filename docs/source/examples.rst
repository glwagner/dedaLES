.. math::

    \newcommand{\b}[1]{\boldsymbol{#1}}
    \newcommand{\r}[1]{\mathrm{#1}}
    \newcommand{\bz}{\b{z}}
    \newcommand{\bu}{\b{u}}
    \newcommand{\bcdot}{\b{\cdot}}
    \newcommand{\d}{\partial}
    \newcommand{\ee}{\r{e}}

    \newcommand{\p}{\, .}
    \newcommand{\c}{\, ,}


.. _dedaLES: https://github.com/glwagner/dedaLES
.. _convection_example.py: https://github.com/glwagner/dedaLES/examples/convection_example.py
.. _Taylor and Ferrari 2010: https://aslopubs.onlinelibrary.wiley.com/doi/epdf/10.4319/lo.2011.56.6.2293

Examples
========

A few example problems with `dedaLES`_.


Convection into a stratified fluid (`convection_example.py`_)
-------------------------------------------------------------

In this initial value problem, a constant flux of buoyancy out of the surface
drives cooling and convection into a stratified fluid. The initial buoyancy 
profile is given by

.. math::

    b(t=0) = N^2_{\infty} (z+h_0) H_d(-z-h_0) 

where :math:`h_0` is the initial mixed layer depth, the buoyancy gradient below
the mixed layer is :math:`\d_z b(z=-Lz) = N^2_{\infty}`, and :math:`H_d(\zeta)`
is a smoothed Heaviside function that goes from
:math:`0` to :math:`1` across :math:`\zeta=0` with transition width :math:`d`:

.. math::

    H_d(\zeta) = \tfrac{1}{2} \left ( 1 + \tanh \left [ \zeta / d \right ] \right ) \p

The lower boundary condition specifies a constant flux such that

.. math::

    \d_z b(z=-L_z) = N^2_{\infty} \p

The upper boundary condition specifies an asymptotically constant cooling rate. 
In terms of the buoyancy flux :math:`F = \kappa \d_z b`, the cooling rate 
in :math:`\r{W \, m^{-2}}` is

.. math::

    Q = \frac{c_P \rho_0 F}{\alpha g}

The buoyancy gradient at the surface associated with a given :math:`Q` is thus

.. math::

    \d_z b(z=0, t \to \infty) = b_{0z} = Q \frac{\alpha g}{c_P \rho_0 \kappa} \p

As in `Taylor and Ferrari 2010`_, we decrease the surface buoyancy flux from 
:math:`0` to its asymptotic value :math:`b_{0z}` gradually, such that
the boundary condition is

.. math::

    \d_z b(z=0, t) = b_{0z} \tanh(t/t_0) \c

where :math:`t_0` is the time-scale over which the buoyancy flux decreases.

We use the following parameters for our examples:

====================    ============================        ==================================
     Parameters                   Values                                  Units
====================    ============================        ==================================
:math:`L_x, L_y`        :math:`200`                         :math:`\r{m}`
:math:`L_z`             :math:`100`                         :math:`\r{m}`
:math:`n_x, n_y`        :math:`64`                          None  
:math:`n_z`             :math:`32`                          None
:math:`Q`               :math:`-0.1, -10, -1000`            :math:`\r{W \, m^{-2}}`
:math:`t_0`             :math:`1`                           :math:`\r{day}`
:math:`N^2_{\infty}`    :math:`9.5 \times 10^{-3}`          :math:`\r{s^{-2}}`
:math:`h_0`             :math:`50`                          :math:`\r{m}`
:math:`d`               :math:`10`                          :math:`\r{m}`
:math:`\alpha`          :math:`2.5 \times 10^{-4}`          :math:`\r{K^{-1}}`
:math:`g`               :math:`9.81`                        :math:`\r{m \, s^{-2}}`
:math:`\rho_0`          :math:`1028`                        :math:`\r{kg \, m^{-3}}`
:math:`c_P`             :math:`3993`                        :math:`\r{J \, kg^{-1} \, K^{-1}}`
:math:`\kappa`          :math:`1.43 \times 10^{-7}`         :math:`\r{m^2 \, s^{-1}}`
====================    ============================        ==================================


