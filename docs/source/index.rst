.. dedaLES documentation master file, created by
   sphinx-quickstart on Mon Dec 10 21:56:40 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _Fluid models: models.rst
.. _Subgrid closures: closures.rst


`dedaLES`_: Large Eddy Simulation with `dedalus`_
=================================================

`dedaLES`_ provides solvers for Large Eddy Simulation (LES) using 
the `dedalus`_ framework for solving partial differential equations 
with spectral methods.

`dedaLES`_ solvers are implemented as Python classes that encapsulate a specific 
equation set for fluid dynamics --- some flavor of Navier-Stokes --- 
in a specific geometry. The turbulence closures are then provided as classes 
which can be added to any of the fluid model equation sets.
We also provide an API for setting boundary conditions, initial conditions,
and simulating specific types of problems clearly and conveniently.

.. _dedalus: http://dedalus-project.org
.. _dedaLES: https://github.com/glwagner/dedaLES

:ref:`Fluid models`
-------------------

We intend to implement solvers for 

* Rotating, stratified Boussinesq flow in a channel
* Compressible flow
* Homogeneous isotropic turbulence 

all with arbitrary tracer fields.


:ref:`Subgrid closures`
-----------------------

We intend to implement the subgrid turbulent closures:

* Constant Smagorinsky (`Constant Smagorinsky wiki`_)
* Modified Constant Smagorinsky for stratified flow (`Pressel et al 2015`_) 
* Anisotropic minimum dissipation (`Rozema et al 2015`_, `Abkar et al 2016`_, `Vreugdenhil and Taylor 2018`_)
* A second-order closure for boundary layer turbulence (`Sullivan et al 1994`_, `Moeng and Wyngaard 1988`_)

.. _Constant Smagorinsky wiki: https://en.wikipedia.org/wiki/Large_eddy_simulation#Smagorinsky–Lilly_model
.. _Moeng and Wyngaard 1988: https://journals.ametsoc.org/doi/pdf/10.1175/1520-0469%281988%29045%3C3573%3ASAOLES%3E2.0.CO%3B2
.. _Sullivan et al 1994: https://link.springer.com/article/10.1007/BF00713741
.. _Pressel et al 2015: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1002/2015MS000496
.. _Rozema et al 2015: https://aip.scitation.org/doi/abs/10.1063/1.4928700
.. _Vreugdenhil and Taylor 2018: https://aip.scitation.org/doi/pdf/10.1063/1.5037039?class=pdf
.. _Abkar et al 2016: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.1.041701


Contents
========

.. toctree::
    :maxdepth: 1

    models.rst
    summaryclosures.rst
    examples.rst
    testcases.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
