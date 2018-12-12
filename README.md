# dedaLES

This python package provides classes and functions for finding numerical 
solutions to the Navier-Stokes equations using 
Direct Numerical Simulation (DNS) and Large Eddy Simulation (LES) 
in [`dedalus`].

This code (and README) are a work in progress.

## Installation

To use `dedaLES`, one must first [install `dedalus`][install dedalus].

Then clone `dedaLES`:

```
git clone https://github.com/glwagner/dedaLES.git
```

## Basic usage

Navigate to `examples/` and do

```
mpiexec python basic_example.py
```

## Documentation

- [**latest**](https://dedales.readthedocs.io/en/latest/)

## Fluid models

* Boussinesq flow in a channel (with optional rotation)

## LES Closures

Currently implemented:

* Constant Smagorinsky ([wikipedia][wiki_const_smag], [Vreugdenhil and Taylor 2018])

Definitely planned:

* Anisotropic Minimum Dissipation ([Rozema et al. 2015], [Vreugdenhil and Taylor 2018])

Aspirational:

* Dynamic Smagorinsky ([wikipedia][wiki_dyn_smag], [Germano et al. 1991])
* Sullivan-McWilliams-Moeng ([Sullivan et al. 1994])


[wiki_const_smag]: https://en.wikipedia.org/wiki/Large_eddy_simulation#Smagorinsky–Lilly_model
[Vreugdenhil and Taylor 2018]: https://aip.scitation.org/doi/pdf/10.1063/1.5037039?class=pdf
[wiki_dyn_smag]: https://en.wikipedia.org/wiki/Large_eddy_simulation#Germano_dynamic_model
[Germano et al. 1991]: https://aip.scitation.org/doi/abs/10.1063/1.857955
[Sullivan et al. 1994]: https://link.springer.com/article/10.1007/BF00713741
[Rozema et al. 2015]: https://aip.scitation.org/doi/abs/10.1063/1.4928700
[install dedalus]: https://dedalus-project.readthedocs.io/en/latest/installation.html
[`dedalus`]: http://dedalus-project.org
