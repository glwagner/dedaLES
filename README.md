# dedaLES

This python package provides classes and functions for performing
Direct Numerical Simulation (DNS) and Large Eddy Simulation (LES) 
of fluid flows using the [dedalus] framework for spectrally solving 
partial differential equations.

This code (and `README`) are a work in progress.

## Installation

To use `dedaLES`, one must first [install dedalus][install dedalus].

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

* [constant Smagorinsky]

Definitely planned:

* [anisotropic minimum dissipation]
* [modified constant Smagorinsky]

Aspirational:

* Dynamic Smagorinsky ([wikipedia][wiki_dyn_smag], [Germano et al. 1991])
* Sullivan-McWilliams-Moeng ([Sullivan et al. 1994])

[constant Smagorinsky]: https://dedales.readthedocs.io/en/latest/closures/constantsmagorinsky.html
[modified constant Smagorinsky]: https://dedales.readthedocs.io/en/latest/closures/constantsmagorinsky.html
[anisotropic minimum dissipation]: https://dedales.readthedocs.io/en/latest/closures/amd.html

[wiki_dyn_smag]: https://en.wikipedia.org/wiki/Large_eddy_simulation#Germano_dynamic_model
[Germano et al. 1991]: https://aip.scitation.org/doi/abs/10.1063/1.857955
[Sullivan et al. 1994]: https://link.springer.com/article/10.1007/BF00713741

[install dedalus]: https://dedalus-project.readthedocs.io/en/latest/installation.html
[dedalus]: http://dedalus-project.org
