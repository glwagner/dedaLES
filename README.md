# DedaLES

This python package provides classes and functions for performing oceanic process studies 
using Direct Numerical Simulation (DNS) and Large Eddy Simulation (LES) in [Dedalus].

## Basic usage

With [Dedalus] installed, navigate to `examples/` and do

```
mpiexec python basic_example.py
```

## LES Closures

Currently impelemented:

* Constant Smagorinsky ([wikipedia][wiki_const_smag], [Vreugdenhil and Taylor 2018])

Definitely planned:

* Anisotropic Minimum Dissipation ([Rozema et al. 2015], [Vreugdenhil and Taylor 2018])

Aspirational:

* Dynamic Smagorinsky ([wikipedia][wiki_dyn_smag], [Germano et al. 1991])
* Sullivan-McWilliams-Moeng ([Sullivan et al. 1994])

[Dedalus]: http://dedalus-project.org

[wiki_const_smag]: https://en.wikipedia.org/wiki/Large_eddy_simulation#Smagorinsky–Lilly_model
[Vreugdenhil and Taylor 2018]: https://aip.scitation.org/doi/pdf/10.1063/1.5037039?class=pdf
[wiki_dyn_smag]: https://en.wikipedia.org/wiki/Large_eddy_simulation#Germano_dynamic_model
[Germano et al. 1991]: https://aip.scitation.org/doi/abs/10.1063/1.857955
[Sullivan et al. 1994]: https://link.springer.com/article/10.1007/BF00713741
[Rozema et al. 2015]: https://aip.scitation.org/doi/abs/10.1063/1.4928700
