# DedaLES

This python package provides classes and functions for performing oceanic process studies 
using Direct Numerical Simulation (DNS) and Large Eddy Simulation (LES) in [Dedalus]().

## Basic usage

With [Dedalus]() installed, navigate to `examples/` and do

```
mpiexec python basic_example.py
```

## LES Closures

Currently impelemented:

* Constant Smagorinsky ([1](), [2]())

Definitely planned:

* Anisotropic Minimum Dissipation ([3]())

Aspirational:

* Dynamic Smagorinsky ([4](), [5]())
* Sullivan-McWilliams-Moeng ([6]())

[Dedalus]: http://dedalus-project.org

[1]: https://en.wikipedia.org/wiki/Large_eddy_simulation#Smagorinsky–Lilly_model
[2]: https://aip.scitation.org/doi/pdf/10.1063/1.5037039?class=pdf
[3]: https://aip.scitation.org/doi/pdf/10.1063/1.5037039?class=pdf
[4]: https://en.wikipedia.org/wiki/Large_eddy_simulation#Germano_dynamic_model
[5]: https://aip.scitation.org/doi/abs/10.1063/1.857955
[6]: https://link.springer.com/article/10.1007/BF00713741
