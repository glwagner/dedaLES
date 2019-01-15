# Notes on numerical methods for LES applications

## Dynamic Smagorinsky

### [Germano et al 1991](https://aip.scitation.org/doi/pdf/10.1063/1.857955?class=pdf)

* Fourier-Chebyshev pseudospectral discretization

## Anisotropic Minimum Dissipation

### [Rozema et al 2015](https://aip.scitation.org/doi/pdf/10.1063/1.4928700?class=pdf)

* "Collocated method for compressible flow"
* "Second- and fourth-order accuracy"
* ?? No information about time-stepping scheme ??
* Compressible freely-decaying grid turbulence and Kelvin-Helmholtz instability
* Incompressible turbulent channel flow

### [Abkar et al 2016]

* Pseudospectral in horizontal (Poincaré constant `Cx = Cy = 1/sqrt(12)`)
* Second-order differences in vertical (Poincaré constant `Cz = 1/sqrt(3)`)
* Second-order Adams-Bashforth time-stepping

### [Abkar and Moin 2017](https://link.springer.com/article/10.1007/s10546-017-0288-4)

* Similar to [Abkar et al 2016]

### [Vreugdenhil and Taylor 2018](https://aip.scitation.org/doi/pdf/10.1063/1.5037039?class=pdf)

* Pseudospectral in horizontal
* Second-order finite volume in vertical
* 3rd-order RK for nonlinear terms, semi-implicit Crank-Nicholson for viscous/diffusive

[Abkar et al 2016]: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.1.041701
