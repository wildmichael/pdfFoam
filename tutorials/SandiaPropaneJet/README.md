# Sandia Propane Jet

## Experimental Data
The experimental data is available from the website of the
[TNF workshop](https://tnfworkshop.org/data-archives/simplejet/propanejet/).

## Case Description
The inflow profiles of the co-flow have been computed from the experimental
data. For the jet, however, no such data is available, and hence the profiles
from the Sandia Flame D have been scaled. The TKE in the co-flow was computed
as

![Turbulent kinetic energy in co-flow](https://latex.codecogs.com/svg.latex?k=\frac{3}{2}{u'}^2)

in the jet it was assumed to be

![Turbulent kinetic energy in jet](https://latex.codecogs.com/svg.latex?k=\frac{1}{2}\left({u'}^2+2{v'}^2\right%29\quad. )

The turbulent dissipation was computed making a turbulence in equilibrium
assumption, resulting in

![Turbulent dissipation](https://latex.codecogs.com/svg.latex?\varepsilon=\sqrt{C_\mu}k\left\|\frac{\partial\tilde{U}}{\partial&space;x}\right\|\quad.)
