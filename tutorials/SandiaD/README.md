# Sandia Flame D

## Experimental Data
The experimental data for the scalars is available from the website of the
[TNF workshop](https://tnfworkshop.org/data-archives/pilotedjet/ch4-air/).
The flow measurements are available upon request from Prof. Andreas Dreizler
from TU Darmstad: `dreizler AT ekt.tu-darmstadt.de`

The data files for the available measurement sets are:

- [Dfav.dat](expdata/Dfav.dat) (scalar data)
- [Dvel.dat](expdata/Dvel.dat) (velocity data)

## Case Description
The inflow profiles have been computed from the file `icbbody.dat`. The TKE
was computed as

![Turbulent kinetic energy](https://latex.codecogs.com/svg.latex?k=\frac{3}{2}{u'}^2)

and the turbulent dissipation was computed making a turbulence in equilibrium
assumption, resulting in

![Turbulent dissipation](https://latex.codecogs.com/svg.latex?\varepsilon=\sqrt{C_\mu}k\left\|\frac{\partial\tilde{U}}{\partial&space;x}\right\|\quad.)
