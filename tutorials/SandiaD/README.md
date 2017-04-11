# Sandia Flame D

## Experimental Data
The experimental data for the scalars is available from the website of the
[TNF workshop](http://www.sandia.gov/TNF/DataArch/FlameD.html).
The flow measurements are available upon request from Prof. Andreas Dreizler
from TU Darmstad: `dreizler AT ekt.tu-darmstadt.de`

The data files for the available measurement sets are:

- [Dfav.dat](expdata/Dfav.dat) (scalar data)
- [Dvel.dat](expdata/Dvel.dat) (velocity data)

## Case Description
The inflow profiles have been computed from the file `icbbody.dat`. The TKE
was computed as

![Turbulent kinetic energy](http://quicklatex.com/cache3/24/ql_788261c63f46efaebb3b3db9590a5824_l3.png)

and the turbulent dissipation was computed making a turbulence in equilibrium
assumption, resulting in

![Turbulent dissipation](http://quicklatex.com/cache3/e5/ql_8a51ac652a13e828a0b58846677919e5_l3.png)
