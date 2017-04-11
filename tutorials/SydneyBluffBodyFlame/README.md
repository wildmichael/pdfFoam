# Sydney Bluff Body Stabilised Flame

## Experimental Data
The experimental data is available from the
[University of Sydney](http://sydney.edu.au/engineering/aeromech/thermofluids/bluff.htm).

The flame labelled _HM1E_ has been simulated. The data files for the various
available measurement sets are:

- [b4f3-b-s1.dat](expdata/b4f3-b-s1.dat) (velocity data, set 1)
- [b4f3-b-s2.dat](expdata/b4f3-b-s1.dat) (velocity data, set 2)
- [b4f3b_Fav.dat](expdata/b4f3b_Fav.dat) (scalar data)

# Case Description
The inflow profiles have been computed from the file `icbbody.dat`. The TKE
was computed as

![Turbulent kinetic energy](http://quicklatex.com/cache3/24/ql_788261c63f46efaebb3b3db9590a5824_l3.png)

and the turbulent dissipation was computed making a turbulence in equilibrium
assumption, resulting in

![Turbulent dissipation](http://quicklatex.com/cache3/e5/ql_8a51ac652a13e828a0b58846677919e5_l3.png)
