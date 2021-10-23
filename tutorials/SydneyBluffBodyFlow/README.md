# Sydney Bluff Body Flow

## Experimental Data
The experimental data is available from the [University of
Sydney](https://web.aeromech.usyd.edu.au/thermofluids/bluff.php).

The flow field has been measured at _61m/s_ bulk jet velocity. The data files
for the three available measurement sets are:

- [b4c1-s1.dat](expdata/b4c1-s1.dat)
- [b4c1-s2.dat](expdata/b4c1-s2.dat)
- [b4c1-s3.dat](expdata/b4c1-s3.dat)


The scalar data is labelled `B4C1`, where three data sets are available:

- [b4c1a-01x.txt](expdata/b4c1a-01x.txt): _50m/s_ bulk jet velocity
- b4c1b-01x.txt: _63m/s_ bulk jet velocity
- b4c1c-01x.txt: _80m/s_ bulk jet velocity

## Case Description
Since the flow field and the scalar fields have been measured at different jet
bulk velocities, two different simulation cases are supplied.

### Flow Field
The case `flow` runs the standard `simpleFoam` solver to obtain a flow
field solution. No particles are involved here. The inflow profiles have been
computed from the file `icbbody.dat` for the jet, for the co-flow a
rectangular velocity profile is assumed as the boundary layer is very small
there. The RMS velocity fluctuations have been assumed to be roughly
_3.144%_ in the co-flow. A jet bulk velocity of _61m/s_ was found to
underpredict the jet velocity rather grossly, and it appears that _63m/s_
are more appropriate. The TKE was computed as

![Turbulent kinetic energy](https://latex.codecogs.com/svg.latex?k=\frac{3}{2}{u'}^2)

In the jet, the turbulent dissipation was computed making a turbulence in
equilibrium assumption, resulting in

![Turbulent dissipation](https://latex.codecogs.com/svg.latex?\varepsilon=\sqrt{C_\mu}k\left\|\frac{\partial\tilde{U}}{\partial&space;x}\right\|\quad.)

In the co-flow, the turbulent dissipation was chosen rather arbitrarily as
_&epsilon;=137m<sup>2</sup>/s<sup>3</sup>_.

### Scalar Fields
The case `scalarsInit` computes, similarly to the `flow` case described in
the previous section, the flow field using `simpleFoam`. The profiles in the
jet have been scaled for the bulk velocity of _50m/s_. This flow field is
then used in the case named `scalars` where the Monte-Carlo simulation is
performed with a frozen FV field, i.e. using one-way coupling.
