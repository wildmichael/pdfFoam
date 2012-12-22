XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
X                    The University of Sydney                           X
X                                                                       X
X	  Department of Mechanical and Mechatronic Engineering          X
X   >>                                                               << X  
X        >>                                                     <<      X
X             >>                   and                    <<            X
X        >>                                                     <<      X 
X   >>                                                               << X
X                     Sandia National Laboratories                      X
X                                                                       X
X                     Combustion Research Facility                      X 
X                                                                       X
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
X                                                                       X
X                Turbulent Nonpremixed Combustion Data                  X
X                                                                       X
X                                 for                                   X
X                                                                       X
X               Piloted and Bluff-Body Stabilised  Flows                X
X                                                                       X
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        IF YOU USE THE DATA PROVIDED HERE PLEASE E-MAIL YOUR NAME, 
      AFFILIATION AND E-MAIL ADDRESS TO A.R. MASRI. WE WILL THEN KEEP 
         YOU INFORMED OF ANY UPDATES OR ADDITIONS TO THE DATABASE	

A.R Masri
Department of Mechanical and Mechatronic Engineering 
The University of Sydney
NSW, 2006 Australia
Tel:   (61) 2 9351 2288
FAX:   (61) 2 9351 7060
Email: masri@mech.eng.usyd.edu.au
http://www.mech.eng.usyd.edu.au/research/energy/#data
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  DISCLAIMER

  No responsibility is assumed by the suppliers of these data for any 
  injury and/or property damage as a matter of products liability, 
  negligence or otherwise, or from any  use or operation of any 
  methods, products, instructions, or ideas based on these data.
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  REFERENCING

  To reference these data use:
  Combustion data base, The University of Sydney and 
  The Combustion Research Facility, Sandia National Laboratories, 
  http://www.mech.eng.usyd.edu.au/research/energy/#data
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                           5 September, 1997

                               Release 2.0

A comprehensive set of data collected in bluff-body stabilised flames and  
nonreacting jets is released. Measurements of flow, mixing, temperature as 
well as composition fields are presented. The composition field measurements 
include species such as CO, CO2, H2, H2O, O2, N2, Hydrocarbon, as well 
as OH and NO. These measurements complement data on pilot-stabilised 
flames which have been available on this web site since 1995. 

If you are interested in the pilot-stabilised flame data, you should 
read the file: Piloted

If you are interested in the bluff-body stabilised flame data, you 
should read the file: Bluffbd

==========================================================================
CNG and LPG fuels
-----------------
A volumetric composition is given for these fuel mixtures which are used 
extensively
CNG: 90.9% CH4, 5.0% C2H6, 1.1%C3H8, 2.4%CO2, balance:C4H10 and N2
LPG: 94.1% C3H8, 5.5% C3H6, 0.4% C4H10

Velocity Measurements
---------------------
Two colour LDV system with frequency shifted beams are used to measure the
horizontal and vertical velocity components. The fuel and the air are seeded
in order to reduce the seeding bias. Uncertainties of the LDV technique are 
mainly associated with the seeding bias due to steep temperature gradients  
and the presence of more than one particle in the probe volume. The error due 
to seeding bias is very hard to quantify and is believed to be small, however, 
the error due to the presence of more than one particle in the measurement 
volume is believed to be 4% for the mean and 7% for the rms fluctuations.

The flowfield data are collected at the University of Sydney and consist 
of radial profiles of mean and rms of fluctuations of the axial and radial 
velocity components at a range of axial locations.  The flowfield data are 
provided for selected jets and flames only. 

Temperature and composition Measurements
----------------------------------------
Temperature and composition data are instantaneous measurements collected 
at the Combustion Research Facility, Sandia National Laboratories, Livermore 
CA. Measurements have been made using the Raman/Rayleigh/LIF technique to 
give instantaneous and simultaneous temperature and the concentration of 
many species at a single point in the flame. The species measured are: 
N2, O2, CH4 (or CH3OH), CO, CO2, H2, H2O. Other species such as OH and NO 
are measured for selected flames only.  A range of fuel mixtures and flame 
velocities ranging from low to close to extinction have been studied. 

The following details give useful information about the processing and 
tabulation of the temperature and composition data:

1/ Each data file contains information about the axial and radial 
   measurement location.

2/ The mixture fraction is obtained using the Bilger formula 
   (Combust. Flame 80:135-149 (1990)) which is given by:

         2(Zc -Zc,o)   + (Zh-Zh,o)   - (Zo-Zo,o)
         -----------     ---------     ---------
              Wc            2Wh            Wo
   zi = ------------------------------------------
         2(Zc,f -Zc,o) + (Zh,f-Zh,o) - (Zo,f-Zo,o)
         -------------  -----------    -----------
              Wc            2Wh            Wo
  
   where Z(i) is a conserved scalar given by the total mass fraction of 
   element (i), and Wi is the molecular weight of elements (carbon, c, 
   hydrogen, h and oxygen, o). Subscripts (f) and (o) refer to the fuel 
   and air streams, respectively.

   Note: Values of Zo,o; Zh,o; Zc,o; Zo,f; Zh,f; and Zc,f used to 
         calculate zi are given for each fuel in the (*.dat) file 
         in the relevant directory.

3/ Mixture fraction ranges from zero to one. Negative values of mixture 
   fraction which may arise due to differential diffusion are not allowed. 
   Users interested in differential diffusion effects can redefine and 
   re-calculate their mixture fraction from the tabulated mass fractions.

4/ The argon contained in air is not accounted for and its mass fraction
   is lumped with that of nitrogen. 

5/ Each file contains Favre and ensemble mean and rms fluctuations for
   the data points contained in the file.

6/ Temperature is obtained either from the Rayleigh signal or from the sum 
   of the species number densities (assuming a mixture of ideal gases).

7/ The percentage mass fraction of species (Y(i)*100) is tabulated except 
   for NO where the percentage mass fraction*100 (Y(NO)*10000)) is given. 
   If the mass fraction of a species is zero for the entire data set in a 
   given directory it means that the species have not been measured. 
   Symbol h-c refers to the parent hydrocarbon fuel (CH4, CH3OH, etc...)

8/ The factor TNDR is defined as:
   TNDR = sum of species number densities measured from Raman and LIF over 
          the total number density obtained from the Rayleigh temperature, 
   or equivalently:
   TNDR = Temp. from Rayleigh / Temp. from sum of species number densities.

   TNDR = 0.0 implies that temperature is obtained from the sum of 
          species number densities.
   Otherwise temperature is obtained from Rayleigh.

   Generally, temperature is obtained from the Rayleigh measurements except 
   in cases where we suspect that the Rayleigh signal is corrupted by Mie  
   scattering. This is the case, for example, in flames of methanol where 
   there may be  scattering from fine droplets. Also, measuring in regions 
   of the flames where there are solid particles or soot particles corrupt 
   the Rayleigh signals.

9/ All the measured mass fractions are normalised such that the sum of the 
   tabulated mass fractions equals one. 

10/When the Rayleigh temperature is used (the factor TNDR is not equal to 
   zero) the original measured, non-normalised species mass fractions may 
   be recovered from the tabulated data as follows:

   Y(i, original)  = TNDR * Y(i, normalised, tabulated)
 
   where i corresponds to any of the tabulated species.
 
11/Departure of the total mass fraction of the measured species from unity 
   is only partly due to the fact that not all species existing in the probe 
   volume are being measured. There are random and systematic sources of 
   error on the measured signals leading to differences between the temperature 
   obtained from Rayleigh and that obtained from the sum of species number 
   densities.  A later section on Accuracy Considerations gives more details 
   about these errors. In cases where the TNDR factor varies significantly 
   from 1.0 the following guidelines are given:
   
   TNDR > 1.0
   ----------
   This generally implies that the species mass fractions are affected by 
   error which remains uncorrected for and which may be due to various 
   interference This leads to artificially high concentrations and hence a
   lower temperature from the sum of species number densities. An estimate 
   of this error on each species is extremely difficult to quantify but a 
   guide to the expected error is given later in the section on Accuracy 
   Considerations.

   TNDR < 1.0
   ----------
   This generally implies that the Rayleigh signal is subject to interference 
   due to Mie scattering leading to artificially lower Rayleigh temperatures 
   and hence a lower value of TNDR.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 A description of the burners, the initial and boundary conditions and 
 the experimental errors are given in the following files:
 Piloted: for the pilot-stabilsied jet flames
 Bluffbd: for the bluff-body stabilsied jets and flames
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Accuracy Considerations
-----------------------
There are a number of sources of error that have to be considered in evaluating 
the overall accuracy of laser-based, instantaneous measurements of species 
concentrations in flames. Only sources of error which may be influencing the 
data are mentioned here: (i)photon noise, (ii)interference error and 
(iii)spatial resolution. Photon noise is associated with the number of 
photons, n collected by a given detector at each laser pulse and it  decreases 
proportionally to 1/(n^0.5). This noise is expected to become significant at 
species mole fractions less than a few percent. The interference error depends 
on the magnitude of the fluorescence or chemiluminescence interference with 
the measured species. The error due to spatial resolution is not considered to 
be substantial and is discussed further below. Other sources of error which are
particular to the Raman system discussed here are (i)calibration drift due to 
the changes in laser lineshape over the lifetime of the dye, (ii)shot-to-shot 
variation in laser lineshape, and (iii)uncertainties in the temperature 
dependent bandwidth calibration factor f(T), especially for intermediate 
temperatures. 

The signal to noise ratio (S/N) gives a measure of the combined shot-to-shot 
random errors which are primarily due to photon statistics. It should be noted 
that, for scalar measurements in uniform steady flows and flames, the ratio of 
the rms fluctuations to the mean is given by the inverse of the S/N ratio. 

A measure of the S/N ratio can be obtained from the calibration data since 
these measurements are made in a uniform field of known temperature and species 
concentration. The reacting and nonreacting calibration data may be used here 
giving the variation of the S/N over a range of temperatures and species 
concentrations. Estimates of the S/N ratios obtained during experiments 
conducted between 1984 and 1992 are given in file stn92.dat. Since then the 
Raman-Rayleigh system has been improved by using better lasers and detectors. 
Estimates of the S/N ratios obtained during experiments conducted in 1995 are 
given in file stn95.dat. S/N ratios are given for typical data samples.  
The results are tabulated versus temperature for the Rayleigh signal and 
versus species number density (molecules/cm^3) for the Raman signals.  

S/N ratios for data collected 1992 and before    stn92.dat     (main directory)
S/N ratios for data collected 1995 and after     stn95.dat     (main directory)

The general trend for the signal to noise ratios presented in stn*.dat is to 
increase with number density. A correlation of the form S/N = A_i * [i]^0.5 
produces an adequate fit for all the scalars shown in stn*.dat. Here A_i 
is a constant and [i] is the number density of species. This square root 
dependence on the number density implies that photon statistics is a major 
contributor to noise on all of the Raman signals. 

Table 1 shows estimates of the percentage errors on various species for 
two typical samples collected in a CH4/H2 flame using the experimental setup.  
Lean and rich sample compositions are obtained from the actual data and are
taken here as illustrations of typical measurement conditions. It is evident 
that the improvements made in 1995 have led to a significant reduction in the 
percentage error on all scalars. The percentage error increases with 
decreasing number density or mole fraction. It should be emphasised that the 
errors reported here do not include the effect of interferences and spatial 
resolution. Raman interferences affect only selected species and are believed 
to have a small contribution to the overall error.  The fluorescence 
interference from soot precursors (mainly in the rich side of the flame) is 
very low in these flames and that improves the signal to noise ratio in all 
the affected Raman signals. Flames with high hydrocarbon fuels are most 
affected and among the Raman signals the CO line suffers the highest 
interference levels.

----------------------------------------------------------------------------
                            Table 1
Estimates of percentage error on typical samples of data 
collected collected in a turbulent methane-hydrogen flame.
----------------------------------------------------------------------------
Sample  Temperature  Species  % Mass    Number      %Error     %Error
                              Fraction  Density     1992       1995
                                                  or Before   or After
----------------------------------------------------------------------------
Lean      1900       CH4        0.0     0.0           -         -   
                     O2         4.0     0.12E18     17.0      10.0
                     N2        75.0     2.63E18      5.0       0.8
                     CO2        8.0     0.18E18     11.1       4.5
                     CO         2.0     0.07E18     16.6       9.0
                     H2         0.5     0.23E18     17.0      12.5
                     H2O       11.0     0.60E18      7.1       5.0
Rich      1400       CH4       18.0     1.09E18     10.0       2.3
                     O2         0.0     0.0           -         -
                     N2        57.0     1.98E18      6.3       1.1
                     CO2        5.5     0.12E18     12.0       5.5
                     CO         5.5     0.19E18     10.0       8.3
                     H2         2.5     1.22E18      6.9       4.0
                     H2O       12.0     0.65E18      7.2       4.0
----------------------------------------------------------------------------

Spatial Resolution
------------------
Spatial resolution issues for these data are discussed elsewhere [12]. The 
length of the measurement probe is 1mm and the diameter is about 0.6mm.  
Typical Kolmogorov length scales in the flame investigated range from 30 to 
150 microns. Using estimates of the length scales, it is found that the 
spatial resolution error ranges from 3% to 16% depending on the axial 
location in the flame and on the jet velocity. More information about 
spatial resolution effects may be found in Mansour, M.S. et al. (Combust. 
Flame 82:411 (1990)).

--------------------------------------------------------------------------
 Here is a list of references which may be used for further information. 
 Reference 12 reviews both the piloted and bluff-body stabilised flame 
 data and should be consulted first.
 Only part of the data availble on this web site are published in these 
 references

 1  Masri, A.R. and Bilger, R.W., `Turbulent Diffusion Flames of 
    Hydrocarbon Fuels Stabilised on a Bluff Body', 
    Twentieth Symposium (International) on Combustion. 
    The Combustion Institute, Pittsburgh, 1985, pp. 319-326.

 2  Dibble, R.W., Masri, A.R. and Bilger, R.W., `The Spontaneous Raman 
    Scattering Technique Applied to Nonpremixed Flames of Methane', 
    Combust. Flame 67:189-206 (1987).

 3  Masri, A.R. and Bilger, R.W., `Turbulent Nonpremixed Flames of 
    Hydrocarbon Fuels Near Extinction: Mean Structure from Probe Measurements', 
    Twenty-first Symposium (International) on Combustion, 
    The Combustion Institute, Pittsburgh, 1988, pp. 1511-5120.

 4  Masri, A.R., Dibble, R.W. and Bilger, R.W., `Turbulent Nonpremixed 
    Flames of Methane Near Extinction: Mean Structure from Raman Measurements', 
    Combust. Flame 71:245-266 (1988).

 5 Masri, A.R., Bilger, R.W. and Dibble, R.W., `Turbulent Nonpremixed 
   Flames of Methane Near Extinction:  Probability Density Functions', 
   Combust. Flame 73:261-285 (1988).

 6 Masri, A.R., Bilger, R.W. and Dibble, R.W., `Conditional Probability 
   Density Functions Measured in Turbulent Nonpremixed Flames of Methane 
   Near Extinction', Combust. Flame 74:267-284 (1988).

 7 Masri, A.R., Bilger, R.W. and Dibble, R.W., `The Local Structure 
   of Turbulent Nonpremixed Flames Near Extinction',  
   Combust. Flame 81:260-276 (1990).

 8 Masri, A.R., Dibble, R.W. and Barlow, R.S., `The Structure of Turbulent 
   Nonpremixed Flames of Methanol over a Range of Mixing Rates', 
   Combust. Flame 89:167-185 (1992).

 9 Masri, A.R., Dibble, R.W. and Barlow, R.S., `Chemical Kinetic Effects 
   in Nonpremixed Flames of $H_{2}/CO_{2}$ Fuel', 
   Combust. Flame 91:285-309 (1992).

10 Masri, A.R., Dibble, R.W. and Barlow, R.S., `Raman-Rayleigh Measurements 
   in Bluff Body Stabilised Flames of Hydrocarbon Fuels', 
   Twenty-fourth Symposium (International) on Combustion, 
   The Combustion Institute, Pittsburgh, 1992, pp. 317-324.

11 Masri, A.R., Dally, B.B., Barlow, R.S. and Carter, C.D., `The Structure 
   of The Recirculation Zone of a Bluff-Body Combustor', 
   Twenty-fifth Symposium (International) on Combustion, 
   The Combustion Institute, Pittsburgh, 1994, pp.1301-1308. 

12 Masri, A.R., Dibble, R.W., and Barlow, R.S., `The Structure of Turbulent 
   Nonpremixed Flames Revealed by Raman-Rayleigh-LIF Measurements', 
   Prog. Energy Combust. Sci., 22:307-362 (1997). 

13 Dally, B.B, Masri, A.R., Barlow, R.S., Fiechtner, G.J., and Fletcher, D.F., 
   `Measurements of NO in Turbulent Nonpremixed Flames Stabilised on a Bluff 
   Body', Twenty-sixth Symposium (International) on Combustion, 
   The Combustion Institute, Pittsburgh, 1996, Vol. 2, pp.2191-2197.


14 Dally, B.B, Masri, A.R., Barlow, R.S., Fiechtner, G.J., and Fletcher, D.F., 
   'Instantaneous and Mean Compositional Structure of Bluff-Body Stabilised 
   Nonpremixed Flames', Combust. Flame 114:119-148 (1998). 

15 Dally, B.B, D.F. Fletcher, and Masri, A.R., 
   'Flow and Mixing Fields of Turbulent Bluff-Body Jets and Flames', 
   Combustion Theory and Modeling 2:193-219 (1998).

