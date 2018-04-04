# Particle Property Extraction

The vector of the extracted, time-averaged quantities is given by

![Extracted variable equations](https://latex.codecogs.com/svg.latex?%5Cbegin%7Balign*%7D%20%5Cvec%7BE%7D%5E%7Bn&plus;1%7D%20%26%3D%20%5Calpha%20%5Cvec%7BE%7D%5En%20&plus;%20%281-%5Calpha%29%20%5Csum%5Chat%20g%28%5Cvec%7Bx%7D%5E*%29m%5E*%5Cvec%7Be%7D%5E*%20%5C%5C%20%26%3D%20%28M%2C%20V%2C%20%5Cvec%7BI%7D%2C%20K%2C%20Z%29%5ET%20%5Cquad%20%5Ctext%7Bwhere%7D%20%5C%5C%20%5Cvec%7Be%7D%20%26%3D%20%281%2C%201/%5Crho%5E*%2C%20%5Cvec%7BU%7D%5E*%2C%20%5Cfrac%7BU_i%5E*%20U_i%5E*%7D%7B2%7D%2C%20z%5E*%29%5ET%20%5Cquad%20.%20%5Cend%7Balign*%7D)

Here, ![Kernel
function](https://latex.codecogs.com/svg.latex?%5Cfn_cm%20%5Chat%20g%28%5Cvec%7Bx%7D%5E*%29)
is a kernel function (such as a top-hat function coinciding with the cell
volume) and &alpha; is the time-averaging parameter which should be close
to 1. The particle density ![Particle
density](https://latex.codecogs.com/svg.latex?\rho^*)
is given by the equation of state (e.g. perfect gas law).

The extracted, time-averaged variables can then be computed as follows:

![Time averaged variable equations](https://latex.codecogs.com/svg.latex?%5Cbegin%7Balign*%7D%20%5Cbar%5Crho%20%26%3D%20%5Cfrac%7BM%7D%7BV%7D%20%5C%5C%20%5Ctilde%7B%5Cvec%7BU%7D%7D%20%26%3D%20%5Cfrac%7B%5Cvec%7BI%7D%7D%7BM%7D%20%5C%5C%20%5Ctilde%20k%20%26%3D%20%5Cfrac%7BK%7D%7BM%7D-%5Cfrac%7B1%7D%7B2%7D%5Ctilde%20U_i%5Ctilde%20U_i%20%5C%5C%20%5Ctilde%20z%20%26%3D%20%5Cfrac%7BZ%7D%7BM%7D%20%5Cend%7Balign*%7D)

Above described mean estimation procedure is implemented in
`Foam::mcParticleCloud::updateCloudPDF()`. The statistical, instantaneous
moments are represented by the `*MomInstant` variables and the time-averaged
moments by the `Foam::mcParticleCloud::*Mom_` member variables.
