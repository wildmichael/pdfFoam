# Particle Property Extraction

The vector of the extracted, time-averaged quantities is given by

![Extracted variable equations](http://quicklatex.com/cache3/f5/ql_190391f7ad38968c2e25c9ae322371f5_l3.png)

Here, ![Kernel
function](http://quicklatex.com/cache3/0e/ql_6bafcab67c452a8d96e8cab7c0de6a0e_l3.png)
is a kernel function (such as a top-hat function coinciding with the cell
volume) and &alpha; is the time-averaging parameter which should be close
to 1. The particle density ![Particle
density](http://quicklatex.com/cache3/5c/ql_054666487a29ec33a4439aac48c4b65c_l3.png)
is given by the equation of state (e.g. perfect gas law).

The extracted, time-averaged variables can then be computed as follows:

![Time averaged variable equations](http://quicklatex.com/cache3/fe/ql_30ea21c71120b67982de7fc71de762fe_l3.png)

Above described mean estimation procedure is implemented in
`Foam::mcParticleCloud::updateCloudPDF()`. The statistical, instantaneous
moments are represented by the `*MomInstant` variables and the time-averaged
moments by the `Foam::mcParticleCloud::*Mom_` member variables.
