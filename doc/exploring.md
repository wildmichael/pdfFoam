# Exploring the Code

On this page, the most important classes are mentioned which should serve as
a starting point when exploring the code base.

## mcThermo
Variable-density solvers interact with the Monte-Carlo particles through the
`Foam::mcThermo` class. Its main task is to provide the mean density to the
finite volume solver. The fact that this density is computed from an ensemble
of stochastic particles is hidden away by this class.

## mcParticleCloud
The `Foam::mcParticleCloud` class is essentially a list of particles with added
functionality. In the current implementation it orchestrates the time-stepping
scheme in the `Foam::mcParticleCloud::evolve()` function, invoking the various
models derived from `Foam::mcModel` at the mid-point. The
`Foam::mcParticleCloud::updateCloudPDF()` function computes the mean variables
from the particle ensemble.

## mcParticle
The Monte-Carlo particles are represented by the `Foam::mcParticle` class. It
is mostly a container for the particle variables, such as the weight,
velocity, local time-stepping parameter, the scalar vector etc. The latter
is represented by `Foam::mcParticle::Phi()` and is managed in a dynamic way,
allowing the user to specify the scalars contained in this vector in the
`constant/thermophysicalProperties::mcThermoCloudProperties::scalarFields`
list. This class is also responsible for the correct _serialisation_ of the
data, i.e. reading and writing the data either from or to a file or a stream.
The latter is used when transferring a particle from one processor to another
in parallel runs. The `Foam::mcParticle::move()` function is responsible for
performing the spatial particle tracking for a given time step.

## mcModel
The `Foam::mcModel` class is the abstract base class for all models used in the
JPDF algorithm. The function `Foam::mcModel::correct()` first invokes the
virtual function `Foam::mcModel::updateInternals()` which should perform any
updates of the internal state that are common to all particles, and then
calls the virtual function `Foam::mcModel::correct(const Foam::mcParticle&)` for
every particle in the cloud. This function should apply the model to the
particle.

## mcVelocityModel
The base class for all particle velocity evolution models is
`Foam::mcVelocityModel`.

## mcMixingModel
`Foam::mcMixingModel` is the base class for all molecular mixing models.

## mcReactionModel
`Foam::mcReactionModel` is the base class for all chemistry models.

## mcOmegaModel
`Foam::mcOmegaModel` is the base class for all models computing the inverse
time scale needed in the stochastic particle equations.

## mcPositionCorrection
All position/density correction models derive from `Foam::mcPositionCorrection`.

## mcLocalTimeStepping
`Foam::mcLocalTimeStepping` is the base class for all local time-stepping
models.

## mcBoundary
For all boundary conditions, `Foam::mcBoundary` is the base class. The function
`Foam::mcBoundary::correct()` is called twice every time-step, once at the
beginning of the time step, and after all particles have finished their
time step. In the first call, the argument `afterMove` is false, in the
second it is true. The function `Foam::mcBoundary::hitPatch()` is invoked when
a particle path intersects with the boundary surface.
