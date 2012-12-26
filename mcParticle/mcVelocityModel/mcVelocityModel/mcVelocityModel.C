/*---------------------------------------------------------------------------*\
                pdfFoam: General Purpose PDF Solution Algorithm
                   for Reactive Flow Simulations in OpenFOAM

 Copyright (C) 2012 Michael Wild, Heng Xiao, Patrick Jenny,
                    Institute of Fluid Dynamics, ETH Zurich
-------------------------------------------------------------------------------
License
    This file is part of pdfFoam.

    pdfFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) version 3 of the same License.

    pdfFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with pdfFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mcVelocityModel.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcVelocityModel, 0);
    defineRunTimeSelectionTable(mcVelocityModel, mcVelocityModel);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcVelocityModel::mcVelocityModel
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcModel(cloud, db, subDictName)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcVelocityModel> Foam::mcVelocityModel::New
(
    mcParticleCloud& cloud,
    const objectRegistry& db
)
{
    word velocityModelType(cloud.thermoDict().lookup("velocityModel"));
    word sd = velocityModelType + "VelocityModelCoeffs";

    mcVelocityModelConstructorTable::iterator cstrIter =
        mcVelocityModelConstructorTablePtr_->find(velocityModelType);

    if (cstrIter == mcVelocityModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcVelocityModel::New(const objectRegistry&, const dictionary&)"
        )   << "Unknown mcVelocityModel type " << velocityModelType
            << endl << endl
            << "Valid mcVelocityModel types are :" << endl
            << mcVelocityModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mcVelocityModel>(cstrIter()(cloud, db, sd));
}

// ************************************************************************* //
