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

#include "mcOmegaModel.H"
#include "mcParticleCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcOmegaModel, 0);
    defineRunTimeSelectionTable(mcOmegaModel, mcOmegaModel);

} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcOmegaModel::mcOmegaModel
(
    mcParticleCloud& cloud,
    const objectRegistry& db,
    const word& subDictName
)
:
    mcModel(cloud, db, subDictName)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mcOmegaModel> Foam::mcOmegaModel::New
(
    mcParticleCloud& cloud,
    const objectRegistry& db
)
{
    word omegaType(cloud.thermoDict().lookup("OmegaModel"));
    word sd = omegaType + "OmegaModelCoeffs";

    mcOmegaModelConstructorTable::iterator cstrIter =
        mcOmegaModelConstructorTablePtr_->find(omegaType);

    if (cstrIter == mcOmegaModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mcOmegaModel::New(const fvMesh&, const dictionary&)"
        )   << "Unknown mcOmegaModel type " << omegaType << endl << endl
            << "Valid mcOmegaModel types are :" << endl
            << mcOmegaModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mcOmegaModel>(cstrIter()(cloud, db, sd));
}

// ************************************************************************* //
