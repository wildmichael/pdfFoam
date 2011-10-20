/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::IOLostParticles

Description
    Keeps track of lost particles

SourceFiles
    IOLostParticlesI.H

\*---------------------------------------------------------------------------*/

#ifndef IOLostParticles_H
#define IOLostParticles_H

#include "Cloud.H"
#include "SLList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class IOLostParticles Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class IOLostParticles :
    public regIOobject
{
    // Private data

        const Cloud<ParticleType>& cloud_;

        SLList<Tuple2<scalar, string> > info_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        IOLostParticles(const IOLostParticles<ParticleType>&);

        //- Disallow default bitwise assignment
        void operator=(const IOLostParticles<ParticleType>&);


public:

    // Static data

        //- Runtime type name information. Derived from cloud type.
        virtual const word& type() const
        {
            return cloud_.type();
        }

    // Constructors

        //- Construct from cloud
        IOLostParticles(const Cloud<ParticleType>& c);


    // Member Functions

        //- Number of lost particles
        label size() const { return info_.size(); }

        //- Add a particle to the lost-particles list
        void add(const ParticleType& p);

        virtual void readData();

        virtual bool write() const;

        virtual bool writeData(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOLostParticlesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
