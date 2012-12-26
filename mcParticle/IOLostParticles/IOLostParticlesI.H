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

#include "IOLostParticles.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::IOLostParticles<ParticleType>::IOLostParticles
(
    const Foam::Cloud<ParticleType>& c
)
:
    regIOobject
    (
        IOobject
        (
            "lostParticles",
            c.time().timeName(),
            c,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    cloud_(c),
    info_()
{
    if (headerOk())
    {
        readData();
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::IOLostParticles<ParticleType>::add(const ParticleType& p)
{
    OStringStream oss;
    oss << p;
    Tuple2<scalar, string> pInfo(cloud_.time().value(), oss.str());
    info_.append(pInfo);
}


template<class ParticleType>
bool Foam::IOLostParticles<ParticleType>::write() const
{
    if (info_.size())
    {
        return regIOobject::write();
    }
    else
    {
        return true;
    }
}


template<class ParticleType>
bool Foam::IOLostParticles<ParticleType>::writeData(Ostream& os) const
{
    typedef SLList<Tuple2<scalar, string> > lostParticlesList;
    os<< info_.size() << nl << token::BEGIN_LIST << nl;
    forAllConstIter(lostParticlesList, info_, iter)
    {
        os << iter() << nl;
    }
    os<< token::END_LIST << endl;
    return os.good();
}


template<class ParticleType>
void Foam::IOLostParticles<ParticleType>::readData()
{
    Istream& is = readStream("");

    token firstToken(is);

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        is.readBeginList("SLList<Tuple2<scalar,string>>");

        for (label i=0; i<s; i++)
        {
            Tuple2<scalar, string> pi;
            is >> pi;
            info_.append(pi);
        }

        // Read end of contents
        is.readEndList("SLList<Tuple2<scalar,string>>");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                "void IOLostParticles<ParticleType>::readData"
                "(SLList<Tuple2<scalar,string> >&, bool)",
                is
            )   << "incorrect first token, '(', found "
                << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);
            // Do not read any fields, position only
            Tuple2<scalar, string> pi;
            is >> pi;
            info_.append(pi);
            is >> lastToken;
        }
    }
    else
    {
        FatalIOErrorIn
        (
            "void IOLostParticles<ParticleType>::readData"
            "(SLList<Tuple2<scalar,string> >&, bool)",
            is
        )   << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    // Check state of IOstream
    is.check
    (
        "void IOLostParticles<ParticleType>::readData"
        "(SLList<Tuple2<scalar,string> >&, bool)"
    );
}

// ************************************************************************* //
