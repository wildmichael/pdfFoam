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

\*---------------------------------------------------------------------------*/

#include "mcParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mcParticle::mcParticle
(
    const Cloud<mcParticle>& cloud,
    Istream& is,
    bool readFields
)
:
    Particle<mcParticle>(cloud, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            m_ = readScalar(is);
            is >> Updf_;
            is >> UParticle_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&m_),
                sizeof(m_) + sizeof(Updf_) + sizeof(UParticle_)
            );
        }
    }

    dt_ = 0.01;
    // Check state of Istream
    is.check("mcParticle::mcParticle(Istream&)");
}


void Foam::mcParticle::readFields(Cloud<mcParticle>& c)
{
    if (!c.size())
    {
        return;
    }
    Particle <mcParticle> :: readFields(c);

    IOField<scalar> m(c.fieldIOobject("m", IOobject::MUST_READ));
    c.checkFieldIOobject(c, m);

    IOField<vector> Updf(c.fieldIOobject("Updf", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Updf);

    IOField<vector> UParticle(c.fieldIOobject("UParticle", IOobject::MUST_READ));
    c.checkFieldIOobject(c, UParticle);

    label i = 0;
    forAllIter(Cloud<mcParticle>, c, iter)
    {
        mcParticle& p = iter();

        p.m_ = m[i];
        p.Updf_ = Updf[i];
        p.UParticle_ = UParticle[i];
        i++;
    }
}


void Foam::mcParticle::writeFields(const Cloud<mcParticle>& c)
{
    Particle<mcParticle>::writeFields(c);

    label np = c.size();

    IOField<scalar> m(c.fieldIOobject("m", IOobject::NO_READ), np);
    IOField<vector> UParticle(c.fieldIOobject("UParticle", IOobject::NO_READ), np);
    IOField<vector> Updf(c.fieldIOobject("Updf", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(Cloud<mcParticle>, c, iter)
    {
        const mcParticle& p = iter();

        m[i] = p.m_;
        Updf[i] = p.Updf_;
        UParticle[i] = p.UParticle_;
        i++;
    }

    m.write();
    Updf.write();
    UParticle.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const mcParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const Particle<mcParticle>&>(p)
            << token::SPACE << p.m_
            << token::SPACE << p.Updf_
            << token::SPACE << p.UParticle_;
    }
    else
    {
        os  << static_cast<const Particle<mcParticle>&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.m_),
            sizeof(p.m_) + sizeof(p.Updf_)  + sizeof(p.UParticle_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const mcParticle&)");

    return os;
}


// ************************************************************************* //
