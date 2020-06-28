/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "chungTransportMixture.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoMixture>
Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture(const dictionary& dict,PtrList<SingleThermoType> &speciesData)
:
    ThermoMixture(dict,speciesData)
{}

template <class ThermoMixture>
inline Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture(
    const word &name,
    PtrList<SingleThermoType> &speciesData,
    const speciesTable &specieNames,
    const dictionary &thermoDict)
    : ThermoMixture(name, speciesData, specieNames, thermoDict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoMixture>
void Foam::chungTransportMixture<ThermoMixture>::chungTransportMixture::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    ThermoMixture::write(os);

    dictionary dict("transport");
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ThermoMixture>
Foam::Ostream& Foam::operator<<(Ostream& os, const chungTransportMixture<ThermoMixture>& ct)
{
    ct.write(os);
    return os;
}


// ************************************************************************* //
