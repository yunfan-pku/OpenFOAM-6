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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Thermo>
inline Foam::chungTransportPH<Thermo>::chungTransportPH(
    const Thermo &t)
    : Thermo(t)
{
}

template <class Thermo>
inline Foam::chungTransportPH<Thermo>::chungTransportPH(
    const word &name,
    const chungTransportPH &ct)
    : Thermo(name, ct)
{
}


template <class Thermo>
inline Foam::autoPtr<Foam::chungTransportPH<Thermo>>
Foam::chungTransportPH<Thermo>::clone() const
{
    return autoPtr<chungTransportPH<Thermo>>(
        new chungTransportPH<Thermo>(*this));
}

template <class Thermo>
inline Foam::autoPtr<Foam::chungTransportPH<Thermo>>
Foam::chungTransportPH<Thermo>::New(
    const dictionary &dict)
{
    return autoPtr<chungTransportPH<Thermo>>(
        new chungTransportPH<Thermo>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
