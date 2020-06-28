/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "PengRobinson.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Specie>
Foam::PengRobinson<Specie>::PengRobinson(
    const dictionary &dict)
    : Specie(dict),
      Tc_(readScalar(dict.subDict("equationOfState").lookup("Tc"))),
      Vc_(readScalar(dict.subDict("equationOfState").lookup("Vc"))),
      Zc_(1.0),
      Pc_(readScalar(dict.subDict("equationOfState").lookup("Pc"))),
      omega_(readScalar(dict.subDict("equationOfState").lookup("omega"))),
      Hig_phase_(dict.subDict("equationOfState").lookup("Hig_phase")),
      Hig2_phase_(dict.subDict("equationOfState").lookup("Hig2_phase"))
{
    Zc_ = Pc_ * Vc_ / (RR * Tc_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Specie>
void Foam::PengRobinson<Specie>::write(Ostream &os) const
{
    Specie::write(os);
    dictionary dict("PengRobinson");
    dict.add("Tc", Tc_);
    dict.add("Vc", Vc_);
    dict.add("Zc", Zc_);
    dict.add("Pc", Pc_);
    dict.add("omega", omega_);
    dict.add("Hig_coef", Hig_phase_);
    dict.add("Hig2_coef", Hig2_phase_);

    os  << indent << dict.dictName() << dict;
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template <class Specie>
Foam::Ostream &Foam::operator<<(
    Ostream &os,
    const PengRobinson<Specie> &pg)
{
    pg.write(os);
    return os;
}

// ************************************************************************* //
