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

#include "ISATmultithermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

template<class ThermoMixture, template<class> class Type>
const Foam::scalar Foam::species::ISATmultithermo<ThermoMixture, Type>::tol_ = 1.0e-4;

template<class ThermoMixture, template<class> class Type>
const int Foam::species::ISATmultithermo<ThermoMixture, Type>::maxIter_ = 100;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoMixture, template<class> class Type>
Foam::species::ISATmultithermo<ThermoMixture, Type>::ISATmultithermo(const dictionary& dict,PtrList<ThermoMixture> &speciesData)
:
    ThermoMixture(dict,speciesData)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoMixture, template<class> class Type>
void Foam::species::ISATmultithermo<ThermoMixture, Type>::write(Ostream& os) const
{
    ThermoMixture::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class ThermoMixture, template<class> class Type>
Foam::Ostream& Foam::species::operator<<
(
    Ostream& os, const ISATmultithermo<ThermoMixture, Type>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
