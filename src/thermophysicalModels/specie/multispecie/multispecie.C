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

#include "multispecie.H"
#include "constants.H"

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(multispecie, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class thermo>
Foam::multispecie<thermo>(const dictionary& thermoDict)
:speciesTable()
{
    wordList s(thermoDict.lookup("species"));
    this->transfer(s);

}

/*
Foam::multispecie::multispecie(const dictionary& dict)
:
    speciesTable()
    //name_(dict.dictName()),
    //Y_(dict.subDict("multispecie").lookupOrDefault("massFraction", 1.0)),
    //molWeight_(readScalar(dict.subDict("multispecie").lookup("molWeight")))
{

    wordList s(dict.lookup("species"));
    species.transfer(s);
}
*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multispecie::write(Ostream& os) const
{
    dictionary dict("multispecie");
    if (Y_ != 1)
    {
        dict.add("massFraction", Y_);
    }
    dict.add("molWeight", molWeight_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const multispecie& st)
{
    st.write(os);
    os.check("Ostream& operator<<(Ostream& os, const multispecie& st)");
    return os;
}


// ************************************************************************* //
