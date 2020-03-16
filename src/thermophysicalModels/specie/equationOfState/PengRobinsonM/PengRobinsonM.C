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

#include "PengRobinsonM.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::PengRobinsonM<Specie>::PengRobinsonM
(
    const speciesTable& specie,
    const dictionary& dict
)
:   single_specie(specie.size()),
    bico(specie.size()),
    Y_(specie.size())
{
    forAll(specie, i)
    {
        bico[i].resize(specie.size());

        single_specie.set
        (
            i,
            new PengRobinsonS<Specie>(dict.subDict(specie[i]))

        );
    }
    
    speciesTable ts(dict.subDict("BinaryInteractionParameter").lookup("specielist"));
    scalarList bilist(dict.subDict("BinaryInteractionParameter").lookup("matrix"));
    int indexi, indexj,len=ts.size();
    forAll(specie, i)
    {
        indexi=ts[specie[i]];
        forAll(specie, j)
        {
            indexj=ts[specie[j]];
            bico[i][j]=bilist[indexi*len+indexj];
            //Info<<bico[i][j]<<" ";
        }
        //Info<<endl;
    }


}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::PengRobinsonM<Specie>::setY(scalarList Yin)
{
    scalar sum=0;
    forAll(Y_,i)
    {
        sum+=Yin[i];
    }
    forAll(Y_,i)
    {
        Y_[i]=Yin[i]/sum;
    }

}
template<class Specie>
void Foam::PengRobinsonM<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PengRobinsonM<Specie>& pg
)
{
    pg.write(os);
    return os;
}


// ************************************************************************* //
