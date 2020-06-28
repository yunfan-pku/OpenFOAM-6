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

Class
    Foam::chungTransportPH

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    chungTransportPHI.H
    chungTransportPH.C

\*---------------------------------------------------------------------------*/

#ifndef chungTransportPH_H
#define chungTransportPH_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class chungTransportPH;
/*
template<class Thermo>
inline chungTransportPH<Thermo> operator+
(
    const chungTransportPH<Thermo>&,
    const chungTransportPH<Thermo>&
);

template<class Thermo>
inline chungTransportPH<Thermo> operator*
(
    const scalar,
    const chungTransportPH<Thermo>&
);
*/
template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const chungTransportPH<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class chungTransportPH Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class chungTransportPH
:
    public Thermo
{
    // Private data


    // Private Member Functions

        //- Construct from components
        inline chungTransportPH
        (
            const Thermo& t
        );


public:

    // Constructors

        //- Construct as named copy
        inline chungTransportPH(const word&, const chungTransportPH&);

        //- Construct from dictionary
        chungTransportPH(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<chungTransportPH> clone() const;

        // Selector from dictionary
        inline static autoPtr<chungTransportPH> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "chung<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]


        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream



    // Member operators


    // Friend operators






    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const chungTransportPH&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chungTransportPHI.H"

#ifdef NoRepository
    #include "chungTransportPH.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
