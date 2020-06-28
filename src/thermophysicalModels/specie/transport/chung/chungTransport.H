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
    Foam::chungTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    chungTransportI.H
    chungTransport.C

\*---------------------------------------------------------------------------*/

#ifndef chungTransport_H
#define chungTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class chungTransport;

template<class Thermo>
inline chungTransport<Thermo> operator+
(
    const chungTransport<Thermo>&,
    const chungTransport<Thermo>&
);

template<class Thermo>
inline chungTransport<Thermo> operator*
(
    const scalar a,
    const chungTransport<Thermo>& b
){FatalErrorInFunction
            << " chungTransport<Thermo> operator*"
            << exit(FatalError);return *(new chungTransport<Thermo>(b));};


template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const chungTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class chungTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class chungTransport
:
    public Thermo
{
    // Private data
    public:
        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number []
        scalar kappa_;


    // Private Member Functions

        //- Construct from components
        inline chungTransport
        (
            const Thermo& t
        );


public:

    // Constructors

        //- Construct as named copy
        inline chungTransport(const word&, const chungTransport&);

        //- Construct from dictionary
        chungTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<chungTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<chungTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "chung<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T,const scalarList *Xinp=NULL) const{return 0;};

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const{return 0;};

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const{return 0;};

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators




    // Friend operators




    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const chungTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chungTransportI.H"

#ifdef NoRepository
    #include "chungTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
