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
    Foam::multispecie

Description
    Base class of the thermophysical property types.

SourceFiles
    multispecieI.H
    multispecie.C

\*---------------------------------------------------------------------------*/

#ifndef multispecie_H
#define multispecie_H

#include "word.H"
#include "scalar.H"
#include "scalarList.H"
#include "dictionary.H"
#include "speciesTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of friend functions and operators
    template <class Thermo>
    class multispecie;

    template <class Thermo>
    Ostream &operator<<(Ostream &, const multispecie<Thermo> &);

    /*---------------------------------------------------------------------------*\
                           Class multispecie Declaration
\*---------------------------------------------------------------------------*/

    template <class Thermo>
    class multispecie
    {
    public:
        // Private data

        //- No of specie
        label N_;

        //- mass fraction
        List<scalar> Y_;

        //- mole fraction
        List<scalar> X_;

        //mutable label changed;
        //- thermo lib
        PtrList<Thermo>& speciesData_;

        //- species name list
        const speciesTable& species_;

    public:
    typedef Thermo SingleThermoType;
        //- Runtime type information
        //ClassName("multispecie");
        static word typeName()
        {
            return  " ";//Thermo::typeName() ;
        }

        // Constructors

        //- Construct from components without name
        //inline multispecie(const scalar Y, const scalar molWeight);

        //- Construct from components with name
        /*
        inline multispecie(
            const word &name,
            const scalar Y,
            const scalar molWeight);
        */
        //- Construct as named copy
        //inline multispecie(const word &name, const multispecie &);

        //- Construct from dictionary
        multispecie(const dictionary &dict,PtrList<Thermo> &);

        inline multispecie(const word &name, PtrList<Thermo> &, const speciesTable &, const dictionary &);


        // Member Functions

        // Access

        //- Name
        //inline const word& name() const;

        //- Molecular weight [kg/kmol]
        //inline scalar W() const;

        //- No of moles of this multispecies in mixture
        inline const scalarList &Y() const;

        inline const scalarList &X() const;
        inline const speciesTable &species() const;
        //- Gas constant [J/(kg K)]
        //inline scalar R() const;


        inline void setX(const scalarList &);
        inline void setY(const scalarList &);
        inline const PtrList<Thermo>& speciesData() const;
        inline const Thermo &operator[](const label) const;
        // I-O

        //- Write to Ostream
        void write(Ostream &os) const;

        // Ostream Operator

        friend Ostream &operator<<<Thermo>(Ostream &, const multispecie &);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multispecieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
