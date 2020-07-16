/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::HYreactingMixture

Description
    Foam::HYreactingMixture

SourceFiles
    HYreactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef HYreactingMixture_H
#define HYreactingMixture_H

#include "speciesTable.H"
#include "chemistryReader.H"
//#include "multiComponentMixture.H"
#include "HYmultiComponentMixture.H"//XYY
#include "thermotable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class HYreactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class HYreactingMixture
:
    public speciesTable,
    public autoPtr<chemistryReader<ThermoType>>,
   // public autoPtr<chemistryReader<thermotable<ThermoType>>>,
    public HYmultiComponentMixture<ThermoType>, //XYY
    public PtrList<Reaction<ThermoType>>
    //public PtrList<Reaction<thermotable<ThermoType>>>
{
    // Private member data

        //- Table of species composition
        speciesCompositionTable speciesComposition_;


public:

    //- The type of thermo package this mixture is instantiated for
    typedef  thermotableH<ThermoType> thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        HYreactingMixture(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        HYreactingMixture(const HYreactingMixture&) = delete;


    //- Destructor
    virtual ~HYreactingMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "HYreactingMixture<" + thermoType::typeName() + '>';
        }

        //- Read dictionary
        void read(const dictionary&);

        label size() const
        {
            return PtrList<Reaction<thermoType>>::size();
        }

        Reaction<thermoType>& operator[](const label i)
        {
            return PtrList<Reaction<thermoType>>::operator[](i);
        }

        const Reaction<thermoType>& operator[](const label i) const
        {
            return PtrList<Reaction<thermoType>>::operator[](i);
        }

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const HYreactingMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HYreactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
