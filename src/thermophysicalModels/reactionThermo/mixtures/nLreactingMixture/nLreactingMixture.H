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
    Foam::nLreactingMixture

Description
    Foam::nLreactingMixture

SourceFiles
    nLreactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef nLreactingMixture_H
#define nLreactingMixture_H

#include "speciesTable.H"
#include "chemistryReader.H"
#include "nonlinearMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class nLreactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType,template<class> class ThermoMixtureType>
class nLreactingMixture
:
    public speciesTable,
    public autoPtr<chemistryReader<ThermoType>>,
    public nonlinearMixture<ThermoType,ThermoMixtureType>,
    public PtrList<Reaction<ThermoType>>
{
    // Private member data

        //- Table of species composition
        speciesCompositionTable speciesComposition_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        nLreactingMixture(const nLreactingMixture&);

        //- Disallow default bitwise assignment
        void operator=(const nLreactingMixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoMixtureType<ThermoType> thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        nLreactingMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~nLreactingMixture()
    {}


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "nLreactingMixture<" + ThermoType::typeName() + '>';
        }

        //- Read dictionary
        void read(const dictionary&);

        label size() const
        {
            return PtrList<Reaction<ThermoType>>::size();
        }

        Reaction<ThermoType>& operator[](const label i)
        {
            return PtrList<Reaction<ThermoType>>::operator[](i);
        }

        const Reaction<ThermoType>& operator[](const label i) const
        {
            return PtrList<Reaction<ThermoType>>::operator[](i);
        }

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nLreactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
