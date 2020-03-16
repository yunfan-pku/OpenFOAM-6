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
    Foam::HYmultiComponentMixture

Description
    Foam::HYmultiComponentMixture

SourceFiles
    HYmultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef HYmultiComponentMixture_H
#define HYmultiComponentMixture_H

#include "basicSpecieMixture.H"
#include "HashPtrTable.H"
#include "thermotable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ThermoType>
class thermotable;

/*---------------------------------------------------------------------------*\
                    Class HYmultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class HYmultiComponentMixture
:
    public basicSpecieMixture
{

public:
       //- The type of thermodynamics this mixture is instantiated for
    typedef thermotable<ThermoType> thermoType;

private:
    // Private Data

        //- Species data
        PtrList<ThermoType> speciesData_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable thermoType mixture_;
        

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable thermoType mixtureVol_;


    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const thermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Construct as copy (not implemented)
        HYmultiComponentMixture(const HYmultiComponentMixture<thermoType>&);


public:

 


    // Constructors

        //- Construct from dictionary, specie names, thermo database,
        //  mesh and phase name
        HYmultiComponentMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const fvMesh&,
            const word&
        );

        //- Construct from dictionary, mesh and phase name
        HYmultiComponentMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~HYmultiComponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "HYmultiComponentMixture<" + thermoType::typeName() + '>';
        }

        const thermoType& cellMixture(const label celli) const;

        const thermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const thermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const thermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Return the raw specie thermodynamic data
        const PtrList<thermoType>& speciesData() const
        {
            return speciesData_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& getLocalThermo(const label speciei) const
        {
            return speciesData_[speciei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HYmultiComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
