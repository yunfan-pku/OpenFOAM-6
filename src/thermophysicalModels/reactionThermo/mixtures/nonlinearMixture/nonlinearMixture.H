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
    Foam::nonlinearMixture

Description
    Foam::nonlinearMixture

SourceFiles
    nonlinearMixture.C

\*---------------------------------------------------------------------------*/

#ifndef nonlinearMixture_H
#define nonlinearMixture_H

#include "basicSpecieMixture.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class nonlinearMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType,template<class> class ThermoMixtureType>
class nonlinearMixture
:
    public basicSpecieMixture
{
    // Private data

        //- Species data
        PtrList<ThermoType> speciesData_;


        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoMixtureType<ThermoType> mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        //mutable ThermoMixtureType<ThermoType> mixtureVol_;


    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Construct as copy (not implemented)
        nonlinearMixture(const nonlinearMixture<ThermoType,ThermoMixtureType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoMixtureType<ThermoType> thermoType;


    // Constructors

        //- Construct from dictionary, specie names, thermo database,
        //  mesh and phase name
        nonlinearMixture
        (
            const dictionary&,
            const speciesTable& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const fvMesh&,
            const word&
        );

        //- Construct from dictionary, mesh and phase name
        nonlinearMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~nonlinearMixture()
    {}


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "nonlinearMixture<" + ThermoType::typeName() + '>';
        }

        const ThermoMixtureType<ThermoType>& cellMixture(const label celli) const;

        const ThermoMixtureType<ThermoType>& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const ThermoMixtureType<ThermoType>& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        inline void newTimeStep();

        const ThermoMixtureType<ThermoType>& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Return the raw specie thermodynamic data
        const PtrList<ThermoType>& speciesData() const
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
    #include "nonlinearMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
