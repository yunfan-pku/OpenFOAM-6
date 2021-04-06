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

#include "nonlinearMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType, template<class> class ThermoMixtureType>
const ThermoType& Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(thermoDict.subDict(species_[i]))
        );
    }

    return speciesData_[0];
}


template<class ThermoType, template<class> class ThermoMixtureType>
void Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0 * Y_[0]);

    for (label n = 1; n < Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType, template<class> class ThermoMixtureType>
Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::nonlinearMixture
(
    const dictionary& thermoDict,
    const speciesTable& specieNames,
    const HashPtrTable<ThermoType>& thermoData,
    const fvMesh& mesh,
    const word& phaseName
)
    :
    basicSpecieMixture(thermoDict, specieNames, mesh, phaseName),
    speciesData_(species_.size()),
    mixture_("mixture", speciesData_, specieNames, thermoDict)//,
    //mixtureVol_("volMixture", speciesData_,specieNames,thermoDict)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*thermoData[species_[i]])
        );
    }

    correctMassFractions();
}


template<class ThermoType, template<class> class ThermoMixtureType>
Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::nonlinearMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
    :
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(thermoDict))//,
    //mixtureVol_("volMixture", speciesData_[0])
{
    correctMassFractions();
}

/*
template<class ThermoType>
Foam::nonlinearMixture<ThermoType>::nonlinearMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const HashPtrTable<ThermoType>& thermoData,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture(thermoDict, specieNames, mesh, phaseName),
    speciesData_( thermoDict),
    mixture_( thermoDict.subDict("nonlinearMixture")),
    mixtureVol_( thermoDict.subDict("nonlinearMixture"))
{
    correctMassFractions();
}


template<class ThermoType>
Foam::nonlinearMixture<ThermoType>::nonlinearMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(thermoDict)),
    mixtureVol_("volMixture", speciesData_[0])
{
    correctMassFractions();
}
*/



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class ThermoType, template<class> class ThermoMixtureType>
inline void Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::newTimeStep()
{
    mixture_.newTimeStep();
}

template<class ThermoType, template<class> class ThermoMixtureType>
const ThermoMixtureType<ThermoType>& Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::cellMixture
(
    const label celli
) const
{
    scalarList Y(Y_.size());
    for (label n = 0; n < Y_.size(); n++)
    {
        Y[n] = Y_[n][celli];
    }
    mixture_.setY(Y);

    return mixture_;
}


template<class ThermoType, template<class> class ThermoMixtureType>
const ThermoMixtureType<ThermoType>& Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{

    scalarList Y(Y_.size());
    for (label n = 0; n < Y_.size(); n++)
    {
        Y[n] = Y_[n].boundaryField()[patchi][facei];
    }
    mixture_.setY(Y);


    return mixture_;
}


template<class ThermoType, template<class> class ThermoMixtureType>
const ThermoMixtureType<ThermoType>& Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalarList Y(Y_.size());
    for (label n = 0; n < Y_.size(); n++)
    {
        Y[n] = Y_[n][celli];
    }
    mixture_.setY(Y);

    return mixture_;
}


template<class ThermoType, template<class> class ThermoMixtureType>
const ThermoMixtureType<ThermoType>& Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalarList Y(Y_.size());
    for (label n = 0; n < Y_.size(); n++)
    {
        Y[n] = Y_[n].boundaryField()[patchi][facei];
    }
    mixture_.setY(Y);


    return mixture_;
}


template<class ThermoType, template<class> class ThermoMixtureType>
void Foam::nonlinearMixture<ThermoType, ThermoMixtureType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_[i] = ThermoType(thermoDict.subDict(species_[i]));
    }
}


// ************************************************************************* //
