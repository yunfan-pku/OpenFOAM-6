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

#include "ISATVLEhePsiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEhePsiThermo<BasicPsiThermo, MixtureType>::calculate()
{
    //static scalar maxvlll = 0;
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
    scalarField& vaporfracCells = this->vaporfrac_.primitiveFieldRef();
    //const typename MixtureType::thermoType& mixture_ = this->cellMixture(celli);
    this->newTimeStep();
    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        //TCells[celli] = mixture_.THE(hCells[celli], pCells[celli], TCells[celli]);
        //psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        // psiCells[celli] = mixture_.psi_a(hCells[celli], pCells[celli], TCells[celli]);
        //psiCells[celli] = mixture_.psi_HP(hCells[celli], pCells[celli], TCells[celli]);
        //psiCells[celli] = mixture_.psi_a(hCells[celli], pCells[celli], TCells[celli]);
        // scalar psi_temp = mixture_.psi_HP(hCells[celli], pCells[celli], TCells[celli]);
       //  scalar Ttemp, psitemp;

        //std::tie(TCells[celli], psiCells[celli]) = mixture_.Tpsi_HP(hCells[celli], pCells[celli], TCells[celli]);
        std::tie(TCells[celli], psiCells[celli], vaporfracCells[celli]) = mixture_.Tpsivf_HP(hCells[celli], pCells[celli], TCells[celli]);

        /*
        if (hCells[celli] > 226106.08 && hCells[celli] < 226106.09 && pCells[celli]>10079795.33 && pCells[celli] < 10079795.34)
            maxvlll = maxvlll;
        psiCells[celli] = mixture_.psi_a(hCells[celli], pCells[celli], TCells[celli]);
        scalar psi_temp = mixture_.psi_HP(hCells[celli], pCells[celli], TCells[celli]);
        scalar psi_temp2 = mixture_.psi(pCells[celli], TCells[celli]);
        //psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        scalar psi_temp = mixture_.psi_HP(hCells[celli], pCells[celli], TCells[celli]);
        scalar psi_temp2 = mixture_.psi_a(pCells[celli], TCells[celli]);
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
 */
 // if (fabs((Ttemp - TCells[celli]) / TCells[celli]) > maxvlll)
 //      maxvlll = fabs((Ttemp - TCells[celli]) / TCells[celli]);
 //   Info << maxvlll << endl;



        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
        //vaporfracCells[celli] = mixture_.vaporfra(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    volScalarField::Boundary& vaporfracBf =
        this->vaporfrac_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];
        fvPatchScalarField& pvaporfrac = vaporfracBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvaporfrac[facei] = mixture_.vaporfra(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                //pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                //ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                //ppsi[facei] = mixture_.psi_HP(phe[facei], pp[facei], pT[facei]);
                //std::tie(pT[facei], ppsi[facei]) = mixture_.Tpsi_HP(phe[facei], pp[facei], pT[facei]);
                std::tie(pT[facei], ppsi[facei], pvaporfrac[facei]) = mixture_.Tpsivf_HP(phe[facei], pp[facei], pT[facei]);

                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                //pvaporfrac[facei] = mixture_.vaporfra(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::ISATVLEhePsiThermo<BasicPsiThermo, MixtureType>::ISATVLEhePsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
    :
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    vaporfrac_
    (
        IOobject
        (
            "vaporfrac",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    )
{
    calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::ISATVLEhePsiThermo<BasicPsiThermo, MixtureType>::~ISATVLEhePsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::ISATVLEhePsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info << "    Finished" << endl;
    }
}


// ************************************************************************* //
