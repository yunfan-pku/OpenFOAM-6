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

/*---------------------------------------------------------------------------*/
/* PengRobinson EOS                                                          */
/* Thermal solver                                                            */
/* Developed by   Ping Yi,  03/01/2019,   UMN                                */
/* mix means each phase mixture; ave means two-phase mixture                 */
/*---------------------------------------------------------------------------*/

#include "PengRobinsonMixture.H"
#include "PengRobinson.H"
#include "specie.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "IFstream.H"
#include "string.H"
#include "tmp.H"
using namespace Foam;
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Specie>
Foam::PengRobinsonMixture<Specie>::PengRobinsonMixture(
    const dictionary &dict)
    : Specie(dict),
      bico(this->N_),
      phaseFlag(1)
{
    forAll(this->X_, i)
    {
        bico[i].resize(this->N_);
    }

    speciesTable ts(dict.subDict("BinaryInteractionParameter").lookup("specielist"));
    scalarList bilist(dict.subDict("BinaryInteractionParameter").lookup("matrix"));
    int indexi, indexj, len = ts.size();
    forAll(this->X_, i)
    {
        indexi = ts[this->species()[i]];
        forAll(this->X_, j)
        {
            indexj = ts[this->species()[j]];
            bico[i][j] = bilist[indexi * len + indexj];
            //Info<<bico[i][j]<<" ";
        }
        //Info<<endl;
    }
}

template <class Specie>
inline Foam::autoPtr<Foam::PengRobinsonMixture<Specie>>
Foam::PengRobinsonMixture<Specie>::clone() const
{
    return autoPtr<PengRobinsonMixture<Specie>>(
        new PengRobinsonMixture<Specie>(*this));
}

template <class Specie>
inline Foam::autoPtr<Foam::PengRobinsonMixture<Specie>>
Foam::PengRobinsonMixture<Specie>::New(
    const dictionary &dict)
{
    return autoPtr<PengRobinsonMixture<Specie>>(
        new PengRobinsonMixture<Specie>(dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::W // kg/mol
    () const
{
    const PengRobinsonMixture<Specie> &This = *this;
    scalar Wmix = 0.0;
    forAll(this->X_, i)
    {
        Wmix += this->X_[i] * This[i].W(); //PY: Y() is mole number of component in mixture
    }
    return Wmix * 1.0e-03;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::rho // PY:phase mixture density
    (
        scalar p,
        scalar T) const
{                                                  //RR: 10e+03*m3 Pa K-1 mol-1
    return W() * p / (Z(p, T) * RR * 1.0e-03 * T); // kg/m3
}
template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::A(
    scalar p,
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<Specie> &This = *this;
    scalar Amix = 0.0;
    forAll(X, i)
    {
        forAll(X, j)
        {
            Amix += X[i] * X[j] * sqrt(This[i].A(p, T) * This[j].A(p, T)) * (1.0 - bico[i][j]);
        }
    }
    return Amix;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::B(
    scalar p,
    scalar T,
    const scalarList *Xinp) const
{
    const scalarList *Xp = &this->X_;
    if (Xinp != NULL)
        Xp = Xinp;
    const scalarList &X = *Xp;

    const PengRobinsonMixture<Specie> &This = *this;
    scalar Bmix = 0.0;
    forAll(X, i)
    {
        Bmix += X[i] * This[i].B(p, T);
    }
    return Bmix;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::cubicSolve(
    scalar a2,
    scalar a1,
    scalar a0,
    label flag)
{

    scalar Q = (3.0 * a1 - a2 * a2) / 9.0;
    scalar Rl = (9.0 * a2 * a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2) / 54.0;
    scalar Q3 = Q * Q * Q;
    scalar D = Q3 + Rl * Rl;
    scalar root = -1;
    if (D <= 0)
    {
        scalar th = ::acos(Rl / sqrt(-Q3));
        scalar qm = 2 * sqrt(-Q);
        scalar r1 = qm * cos(th / 3) - a2 / 3;
        scalar r2 =
            qm * cos((th + 2 * constant::mathematical::pi) / 3) - a2 / 3;
        scalar r3 =
            qm * cos((th + 4 * constant::mathematical::pi) / 3) - a2 / 3;

        // three root;
        if (flag == 1) //gas
        {
            root = max(r1, max(r2, r3));
        }
        else //liquid
        {
            root = max(r1, max(r2, r3));
            root = min(r1 > 1e-3 ? r1 : root, root);
            root = min(r2 > 1e-3 ? r2 : root, root);
            root = min(r3 > 1e-3 ? r3 : root, root);
        }
    }
    else
    {
        // One root is real: in supercritical zone
        scalar D05 = sqrt(D);
        scalar S = sign(Rl + D05) * pow(mag(Rl + D05), 1.0 / 3.0);
        scalar Tl = sign(Rl - D05) * pow(mag(Rl - D05), 1.0 / 3.0);
        /*
        scalar RlD = Rl + D05;
        scalar S = 0.0;
        if (RlD > 0.0)
        {
            S = sign(Rl + D05) pow(mag(Rl + D05), 1.0 / 3.0);
        }
        else
        {
            S = -pow(-RlD, 1.0 / 3.0);
        }
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0 / 3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0 / 3.0);
        }
        */
        root = S + Tl - a2 / 3;
    }

    if (root <= 0.0)
        FatalErrorInFunction << "Z is negative! Z = " << root << exit(FatalError);
    return root;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::Z(
    scalar p,
    scalar T,
    const label *pf,
    const scalarList *Xinp) const
{
    scalar A = this->A(p, T, Xinp), B = this->B(p, T, Xinp);
    scalar a2 = B - 1.0;
    scalar a1 = A - 2.0 * B - 3.0 * B * B;
    scalar a0 = -A * B + B * B + B * B * B;
    if (pf != NULL)
        return cubicSolve(a2, a1, a0, *pf);
    else
        return cubicSolve(a2, a1, a0, phaseFlag);
}

template <class Specie>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<Specie>::A2(
    scalar p, //inp
    scalar T  //inp
) const
{
    scalarList &Amix2 = *(new scalarList(this->N_));
    
    forAll(this->X_, i)
    {
        Amix2[i] = 0.0;
        forAll(this->X_, j)
        {
            Amix2[i] += this->X_[j] * sqrt((*this)[i].A(p, T) * (*this)[j].A(p, T)) * (1.0 - bico[i][j]);
        }
    }
    return autoPtr<scalarList>(&Amix2);
}

// PY: fugacity coefficient of species in each phase
template <class Specie>
inline Foam::autoPtr<scalarList> Foam::PengRobinsonMixture<Specie>::fugacityCoefficient(
    scalar p, //input
    scalar T) const
{
    scalarList &fugcoef = *(new scalarList(this->N_));
    scalar Amixz = this->A(p, T);
    scalar Bmixz = this->B(p, T);
    scalar Zmix = 0.0;
    autoPtr<scalarList> pAmix2(A2(p, T));
    const scalarList &Amix2 = pAmix2();
    const scalar small8 = 1.e-08; //PY: need to be noticed

    //this->ABZmix(p, T, comp_inp, Amixz, Bmixz, Amix2, Bz);

    Zmix = this->Z(p, T);

    label temp = 1;
    scalar Zmix0 = this->Z(p, T, &temp);

    scalar comptrifile = 0.7;
    scalarList Xs(this->X_.size());

    if (phaseFlag == 0)
    {
        while (Zmix < 0.0)
        {
            Xs[0] = comptrifile;
            for (int i = 1; i < this->X_.size(); i++)
                Xs[i] = (1.0 - comptrifile) / (Xs.size() - 1);
            Zmix = this->Z(p, T, &phaseFlag, &Xs);
            comptrifile = comptrifile + 0.01;
            if (comptrifile >= 1.0)
            {
                Zmix = Zmix0;
            }
        }
    }

    scalar cc0 = 2.0 * sqrt(2.0);
    scalar cc1 = 1.0 + sqrt(2.0);
    scalar cc2 = 1.0 - sqrt(2.0);

    scalar zfactor00 = Zmix;
    scalar ZMB = zfactor00 - Bmixz;

    if (ZMB < 0.0)
    {
        zfactor00 = Bmixz + small8;
        ZMB = small8;
    }
    forAll(this->X_, spid)
    {
        scalar b = (*this)[spid].B(p, T);
        fugcoef[spid] = exp(b / Bmixz * (zfactor00 - 1.0) - log(ZMB) - Amixz / (cc0 * Bmixz) * (2.0 * Amix2[spid] / Amixz - b / Bmixz) * log((zfactor00 + cc1 * Bmixz) / (zfactor00 + cc2 * Bmixz)));
    }
    return autoPtr<scalarList>(&fugcoef);
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::H(scalar p, scalar T) const
{
    return 0;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}

template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::S(
    scalar p,
    scalar T) const
{
    return 0;
}
template <class Specie>
inline Foam::scalar Foam::PengRobinsonMixture<Specie>::psi(
    scalar p,
    scalar T) const
{
    const scalar Z = this->Z(p, T);

    return 1.0 / (Z * this->R() * T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// ************************************************************************* //
